"""The front end (:mod:`t1sim.frontend`): what t1sim adds to the segmentation.

The segmentation itself is CAT-Surface's and T1Prep's and is tested there.
Here: the decomposition of the normalised image into tissue fractions, the
tissue intensities read from the mixture of Gaussians, and how the three
steps are chained.
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import frontend  # noqa: E402


class TestTheFractions(unittest.TestCase):
    def fractions(self, *values):
        return frontend.fractions_from_label(np.array(values, dtype=np.float32).reshape(-1, 1, 1))[:, 0, 0]

    def test_pure_tissues(self):
        np.testing.assert_array_equal(self.fractions(1, 2, 3), np.eye(3, dtype=np.float32))

    def test_a_mixture_is_linear_between_two_tissues(self):
        np.testing.assert_allclose(self.fractions(1.25, 2.5), [[.75, .25, 0], [0, .5, .5]])

    def test_darker_than_csf_is_csf_and_brighter_than_wm_is_wm(self):
        np.testing.assert_allclose(self.fractions(0.3, 3.6), [[1, 0, 0], [0, 0, 1]], atol=1e-6)

    def test_nothing_and_far_too_bright_is_no_brain(self):
        np.testing.assert_array_equal(self.fractions(0, 4, 5), np.zeros((3, 3), np.float32))

    def test_they_sum_to_one_inside_the_brain(self):
        label = np.linspace(0.01, 3.99, 500, dtype=np.float32).reshape(-1, 1, 1)
        np.testing.assert_allclose(frontend.fractions_from_label(label).sum(-1), 1, atol=1e-6)

    def test_the_label_comes_back_between_csf_and_wm(self):
        label = np.linspace(1, 3, 201, dtype=np.float32).reshape(-1, 1, 1)
        back = frontend.label_from_fractions(frontend.fractions_from_label(label))
        np.testing.assert_allclose(back, label, atol=1e-6)


class TestTheTissueMeans(unittest.TestCase):
    params = {
        "classes": ["GM", "WM", "CSF", "bone", "soft tissue", "air"],
        "lkp": [0, 1, 1, 2, 2, 3, 4, 5],
        "mg": [1.0, 0.75, 0.25, 0.6, 0.4, 1.0, 1.0, 1.0],
        "mn": [100.0, 160.0, 200.0, 30.0, 90.0, 10.0, 120.0, 1.0],
    }

    def test_gm_and_wm_are_weighted_by_the_mixing_proportions(self):
        means = frontend.tissue_means(self.params)
        self.assertEqual(means["GM"], 100.0)
        self.assertAlmostEqual(means["WM"], 0.75 * 160 + 0.25 * 200)

    def test_csf_is_its_darkest_gaussian(self):
        self.assertEqual(frontend.tissue_means(self.params)["CSF"], 30.0)

    def test_a_class_without_a_gaussian_is_an_error(self):
        params = dict(self.params, lkp=[0, 1, 1, 3, 3, 3, 4, 5])
        with self.assertRaisesRegex(ValueError, "CSF"):
            frontend.tissue_means(params)


class TestTheChain(unittest.TestCase):
    """SANLM, the segmentation and LAS, with all three replaced by stand-ins."""

    shape = (6, 7, 8)

    def run_segment(self, scan):
        grid = nib.Nifti1Image(np.zeros(self.shape, np.float32), np.eye(4))
        mask = np.zeros(self.shape, bool)
        mask[1:5, 1:6, 1:7] = True
        first_label = np.where(mask, 2.0, 0.0).astype(np.float32)
        las = np.where(mask, 2.5 / 3, 0.9).astype(np.float32)   # 0.9 outside: must not survive
        seg = {
            "prob": "posteriors", "label": first_label,
            "corrected": np.full(self.shape, 7.0, np.float32), "grid": grid,
            "affine": np.eye(4), "world": np.eye(4), "native_mat": np.eye(4),
            "native_shape": scan.shape[:3], "params": dict(TestTheTissueMeans.params),
        }
        cat_surf = mock.Mock()
        cat_surf.vol_sanlm.side_effect = lambda data: data + 1
        unified_seg = mock.Mock()
        unified_seg.segment.return_value = seg
        unified_seg.brain_mask.return_value = mask
        unified_seg.working_grid_maps.return_value = ("affine", "brain", "first label")
        apply_las = mock.Mock(return_value=nib.Nifti1Image(las, np.eye(4)))
        with mock.patch.object(frontend, "_t1prep", return_value=(cat_surf, apply_las, unified_seg)):
            out = frontend.segment(scan, template_dir="/templates", space="template")
        return out, mask, cat_surf, apply_las, unified_seg

    def test_the_label_is_the_normalised_image_inside_the_mask(self):
        scan = nib.Nifti1Image(np.ones(self.shape, np.float32), np.eye(4))
        out, mask, _, apply_las, _ = self.run_segment(scan)
        np.testing.assert_allclose(out["label"][mask], 2.5, atol=1e-6)
        np.testing.assert_array_equal(out["label"][~mask], 0)
        np.testing.assert_allclose(out["fractions"][mask], [[0, .5, .5]] * mask.sum(), atol=1e-6)
        self.assertEqual(out["means"]["CSF"], 30.0)

    def test_las_gets_the_gm_mean_of_the_mixture_and_estimates_the_others(self):
        scan = nib.Nifti1Image(np.ones(self.shape, np.float32), np.eye(4))
        out, _, _, apply_las, _ = self.run_segment(scan)
        apply_las.assert_called_once_with("brain", "first label", verbose=False,
                                          peaks=(None, out["means"]["GM"], None))
        self.assertEqual(out["means"]["GM"], 100.0)

    def test_the_segmentation_gets_the_denoised_and_the_raw_scan_in_ras_order(self):
        data = np.arange(np.prod(self.shape), dtype=np.float32).reshape(self.shape)
        flipped = nib.Nifti1Image(data, np.diag([-1.0, 1, 1, 1]))          # LAS storage
        _, _, cat_surf, _, unified_seg = self.run_segment(flipped)
        denoised, raw = unified_seg.segment.call_args.args
        self.assertEqual(nib.aff2axcodes(raw.affine), ("R", "A", "S"))
        np.testing.assert_array_equal(raw.get_fdata(), data[::-1])
        np.testing.assert_array_equal(denoised.get_fdata(), data[::-1] + 1)
        self.assertEqual(unified_seg.segment.call_args.kwargs["template_dir"], "/templates")

    def test_the_denoised_scan_comes_back_on_the_grid_it_is_stored_on(self):
        data = np.arange(np.prod(self.shape), dtype=np.float32).reshape(self.shape)
        stored = nib.Nifti1Image(data, np.diag([-0.8, 0.9, 1.1, 1.0]))     # LAS storage, as many scans
        out, _, _, _, _ = self.run_segment(stored)
        self.assertEqual(out["scan"].shape, stored.shape)
        np.testing.assert_allclose(out["scan"].affine, stored.affine)
        np.testing.assert_array_equal(np.asarray(out["scan"].dataobj), data + 1)   # the stand-in filter adds 1
        # the same for axes in another order
        swapped = nib.Nifti1Image(np.transpose(data, (1, 2, 0)),
                                  np.array([[0, 0, 0.8, 0], [0.9, 0, 0, 0], [0, 1.1, 0, 0], [0, 0, 0, 1.0]]))
        out, _, _, _, _ = self.run_segment(swapped)
        self.assertEqual(out["scan"].shape, swapped.shape)
        np.testing.assert_array_equal(np.asarray(out["scan"].dataobj), np.asarray(swapped.dataobj) + 1)

    def test_the_deformation_is_handed_on_when_asked_for(self):
        scan = nib.Nifti1Image(np.ones(self.shape, np.float32), np.eye(4))
        grid = nib.Nifti1Image(np.zeros(self.shape, np.float32), np.eye(4))
        mask = np.ones(self.shape, bool)
        calls = {}

        def segment(denoised, raw, template_dir=None, verbose=False, deformation=False):
            calls["deformation"] = deformation
            return {"prob": None, "label": None, "corrected": np.ones(self.shape, np.float32), "grid": grid,
                    "affine": np.eye(4), "world": np.eye(4), "native_mat": np.eye(4), "native_shape": self.shape,
                    "params": dict(TestTheTissueMeans.params), "iy": "iy", "template_mat": "mat"}

        unified_seg = mock.Mock(segment=segment)
        unified_seg.brain_mask.return_value = mask
        unified_seg.working_grid_maps.return_value = (None, None, None)
        apply_las = mock.Mock(return_value=nib.Nifti1Image(np.full(self.shape, 2 / 3, np.float32), np.eye(4)))
        cat_surf = mock.Mock()
        cat_surf.vol_sanlm.side_effect = lambda data: data
        with mock.patch.object(frontend, "_t1prep", return_value=(cat_surf, apply_las, unified_seg)):
            plain = frontend.segment(scan, template_dir="/t", space="template")
            self.assertFalse(calls["deformation"])
            self.assertNotIn("iy", plain)
            out = frontend.segment(scan, template_dir="/t", deformation=True, space="template")
        self.assertTrue(calls["deformation"])
        self.assertEqual((out["iy"], out["template_mat"]), ("iy", "mat"))

    def test_a_brain_voxel_that_las_maps_to_zero_is_csf_and_no_hole(self):
        scan = nib.Nifti1Image(np.ones(self.shape, np.float32), np.eye(4))
        grid = nib.Nifti1Image(np.zeros(self.shape, np.float32), np.eye(4))
        mask = np.ones(self.shape, bool)
        mask[0] = False                                         # outside the brain
        las = np.full(self.shape, 2 / 3, np.float32)
        dark = (1, 1, 1)
        las[dark] = 0.0                                         # darker than everything LAS knows
        unified_seg = mock.Mock()
        unified_seg.segment.return_value = {
            "prob": None, "label": None, "corrected": np.ones(self.shape, np.float32), "grid": grid,
            "affine": np.eye(4), "world": np.eye(4), "native_mat": np.eye(4), "native_shape": self.shape,
            "params": dict(TestTheTissueMeans.params)}
        unified_seg.brain_mask.return_value = mask
        unified_seg.working_grid_maps.return_value = (None, None, None)
        apply_las = mock.Mock(return_value=nib.Nifti1Image(las, np.eye(4)))
        cat_surf = mock.Mock()
        cat_surf.vol_sanlm.side_effect = lambda data: data
        with mock.patch.object(frontend, "_t1prep", return_value=(cat_surf, apply_las, unified_seg)):
            out = frontend.segment(scan, template_dir="/t", space="template")
        np.testing.assert_array_equal(out["fractions"][dark], [1, 0, 0])
        self.assertEqual(float(out["label"][dark]), 1.0)
        self.assertGreater(float(out["normalised"][dark]), 0)
        # outside the brain there is nothing
        self.assertEqual(float(np.abs(out["fractions"][0]).max()), 0.0)
        self.assertEqual(float(out["normalised"][0].max()), 0.0)
        np.testing.assert_allclose(out["fractions"][mask].sum(-1), 1.0, atol=1e-6)

    def test_a_t1prep_without_the_deformation_is_named(self):
        scan = nib.Nifti1Image(np.ones(self.shape, np.float32), np.eye(4))

        def segment(denoised, raw, template_dir=None, verbose=False):
            return {}

        unified_seg = mock.Mock(segment=segment)
        with mock.patch.object(frontend, "_t1prep", return_value=(mock.Mock(vol_sanlm=lambda d: d), None, unified_seg)):
            with self.assertRaisesRegex(RuntimeError, "deformation"):
                frontend.segment(scan, template_dir="/t", deformation=True, space="template")

    def test_a_series_is_refused(self):
        scan = nib.Nifti1Image(np.ones(self.shape + (2,), np.float32), np.eye(4))
        with self.assertRaisesRegex(ValueError, "single 3-D volume"):
            self.run_segment(scan)


class TestTheGridOfTheScan(unittest.TestCase):
    """The scan cut into voxels of 0.5 mm along its own axes, where the working grid covers it."""

    @staticmethod
    def cover(size=(400, 400, 400)):
        """A working grid in template millimetres, larger than the scans below unless *size* says otherwise."""
        affine = np.diag([0.5, 0.5, 0.5, 1.0])
        affine[:3, 3] = -0.25 * (np.array(size) - 1)
        return affine, size

    def test_a_scan_of_1_mm_is_cut_into_two_by_two_by_two(self):
        mat = np.diag([1.0, 1.0, 1.0, 1.0])
        mat[:3, 3] = [-20, -25, -15]
        affine, shape = frontend.scan_grid(mat, (40, 50, 30), *self.cover(), np.eye(4))
        self.assertEqual(shape, (80, 100, 60))
        np.testing.assert_allclose(affine[:3, :3], 0.5 * np.eye(3))
        # the first fine voxel fills the first quarter of the first scan voxel: the faces coincide
        np.testing.assert_allclose(affine[:3, 3], mat[:3, 3] - 0.25)

    def test_a_scan_of_0_75_mm_and_another_storage_order(self):
        mat = np.diag([-0.75, 0.75, 0.75, 1.0])                         # first axis to the left
        mat[:3, 3] = [30, -40, -30]
        affine, shape = frontend.scan_grid(mat, (40, 48, 40), *self.cover(), np.eye(4))
        self.assertEqual(shape, (60, 72, 60))
        np.testing.assert_allclose(affine[:3, :3], np.diag([-0.5, 0.5, 0.5]))
        # the field of view is the same: from the outer face of the first voxel to that of the last
        for axis in range(3):
            self.assertAlmostEqual(affine[axis, 3] - affine[axis, axis] / 2, mat[axis, 3] - mat[axis, axis] / 2)

    def test_only_the_part_the_working_grid_covers_is_kept(self):
        mat = np.diag([1.0, 1.0, 1.0, 1.0])
        mat[:3, 3] = [-50, -50, -80]                                      # 100 x 100 x 160 mm, much neck
        world = np.diag([1.1, 1.1, 1.1, 1.0])                             # the template is larger than the head
        cover, size = self.cover(size=(200, 200, 200))                    # 100 mm in template mm: 91 mm of the scan
        affine, shape = frontend.scan_grid(mat, (100, 100, 160), cover, size, world)
        extent = 100 / 1.1
        for axis in range(3):
            self.assertGreaterEqual(shape[axis] * 0.5, extent)
            self.assertLessEqual(shape[axis] * 0.5, extent + 2 * (frontend.MARGIN + 1) * 0.5)
        np.testing.assert_allclose(affine[:3, :3], 0.5 * np.eye(3))
        # still on the lattice of the scan: every fourth face of the fine grid is a face of a scan voxel
        offset = (affine[:3, 3] - 0.25 - (mat[:3, 3] - 0.5)) / 0.5
        np.testing.assert_allclose(offset, np.round(offset), atol=1e-9)

    def test_a_sheared_header_gets_a_grid_at_right_angles(self):
        mat = np.array([[1.0, 0.1, 0, -20], [0, 1.0, 0, -25], [0, 0, 1.0, -15], [0, 0, 0, 1]])
        affine, _ = frontend.scan_grid(mat, (40, 50, 30), *self.cover(), np.eye(4))
        directions = affine[:3, :3] / 0.5
        np.testing.assert_allclose(directions.T @ directions, np.eye(3), atol=1e-9)

    def test_a_sheared_header_still_covers_what_the_working_grid_covers(self):
        # 1 mm voxels whose second axis leans by 6 degrees; the working grid lies inside the scan
        mat = np.array([[1.0, 0.1, 0, -60], [0, 1.0, 0, -70], [0, 0, 1.0, -50], [0, 0, 0, 1]])
        shape = (120, 140, 100)
        cover_affine, size = self.cover(size=(160, 180, 140))       # 80 x 90 x 70 mm about the origin
        world = np.diag([1.1, 1.1, 1.1, 1.0])
        affine, fine_shape = frontend.scan_grid(mat, shape, cover_affine, size, world)
        corners = np.array([[i, j, k, 1.0] for i in (-0.5, size[0] - 0.5) for j in (-0.5, size[1] - 0.5)
                            for k in (-0.5, size[2] - 0.5)])
        inside = np.linalg.solve(affine, np.linalg.solve(world, cover_affine @ corners.T)).T[:, :3]
        # every corner of the working grid lies inside the fine grid, with the margin
        self.assertGreaterEqual(float(inside.min()), frontend.MARGIN - 1.5)
        self.assertLessEqual(float((inside - (np.array(fine_shape) - 1)).max()), -(frontend.MARGIN - 1.5))

    def t1prep(self, world, cover_affine, cover_size):
        """A T1Prep that segments nothing, and what its unified segmentation was given."""
        calls = {}

        def segment(denoised, raw, template_dir=None, verbose=False, working_grid=None, world=None):
            calls["grid"], calls["world"] = working_grid, world
            n = working_grid.shape
            return {"prob": None, "label": np.zeros(n, np.float32), "corrected": np.ones(n, np.float32),
                    "grid": working_grid, "affine": np.eye(4), "world": world, "native_mat": raw.affine,
                    "native_shape": raw.shape, "params": dict(TestTheTissueMeans.params)}

        unified_seg = mock.Mock(segment=segment)
        unified_seg.unified_templates.return_value = ("templates", np.eye(4))
        unified_seg.affine_registration.return_value = (world, {})
        unified_seg.default_working_grid.return_value = nib.Nifti1Image(
            np.zeros(cover_size, np.uint8), cover_affine)
        unified_seg.brain_mask.side_effect = lambda prob: np.ones(calls["grid"].shape, bool)
        unified_seg.working_grid_maps.side_effect = lambda seg, mask: (
            None, nib.Nifti1Image(np.ones(seg["grid"].shape, np.float32), seg["grid"].affine),
            nib.Nifti1Image(np.ones(seg["grid"].shape, np.float32), seg["grid"].affine))
        apply_las = mock.Mock(side_effect=lambda brain, label, verbose=False, peaks=None: nib.Nifti1Image(
            np.full(brain.shape, 2 / 3, np.float32), brain.affine))
        cat_surf = mock.Mock()
        cat_surf.vol_sanlm.side_effect = lambda data: data
        return (cat_surf, apply_las, unified_seg), calls

    def test_the_front_end_makes_the_anatomy_on_it(self):
        shape = (12, 14, 10)
        stored = nib.Nifti1Image(np.ones(shape, np.float32), np.diag([1.0, 1.0, 1.0, 1.0]))
        world = np.diag([1.1, 1.1, 1.1, 1.0])
        cover_affine, cover_size = self.cover(size=(80, 80, 80))
        modules, calls = self.t1prep(world, cover_affine, cover_size)
        _, apply_las, unified_seg = modules
        with mock.patch.object(frontend, "_t1prep", return_value=modules):
            out = frontend.segment(stored, template_dir="/t")
        planned, _ = frontend.scan_grid(stored.affine, shape, cover_affine, cover_size, world)
        self.assertEqual(out["space"], "scan")
        self.assertEqual(out["grid"].shape, (24, 28, 20))
        np.testing.assert_allclose(out["grid"].affine, planned)
        # T1Prep got the same grid in template millimetres, and the front end says so
        np.testing.assert_allclose(calls["grid"].affine, world @ planned)
        np.testing.assert_allclose(out["template_affine"], world @ planned)
        self.assertEqual(out["template_grid"]["shape"], cover_size)
        # the registration is run once and handed on
        self.assertEqual(unified_seg.affine_registration.call_count, 1)
        np.testing.assert_array_equal(calls["world"], world)
        # LAS measured its millimetres on the scan: voxels of 0.5 mm
        brain = apply_las.call_args.args[0]
        np.testing.assert_allclose(brain.header.get_zooms()[:3], (0.5, 0.5, 0.5))
        np.testing.assert_allclose(out["fractions"][..., 1], 1.0, atol=1e-6)          # LAS said 2: grey matter
        with self.assertRaises(ValueError):
            frontend.segment(stored, template_dir="/t", space="mni")

        # a T1Prep that cannot take the registration registers the scan again, here a little elsewhere
        moved = world.copy()
        moved[0, 3] += 0.2
        segment = unified_seg.segment

        def older(denoised, raw, template_dir=None, verbose=False, working_grid=None):
            seg = segment(denoised, raw, template_dir, verbose, working_grid)
            seg["world"] = moved
            return seg

        unified_seg.segment = older
        with mock.patch.object(frontend, "_t1prep", return_value=modules):
            out = frontend.segment(stored, template_dir="/t")
        # the grid is then where T1Prep really put it in the scan, not where it was planned
        np.testing.assert_allclose(out["grid"].affine, np.linalg.inv(moved) @ world @ planned)
        self.assertGreater(float(np.abs(out["grid"].affine - planned).max()), 0.1)

    def test_a_scan_stored_in_another_order(self):
        # 0.75 mm, stored P, I, L, a scanner code: the grid follows RAS, the scan keeps its order and code
        shape = (16, 12, 14)
        affine = np.array([[0, 0, -0.75, 5.0], [-0.75, 0, 0, 6.0], [0, -0.75, 0, 4.0], [0, 0, 0, 1]])
        data = np.random.default_rng(3).random(shape).astype(np.float32)
        stored = nib.Nifti1Image(data, affine)
        stored.header.set_sform(affine, code=1)
        stored.header.set_qform(affine, code=1)
        world = np.diag([1.1, 1.1, 1.1, 1.0])
        cover_affine, cover_size = self.cover(size=(80, 80, 80))
        modules, _ = self.t1prep(world, cover_affine, cover_size)
        with mock.patch.object(frontend, "_t1prep", return_value=modules):
            out = frontend.segment(stored, template_dir="/t")
        canonical = nib.as_closest_canonical(stored)
        planned, planned_shape = frontend.scan_grid(canonical.affine, canonical.shape, cover_affine, cover_size,
                                                    world)
        self.assertEqual(out["grid"].shape, planned_shape)
        np.testing.assert_allclose(out["grid"].affine, planned)
        self.assertTrue(np.all(np.diag(out["grid"].affine)[:3] > 0))          # RAS
        # the first fine voxel fills the first corner of the field of view, two thirds of a scan voxel
        corner = canonical.affine @ np.array([-0.5, -0.5, -0.5, 1.0])
        np.testing.assert_allclose(planned @ np.array([-0.5, -0.5, -0.5, 1.0]), corner, atol=1e-9)
        np.testing.assert_array_equal(np.asarray(out["scan"].dataobj), data)
        np.testing.assert_allclose(out["scan"].affine, affine)
        self.assertEqual(int(out["scan"].header["sform_code"]), 1)


class TestTheScanInItsStorageOrder(unittest.TestCase):
    """The denoised scan goes back to the order, the affine and the NIfTI code it was stored with."""

    def test_order_affine_and_code_are_those_of_the_stored_scan(self):
        data = np.arange(4 * 5 * 6, dtype=np.float32).reshape(4, 5, 6)
        affine = np.array([[0, 0, -1.0, 10], [-1.0, 0, 0, 20], [0, -1.0, 0, 30], [0, 0, 0, 1]])   # P, I, L
        for sform, qform, expected in ((1, 1, 1), (2, 2, 2), (0, 1, 1), (0, 0, 1), (2, 1, 2)):
            with self.subTest(sform=sform, qform=qform):
                stored = nib.Nifti1Image(data, affine)
                stored.header.set_sform(affine, code=sform)
                stored.header.set_qform(affine, code=qform)
                canonical = nib.as_closest_canonical(stored)
                self.assertFalse(np.allclose(canonical.affine, affine))
                back = frontend._in_storage_order(np.asarray(canonical.dataobj), canonical.affine, stored)
                np.testing.assert_array_equal(np.asarray(back.dataobj), data)
                np.testing.assert_allclose(back.affine, stored.affine)
                self.assertEqual(int(back.header["sform_code"]), expected)
                self.assertEqual(frontend.world_code(back), expected)

    def test_an_image_without_a_nifti_header_is_in_the_scanner(self):
        self.assertEqual(frontend.world_code(object()), 1)


class TestWhatItNeeds(unittest.TestCase):
    def modules(self, apply_las):
        t1prep = mock.Mock(apply_LAS=apply_las, unified_seg=mock.Mock())
        return mock.patch.dict(sys.modules, {"cat_surf": mock.Mock(), "t1prep": t1prep})

    def test_a_t1prep_whose_las_takes_no_peaks_is_named(self):
        def apply_las(t1, label, verbose=False):
            return None

        with self.modules(apply_las):
            with self.assertRaisesRegex(RuntimeError, "peaks"):
                frontend._t1prep()

    def test_one_that_does_is_taken(self):
        def apply_las(t1, label, verbose=False, peaks=None):
            return None

        with self.modules(apply_las):
            self.assertIs(frontend._t1prep()[1], apply_las)

    def test_no_t1prep_at_all_is_named(self):
        with mock.patch.dict(sys.modules, {"cat_surf": mock.Mock(), "t1prep": None}):
            with self.assertRaisesRegex(RuntimeError, "t1prep.unified_seg"):
                frontend._t1prep()


if __name__ == "__main__":
    unittest.main()
