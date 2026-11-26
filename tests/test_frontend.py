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
            out = frontend.segment(scan, template_dir="/templates")
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
            plain = frontend.segment(scan, template_dir="/t")
            self.assertFalse(calls["deformation"])
            self.assertNotIn("iy", plain)
            out = frontend.segment(scan, template_dir="/t", deformation=True)
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
            out = frontend.segment(scan, template_dir="/t")
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
                frontend.segment(scan, template_dir="/t", deformation=True)

    def test_a_series_is_refused(self):
        scan = nib.Nifti1Image(np.ones(self.shape + (2,), np.float32), np.eye(4))
        with self.assertRaisesRegex(ValueError, "single 3-D volume"):
            self.run_segment(scan)


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
