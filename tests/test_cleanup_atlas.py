"""Cleanup of the ground truth (:mod:`t1sim.cleanup`) and maps of the template space (:mod:`t1sim.atlas`)."""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import atlas, cleanup  # noqa: E402

try:
    import cat_surf  # noqa: F401
    HAVE_CAT_SURF = True
except ImportError:
    HAVE_CAT_SURF = False

AFFINE = np.eye(4)
N = 64


def _radius():
    axis = np.arange(N) - (N - 1) / 2
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.sqrt(x ** 2 + y ** 2 + z ** 2), x


class TestApply(unittest.TestCase):
    def setUp(self):
        self.label = np.array([2.0, 2.4, 2.0, 2.6], np.float32).reshape(4, 1, 1)
        self.fractions = np.array([[0, 1, 0], [0, .6, .4], [0, 1, 0], [0, .4, .6]], np.float32).reshape(4, 1, 1, 3)

    def test_a_vessel_becomes_csf(self):
        to_csf = np.array([True, False, False, False]).reshape(4, 1, 1)
        fractions, label = cleanup.apply(self.fractions, self.label, to_csf)
        np.testing.assert_array_equal(fractions[0, 0, 0], [1, 0, 0])
        self.assertEqual(label[0, 0, 0], 1.0)
        np.testing.assert_array_equal(fractions[1:], self.fractions[1:])

    def test_periventricular_voxels_become_the_mixture_of_csf_and_wm_with_their_intensity(self):
        pve = np.array([False, True, True, False]).reshape(4, 1, 1)
        means = np.array([150.0, 380.0, 500.0], np.float32)      # GM is not half way
        fractions, label = cleanup.apply(self.fractions, self.label, None, pve, means=means)
        np.testing.assert_array_equal(fractions[..., 1][pve], 0)
        np.testing.assert_allclose(fractions.sum(-1), 1, atol=1e-6)
        # the image is rendered from the fractions as they were: the new ones give the same intensity
        np.testing.assert_allclose(fractions @ means, self.fractions @ means, rtol=1e-5)
        np.testing.assert_allclose(fractions[2, 0, 0], [120 / 350, 0, 230 / 350], atol=1e-6)
        # the label is that of the new fractions, no longer 2
        np.testing.assert_allclose(label, fractions[..., 0] + 2 * fractions[..., 1] + 3 * fractions[..., 2], atol=1e-6)
        self.assertAlmostEqual(float(label[2, 0, 0]), 3 - 2 * 120 / 350, places=5)
        np.testing.assert_array_equal(fractions[~pve], self.fractions[~pve])
        with self.assertRaisesRegex(ValueError, "means"):
            cleanup.apply(self.fractions, self.label, None, pve)

    def test_a_voxel_without_fractions_becomes_csf(self):
        # the brightest vessels have no label and no fractions before the cleanup
        fractions, label = self.fractions.copy(), self.label.copy()
        fractions[0], label[0] = 0, 0
        to_csf = np.array([True, False, False, False]).reshape(4, 1, 1)
        fractions, label = cleanup.apply(fractions, label, to_csf)
        np.testing.assert_array_equal(fractions[0, 0, 0], [1, 0, 0])
        self.assertEqual(label[0, 0, 0], 1.0)

    def test_the_input_is_not_changed(self):
        before = self.fractions.copy()
        cleanup.apply(self.fractions, self.label, np.ones((4, 1, 1), bool))
        np.testing.assert_array_equal(self.fractions, before)


@unittest.skipUnless(HAVE_CAT_SURF, "cat_surf is needed for the path through the tissue")
class TestDura(unittest.TestCase):
    """A ball of WM in 3 mm of cortex; beyond 4 mm of CSF a sheet as bright as cortex."""

    def setUp(self):
        r, self.x = _radius()
        self.r = r
        label = np.zeros((N, N, N), np.float32)
        label[r < 26] = 1
        label[(r >= 19) & (r < 21)] = 2          # the sheet
        label[r < 15] = 2
        label[r < 12] = 3
        self.label = label
        self.sheet = (r >= 19) & (r < 21)
        self.cortex = (r >= 12) & (r < 15)
        self.regions = np.zeros((N, N, N), np.uint8)

    def detect(self, **kwargs):
        options = dict(vessels=0, dura=4.0, pve=False)
        options.update(kwargs)
        return cleanup.detect(self.label, self.label, AFFINE, self.regions, **options)

    def test_the_sheet_beyond_the_csf_is_dura_and_the_cortex_is_not(self):
        found = self.detect()["dura"]
        self.assertGreater(float(found[self.sheet].mean()), 0.99)
        self.assertEqual(int(found[self.cortex].sum()), 0)
        self.assertEqual(int(found[~self.sheet].sum()), 0)

    def test_a_bridge_of_tissue_makes_the_near_part_cortex(self):
        self.label[(self.x > 0) & (self.r >= 15) & (self.r < 19) & (np.abs(_radius()[0] - self.x) < 2)] = 2
        centre = N // 2
        # along the bridge the sheet begins 7 mm from the white matter
        self.assertTrue(self.detect(dura=4.0)["dura"][centre + 19, centre, centre])
        found = self.detect(dura=8.0)["dura"]
        self.assertFalse(found[centre + 19, centre, centre])
        self.assertTrue(found[centre - 20, centre, centre])       # the far side has no bridge

    def test_a_cerebellum_is_spared(self):
        self.regions[self.x > 0] = atlas.CEREBELLUM[0]
        found = self.detect()["dura"]
        self.assertEqual(int(found[self.x > 0].sum()), 0)
        self.assertGreater(int(found[self.x < 0].sum()), 0)

    def test_it_can_be_switched_off_and_a_negative_length_is_refused(self):
        self.assertEqual(int(self.detect(dura=0)["dura"].sum()), 0)
        with self.assertRaises(ValueError):
            self.detect(dura=-1)


@unittest.skipUnless(HAVE_CAT_SURF, "cat_surf is needed for the growth from the brightest voxels")
class TestBrightVessels(unittest.TestCase):
    """A vessel in the cortex whose core is too bright to have a label, and a bright spot deep in the white matter."""

    def setUp(self):
        r, x = _radius()
        label = np.zeros((N, N, N), np.float32)
        label[r < 26] = 1
        label[r < 15] = 2
        label[r < 12] = 3
        normalised = label.copy()
        c = N // 2
        self.shell = np.zeros(label.shape, bool)
        self.shell[c + 12:c + 15, c - 1:c + 2, c - 1:c + 2] = True          # in the cortex
        self.core = np.zeros(label.shape, bool)
        self.core[c + 13, c, c] = True
        self.spot = np.zeros(label.shape, bool)
        self.spot[c - 1:c + 2, c - 1:c + 2, c - 1:c + 2] = True              # in the middle of the white matter
        normalised[self.shell | self.spot] = 3.7
        label[self.shell] = 3                       # brighter than white matter is white matter in the label
        normalised[self.core] = 4.5
        label[self.core] = 0                        # and from 4 on there is no label
        self.label, self.normalised = label, normalised
        # T1Prep's own detection finds nothing here: only the rule of the brightness is tested
        fake = mock.MagicMock()
        fake.vessel_weight.return_value = np.zeros(label.shape, np.float32)
        package = mock.MagicMock(vessels=fake)
        patcher = mock.patch.dict(sys.modules, {"t1prep": package, "t1prep.vessels": fake})
        patcher.start()
        self.addCleanup(patcher.stop)

    def detect(self, **kwargs):
        options = dict(vessels=0.5, dura=0, pve=False)
        options.update(kwargs)
        return cleanup.detect(self.label, self.normalised, AFFINE, np.zeros(self.label.shape, np.uint8),
                              **options)["vessels"]

    def test_the_vessel_is_found_with_its_core_and_the_white_matter_is_left(self):
        found = self.detect()
        self.assertTrue(found[self.shell].all())
        self.assertTrue(found[self.core].all())      # the hole in the ground truth
        self.assertEqual(int(found[self.spot].sum()), 0)
        self.assertEqual(int(found[~self.shell].sum()), 0)

    def test_the_cleaned_ground_truth_has_no_hole(self):
        fractions = np.zeros(self.label.shape + (3,), np.float32)
        for k in range(3):
            fractions[..., k] = self.label == k + 1
        cleaned, label = cleanup.apply(fractions, self.label, self.detect())
        np.testing.assert_array_equal(cleaned[self.core], [[1, 0, 0]])
        np.testing.assert_array_equal(label[self.shell], 1.0)

    def test_it_can_be_switched_off(self):
        self.assertEqual(int(self.detect(vessels=0).sum()), 0)


class TestPeriventricular(unittest.TestCase):
    """A ventricle inside white matter whose rim has the intensity of grey matter."""

    def setUp(self):
        r, _ = _radius()
        label = np.zeros((N, N, N), np.float32)
        label[r < 28] = 1
        label[r < 24] = 2                         # cortex, 4 mm
        label[r < 20] = 3
        label[r < 7] = 2                          # the rim: one voxel between WM and CSF
        label[r < 6] = 1                          # the ventricle
        self.label, self.r = label, r
        self.regions = np.where(r < 6, atlas.VENTRICLES[2], 0).astype(np.uint8)

    def test_the_rim_is_found_and_the_cortex_is_not(self):
        found = cleanup.detect(self.label, self.label, AFFINE, self.regions, vessels=0, dura=0, pve=True)["pve"]
        rim = (self.r >= 6) & (self.r < 7)
        self.assertGreater(float(found[rim].mean()), 0.9)
        self.assertEqual(int(found[(self.r >= 20) & (self.r < 24)].sum()), 0)

    def test_not_inside_the_deep_grey_matter(self):
        regions = self.regions.copy()
        regions[(self.r >= 6) & (self.r < 12)] = atlas.THALAMUS[0]
        found = cleanup.detect(self.label, self.label, AFFINE, regions, vessels=0, dura=0, pve=True)["pve"]
        self.assertEqual(int(found.sum()), 0)


class TestClean(unittest.TestCase):
    def test_the_deformation_is_required(self):
        with self.assertRaisesRegex(ValueError, "deformation"):
            cleanup.clean({"fractions": None, "label": None})


class TestToGrid(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.data = rng.integers(1, 200, (10, 12, 8)).astype(np.uint8)
        self.image = nib.Nifti1Image(self.data, np.diag([2.0, 2.0, 2.0, 1.0]))
        grid = np.stack(np.meshgrid(np.arange(10.0), np.arange(12.0), np.arange(8.0), indexing="ij"), -1)
        self.iy = grid.astype(np.float32)

    def test_the_identity_returns_the_map(self):
        out = atlas.to_grid(self.image, self.iy, self.image.affine)
        np.testing.assert_array_equal(out, self.data)
        self.assertEqual(out.dtype, np.uint8)

    def test_the_grids_of_map_and_templates_may_differ(self):
        # templates with 4 mm voxels: template voxel t lies at map voxel 2 t
        out = atlas.to_grid(self.image, self.iy[:5, :6, :4], np.diag([4.0, 4.0, 4.0, 1.0]))
        np.testing.assert_array_equal(out, self.data[::2, ::2, ::2])

    def test_where_the_deformation_is_undefined_or_the_map_ends_is_the_fill_value(self):
        iy = self.iy.copy()
        iy[0, 0, 0] = np.nan
        iy[1, 0, 0] = (50, 0, 0)
        out = atlas.to_grid(self.image, iy, self.image.affine)
        self.assertEqual((out[0, 0, 0], out[1, 0, 0]), (0, 0))
        self.assertEqual(out[2, 0, 0], self.data[2, 0, 0])

    def test_a_continuous_map_is_interpolated(self):
        ramp = nib.Nifti1Image(np.arange(10, dtype=np.float32)[:, None, None] * np.ones((10, 12, 8), np.float32),
                               np.eye(4))
        iy = self.iy.copy()
        iy[..., 0] += 0.25
        out = atlas.to_grid(ramp, iy, np.eye(4), order=1, fill=np.nan)
        np.testing.assert_allclose(out[:9, 3, 3], np.arange(9) + 0.25, atol=1e-5)
        self.assertTrue(np.isnan(out[9, 3, 3]))                 # beyond the last voxel

    def test_regions_fill_the_voxels_without_one(self):
        regions = np.zeros((9, 9, 9), np.uint8)
        regions[1, 4, 4], regions[7, 4, 4] = 5, 9
        out = atlas.fill_from_nearest(regions, 1.0)
        self.assertEqual((out[0, 4, 4], out[3, 4, 4], out[6, 4, 4]), (5, 5, 9))
        hole = np.zeros((9, 9, 9), bool)
        hole[1, 4, 4] = True
        self.assertTrue(np.all(atlas.fill_from_nearest(regions, 1.0, hole=hole) == 9))


class TestTheTemplateFolder(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.folder = Path(self._tmp.name)
        nib.save(nib.Nifti1Image(np.zeros((2, 2, 2), np.uint8), np.eye(4)), self.folder / "Hammers.nii.gz")
        (self.folder / "Hammers.csv").write_text("ROIid;ROIname;Vgm\n28;FL middle frontal gyrus L;1\n29;other;1\n")
        patch = mock.patch.object(atlas, "template_folder", return_value=self.folder)
        patch.start()
        self.addCleanup(patch.stop)

    def test_an_atlas_is_found_by_its_name_whatever_the_case(self):
        self.assertEqual(atlas.find("hammers"), self.folder / "Hammers.nii.gz")
        self.assertEqual(atlas.find(self.folder / "Hammers.nii.gz"), self.folder / "Hammers.nii.gz")
        with self.assertRaises(FileNotFoundError):
            atlas.find("nothing")

    def test_its_regions_come_from_the_table_next_to_it(self):
        self.assertEqual(atlas.regions("Hammers"), {28: "FL middle frontal gyrus L", 29: "other"})

    def test_the_parcels_of_the_lobes_are_cortical(self):
        for numbers in (atlas.OCCIPITAL, atlas.FRONTAL):
            self.assertTrue(all(n >= atlas.CORTEX_FROM for n in numbers))
        self.assertFalse(set(atlas.OCCIPITAL) & set(atlas.FRONTAL))


if __name__ == "__main__":
    unittest.main()
