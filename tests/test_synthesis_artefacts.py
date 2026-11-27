"""Synthesis (:mod:`t1sim.synthesis`) and the artefacts of the image domain (:mod:`t1sim.artefacts`)."""

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

from t1sim import artefacts, synthesis  # noqa: E402

MEANS = {"CSF": 100.0, "GM": 300.0, "WM": 400.0}


class TestSynthesis(unittest.TestCase):
    def render(self, fractions, corrected, **kwargs):
        fractions = np.asarray(fractions, np.float32).reshape(-1, 1, 1, 3)
        corrected = np.asarray(corrected, np.float32).reshape(-1, 1, 1)
        return synthesis.render(fractions, MEANS, corrected, clip_percentile=None, **kwargs)[:, 0, 0]

    def test_a_pure_tissue_has_its_intensity(self):
        np.testing.assert_allclose(self.render(np.eye(3), [7, 7, 7]), [100, 300, 400])

    def test_a_mixture_is_the_weighted_mean(self):
        np.testing.assert_allclose(self.render([[0, .5, .5], [.25, .75, 0]], [9, 9]), [350, 250])

    def test_what_is_not_brain_keeps_the_scan(self):
        np.testing.assert_allclose(self.render([[0, 0, 0]], [222.0]), [222.0])

    def test_the_edge_of_the_brain_blends_the_two(self):
        np.testing.assert_allclose(self.render([[0.5, 0, 0]], [200.0]), [0.5 * 100 + 0.5 * 200])

    def test_a_voxel_without_data_is_nearly_zero(self):
        np.testing.assert_allclose(self.render([[0, 1, 0], [0, 0, 0]], [0.0, np.nan]), [1e-3, 1e-3])

    def test_a_negative_intensity_of_the_scan_becomes_zero(self):
        # the resampling of the scan leaves them in the background; a magnitude image has none
        np.testing.assert_allclose(self.render([[0, 0, 0], [0.5, 0, 0]], [-80.0, -80.0]), [0.0, 50.0])

    def test_the_brightest_voxels_are_limited(self):
        fractions = np.zeros((100, 1, 1, 3), np.float32)
        corrected = np.arange(1, 101, dtype=np.float32).reshape(100, 1, 1)
        out = synthesis.render(fractions, MEANS, corrected, clip_percentile=90)
        self.assertEqual(float(out.max()), 90.5)

    def test_three_means_in_the_order_of_the_fractions(self):
        fractions = np.eye(3, dtype=np.float32).reshape(3, 1, 1, 3)
        out = synthesis.render(fractions, [1, 2, 3], np.ones((3, 1, 1)), clip_percentile=None)
        np.testing.assert_allclose(out[:, 0, 0], [1, 2, 3])
        with self.assertRaises(ValueError):
            synthesis.render(fractions[..., :2], MEANS, np.ones((3, 1, 1)))


class TestTheSimulatedBiasField(unittest.TestCase):
    def test_it_is_the_field_of_the_matlab_tool(self):
        # add_simulated_bias_field of matlab/mri_simulate.m, rf.type = [2 0] and [3 7]
        field = artefacts.simulated_bias_field((60, 72, 50), 1.5, percent=30, strength=2, seed=0)
        got = [field[i] for i in ((0, 0, 0), (30, 36, 25), (59, 71, 49), (10, 60, 5), (45, 8, 40))]
        np.testing.assert_allclose(got, [1.04916, 0.992169, 1.025786, 1.125057, 0.912374], atol=2e-6)
        field = artefacts.simulated_bias_field((96, 80, 64), (0.5, 0.75, 1.0), percent=-20, strength=3, seed=7)
        np.testing.assert_allclose([field[i] for i in ((0, 0, 0), (48, 40, 32), (95, 79, 63))],
                                   [1.021693, 1.044086, 1.038065], atol=2e-6)

    def test_mean_one_and_the_amplitude_asked_for(self):
        field = artefacts.simulated_bias_field((40, 44, 36), 2.0, percent=40, strength=3, seed=5)
        self.assertAlmostEqual(float(field.mean()), 1.0, places=5)
        self.assertAlmostEqual(float(field.max() - field.min()), 0.40, places=5)

    def test_a_negative_amplitude_swaps_bright_and_dark(self):
        plus = artefacts.simulated_bias_field((30, 30, 30), 2.0, percent=30, seed=1)
        minus = artefacts.simulated_bias_field((30, 30, 30), 2.0, percent=-30, seed=1)
        np.testing.assert_allclose(plus - 1, -(minus - 1), atol=1e-6)

    def test_other_seed_other_field_and_a_strength_below_one_is_refused(self):
        a = artefacts.simulated_bias_field((24, 24, 24), 2.0, seed=1)
        b = artefacts.simulated_bias_field((24, 24, 24), 2.0, seed=2)
        self.assertGreater(float(np.abs(a - b).max()), 0.01)
        with self.assertRaises(ValueError):
            artefacts.simulated_bias_field((24, 24, 24), 2.0, strength=0)


class TestThePredefinedBiasField(unittest.TestCase):
    """A stored field in template space, with the defects of the fields that come with t1sim."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        shape = (12, 14, 10)
        x = np.arange(shape[0], dtype=np.float32)[:, None, None]
        self.stored = np.broadcast_to(5 + 0.05 * x, shape).copy()
        self.stored[:, :, 0] *= 0.5                 # an edge that was interpolated against zeros
        self.path = Path(self._tmp.name) / "rf100_X.nii"
        nib.save(nib.Nifti1Image(self.stored, np.eye(4)), self.path)
        self.iy = np.stack(np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij"), -1)

    def test_the_halved_edge_is_repaired(self):
        field = artefacts.template_bias_field("X", self.iy, np.eye(4), 30.0, folder=self._tmp.name)
        np.testing.assert_allclose(field[:, :, 0], field[:, :, 1], rtol=1e-6)
        self.assertAlmostEqual(float(field.mean()), 1.0, places=5)
        # 30% of a field that spans 0.55
        self.assertAlmostEqual(float(field.max() - field.min()), 0.3 * 0.55, places=4)
        np.testing.assert_array_equal(artefacts._repair_edges(self.stored)[:, :, 1:], self.stored[:, :, 1:])

    def test_where_the_deformation_is_undefined_the_field_goes_on(self):
        iy = self.iy.copy()
        iy[:2] = np.nan
        field = artefacts.template_bias_field(self.path, iy, np.eye(4), 30.0)
        self.assertTrue(np.isfinite(field).all())
        np.testing.assert_allclose(field[0], field[2], rtol=1e-6)       # the nearest value, not a step to 1
        np.testing.assert_allclose(field[1], field[2], rtol=1e-6)

    def test_a_negative_amplitude_inverts_it_and_a_missing_field_is_named(self):
        up = artefacts.template_bias_field(self.path, self.iy, np.eye(4), 30.0)
        down = artefacts.template_bias_field(self.path, self.iy, np.eye(4), -30.0)
        np.testing.assert_allclose(up + down, 2.0, atol=1e-5)
        with self.assertRaisesRegex(FileNotFoundError, "rf100_Q"):
            artefacts.template_bias_field("Q", self.iy, np.eye(4), 30.0, folder=self._tmp.name)


class TestContrastAndNoise(unittest.TestCase):
    def test_the_contrast_change_keeps_white_matter_and_raises_the_ratios(self):
        means = np.array([150.0, 380.0, 500.0, 800.0], np.float32)          # CSF, GM, WM and fat
        low = artefacts.change_contrast(means, 0.5, 500.0)
        high = artefacts.change_contrast(means, 1.5, 500.0)
        self.assertAlmostEqual(float(low[2]), 500.0, places=3)
        self.assertAlmostEqual(float(high[2]), 500.0, places=3)
        np.testing.assert_allclose(low / low[2], (means / 500) ** 0.5, rtol=1e-5)
        np.testing.assert_allclose(high / high[2], (means / 500) ** 1.5, rtol=1e-5)
        self.assertGreater(low[1] / low[2], means[1] / means[2])             # less contrast
        self.assertLess(high[1] / high[2], means[1] / means[2])              # more
        np.testing.assert_allclose(artefacts.change_contrast(means, 1, 500.0), means, rtol=1e-6)

    def test_the_contrast_change_of_one_intensity_and_of_what_cannot_be_mapped(self):
        self.assertAlmostEqual(artefacts.change_contrast(125.0, 0.5, 500.0), 250.0, places=4)
        self.assertIsInstance(artefacts.change_contrast(125.0, 0.5, 500.0), float)
        # below 0 counts as 0, where a power would not be a number
        np.testing.assert_array_equal(artefacts.change_contrast(np.array([-5.0, 0.0]), 0.5, 500.0), [0, 0])
        for exponent, white in ((0, 500.0), (-1, 500.0), (0.5, 0.0)):
            with self.assertRaises(ValueError):
                artefacts.change_contrast(np.ones(3), exponent, white)

    def test_rician_noise_has_the_standard_deviation_at_a_high_snr(self):
        image = np.full((40, 40, 40), 400.0, np.float32)
        noisy = artefacts.add_noise(image, 10.0, seed=1, ceiling=1e4)
        self.assertAlmostEqual(float(noisy.std()), 10.0, delta=0.15)
        self.assertAlmostEqual(float(noisy.mean()), 400.0 + 10.0 ** 2 / 800, delta=0.15)

    def test_the_background_gets_the_rayleigh_floor(self):
        noisy = artefacts.add_noise(np.zeros((40, 40, 40), np.float32), 10.0, seed=1, ceiling=1e4)
        self.assertAlmostEqual(float(noisy.mean()), 10.0 * np.sqrt(np.pi / 2), delta=0.15)
        self.assertGreaterEqual(float(noisy.min()), 0)

    def test_gaussian_noise_is_added_without_a_bias(self):
        image = np.full((30, 30, 30), 100.0, np.float32)
        noisy = artefacts.add_noise(image, 5.0, rician=False, seed=2)
        self.assertAlmostEqual(float(noisy.mean()), 100.0, delta=0.1)
        self.assertAlmostEqual(float(noisy.std()), 5.0, delta=0.1)

    def test_the_brightest_tissue_keeps_all_of_its_noise(self):
        # the image maximum is no ceiling: half of the noise lies above it
        image = np.full((40, 40, 40), 400.0, np.float32)
        noisy = artefacts.add_noise(image, 16.0, seed=3)
        self.assertAlmostEqual(float(noisy.std()), 16.0, delta=0.25)
        self.assertAlmostEqual(float(np.mean(noisy > 400.0)), 0.5, delta=0.02)
        limited = artefacts.add_noise(image, 16.0, seed=3, ceiling=400.0)
        self.assertLess(float(limited.std()), 0.7 * 16.0)       # what a ceiling at the maximum does

    def test_the_seed_decides_and_no_noise_changes_nothing(self):
        image = np.full((8, 8, 8), 50.0, np.float32)
        np.testing.assert_array_equal(artefacts.add_noise(image, 3.0, seed=4), artefacts.add_noise(image, 3.0, seed=4))
        self.assertFalse(np.array_equal(artefacts.add_noise(image, 3.0, seed=4), artefacts.add_noise(image, 3.0, seed=5)))
        np.testing.assert_array_equal(artefacts.add_noise(image, 0.0), image)


class TestTheAcquisition(unittest.TestCase):
    def setUp(self):
        self.affine = np.array([[0.5, 0, 0, -10.0], [0, 0.5, 0, -12.0], [0, 0, 0.5, -8.0], [0, 0, 0, 1]])
        self.shape = (40, 48, 32)

    def test_the_voxel_size_is_exact_and_the_grid_is_centred(self):
        acq = artefacts.Acquisition(self.shape, self.affine, (1.0, 1.5, 0.5))
        self.assertEqual(acq.shape, (20, 16, 32))
        np.testing.assert_allclose(np.sqrt((acq.affine[:3, :3] ** 2).sum(0)), [1.0, 1.5, 0.5])
        # both grids have the same centre in the world
        centre = lambda affine, shape: affine @ np.append((np.array(shape) - 1) / 2, 1)
        np.testing.assert_allclose(centre(acq.affine, acq.shape), centre(self.affine, self.shape))

    def test_a_voxel_size_that_does_not_divide_the_field_of_view(self):
        acq = artefacts.Acquisition(self.shape, self.affine, 0.7)
        self.assertEqual(acq.shape, (28, 34, 22))
        np.testing.assert_allclose(np.sqrt((acq.affine[:3, :3] ** 2).sum(0)), 0.7)

    def test_the_average_is_the_mean_over_the_box_of_a_voxel(self):
        rng = np.random.default_rng(0)
        volume = rng.random(self.shape).astype(np.float32)
        acq = artefacts.Acquisition(self.shape, self.affine, 1.0)
        blocks = volume.reshape(20, 2, 24, 2, 16, 2).mean((1, 3, 5))
        np.testing.assert_allclose(acq.average(volume), blocks, atol=1e-6)

    def test_the_average_keeps_the_volume_and_the_classes(self):
        rng = np.random.default_rng(1)
        fractions = rng.dirichlet((1, 1, 1), self.shape).astype(np.float32)
        acq = artefacts.Acquisition(self.shape, self.affine, (1.0, 1.0, 1.25))
        low = acq.average(fractions)
        self.assertEqual(low.shape, acq.shape + (3,))
        np.testing.assert_allclose(low.sum(-1), 1.0, atol=1e-5)
        # the part of the fine grid the new grid covers: here all of x and y, 30 of 32 planes of z
        covered = fractions[:, :, 1:31].sum((0, 1, 2), dtype=np.float64) * 0.5 ** 3
        np.testing.assert_allclose(low.sum((0, 1, 2), dtype=np.float64) * np.prod(acq.voxel), covered, rtol=1e-5)

    def test_acquiring_removes_what_the_new_voxels_cannot_hold(self):
        x = np.arange(self.shape[0])[:, None, None] * np.ones(self.shape)
        slow = np.cos(2 * np.pi * x / 20).astype(np.float32)             # 10 mm period
        fast = np.cos(2 * np.pi * x / 2.5).astype(np.float32)            # 1.25 mm period
        acq = artefacts.Acquisition(self.shape, self.affine, (1.0, 0.5, 0.5))
        # gone, but for what the end of the volume leaves of it next to a face
        self.assertLess(float(np.abs(acq.acquire(fast))[:-4].max()), 0.05)
        got = acq.acquire(slow)
        expected = np.cos(2 * np.pi * (acq.first[0] + 2 * np.arange(acq.shape[0])) / 20)
        # a slow wave stays, to the little its mirror image at the end of the volume adds
        np.testing.assert_allclose(got[:, 5, 5], expected, atol=0.015)

    def test_what_touches_one_end_of_the_volume_does_not_ring_into_the_other(self):
        # a Fourier transform takes the volume as periodic: the neck would show up above the head
        image = np.zeros((96, 32, 32), np.float32)
        image[:10] = 400.0
        acq = artefacts.Acquisition(image.shape, np.diag([0.5, 0.5, 0.5, 1.0]), (1.0, 0.5, 0.5))
        far_end = acq.acquire(image)[-4:]
        self.assertLess(float(np.abs(far_end).max()), 0.01 * 400)

    def test_a_smooth_image_survives_the_way_there_and_back(self):
        x, y, z = np.meshgrid(*[np.linspace(-1, 1, n) for n in self.shape], indexing="ij")
        image = (100 * np.exp(-(x ** 2 + y ** 2 + z ** 2) * 2)).astype(np.float32)
        acq = artefacts.Acquisition(self.shape, self.affine, 1.0)
        back = acq.to_fine(acq.acquire(image))
        self.assertEqual(back.shape, self.shape)
        self.assertLess(float(np.abs(back - image)[4:-4, 4:-4, 4:-4].max()), 1.0)

    def test_the_same_voxel_size_changes_nothing(self):
        image = np.random.default_rng(2).random(self.shape).astype(np.float32)
        acq = artefacts.Acquisition(self.shape, self.affine, 0.5)
        np.testing.assert_allclose(acq.acquire(image), image, atol=1e-6)
        np.testing.assert_allclose(acq.affine, self.affine)

    def test_a_volume_of_another_size_is_refused(self):
        acq = artefacts.Acquisition(self.shape, self.affine, 1.0)
        with self.assertRaises(ValueError):
            acq.acquire(np.zeros((4, 4, 4), np.float32))


def rotation_about_x(angle):
    """A rotation about the first axis, by *angle* radians."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])


class TestTheObliqueAcquisition(unittest.TestCase):
    """A grid of other voxels, another orientation and a larger field of view than that of the image."""

    @classmethod
    def setUpClass(cls):
        n = 64
        axis = np.arange(n) - (n - 1) / 2
        x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
        r = np.sqrt(x ** 2 + y ** 2 + (1.2 * z) ** 2) * 0.5              # mm, 0.5 mm voxels
        ramp = lambda radius: np.clip((radius - r) / 0.5 + 0.5, 0, 1)
        gm, wm = ramp(11.0) - ramp(8.0), ramp(8.0)
        cls.fractions = np.stack([ramp(14.0) - ramp(11.0), gm, wm], -1).astype(np.float32)
        cls.image = (cls.fractions @ np.float32([100, 300, 400])).astype(np.float32)
        cls.affine = np.diag([0.5, 0.5, 0.5, 1.0])
        cls.affine[:3, 3] = -0.5 * (n - 1) / 2
        # 0.8 x 0.85 x 0.9 mm, turned about all three axes, shifted, and larger than the image
        a, b, c = 0.2, -0.1, 0.15
        turn = (np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
                @ np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
                @ np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]]))
        cls.shape = (56, 60, 52)
        cls.target = np.eye(4)
        cls.target[:3, :3] = turn @ np.diag([0.8, 0.85, 0.9])
        cls.target[:3, 3] = -cls.target[:3, :3] @ ((np.array(cls.shape) - 1) / 2) + [1.3, -0.7, 0.4]
        cls.oblique = artefacts.ObliqueAcquisition(cls.image.shape, cls.affine, cls.shape, cls.target)

    def test_on_a_grid_along_the_axes_it_is_the_acquisition_of_that_grid(self):
        straight = artefacts.Acquisition(self.image.shape, self.affine, 1.0)
        oblique = artefacts.ObliqueAcquisition(self.image.shape, self.affine, straight.shape, straight.affine)
        inner = (slice(2, -2),) * 3                                  # the two treat the border differently
        np.testing.assert_allclose(oblique.acquire(self.image)[inner], straight.acquire(self.image)[inner], atol=0.01)
        np.testing.assert_allclose(oblique.average(self.fractions), straight.average(self.fractions), atol=1e-5)
        np.testing.assert_allclose(oblique.voxel, 1.0)

    def test_the_average_keeps_the_volumes_and_the_fractions(self):
        averaged = self.oblique.average(self.fractions)
        self.assertEqual(averaged.shape, self.shape + (3,))
        fine = self.fractions.sum((0, 1, 2), dtype=np.float64) * 0.5 ** 3
        coarse = averaged.sum((0, 1, 2), dtype=np.float64) * abs(np.linalg.det(self.target[:3, :3]))
        np.testing.assert_allclose(coarse, fine, rtol=2e-4)
        self.assertGreaterEqual(float(averaged.min()), 0.0)
        self.assertLessEqual(float(averaged.sum(-1).max()), 1.0 + 1e-5)
        # one volume without classes is averaged the same way
        np.testing.assert_allclose(self.oblique.average(self.fractions[..., 1]), averaged[..., 1], atol=1e-6)

    @staticmethod
    def integral(lattice, volume, slices):
        """The mean over every new voxel, from the boxes of *slices* slices along every axis of it.

        Slices much smaller than a fine voxel lie in one fine voxel nearly
        all, and their boxes then give the value there: the result converges
        to the integral, here to 0.0003 with 6 slices for voxels near cubes.
        """
        def many(linear):
            lengths = np.sqrt((linear ** 2).sum(0))
            return slices * np.maximum(1, np.round(lengths / lengths.min())).astype(int)

        with mock.patch.object(artefacts, "_box_counts", many):
            return lattice.average(volume)

    def test_the_average_is_the_mean_over_the_voxel_and_not_its_centre(self):
        from scipy import ndimage

        gm = self.fractions[..., 1]
        exact = self.integral(self.oblique, gm, 6)
        averaged = self.oblique.average(gm)
        boundary = (exact > 0.02) & (exact < 0.98)
        self.assertLess(float(np.abs(averaged - exact).max()), 0.008)                # 0.0055
        self.assertLess(float(np.sqrt(np.mean((averaged - exact)[boundary] ** 2))), 0.002)   # 0.0014
        # the fine voxel at the centre of every new voxel is far from it
        centre = ndimage.affine_transform(gm, self.oblique.linear, offset=self.oblique.offset, order=0,
                                          output_shape=self.shape)
        self.assertGreater(float(np.abs(centre - exact).max()), 0.3)

    def test_a_turned_voxel_four_times_as_long_as_wide_is_cut_into_more_slices_along_its_length(self):
        a, b = 0.4, 0.25
        turn = (np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
                @ np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]]))
        shape = (60, 60, 15)
        target = np.eye(4)
        target[:3, :3] = turn @ np.diag([0.5, 0.5, 2.0])
        target[:3, 3] = -target[:3, :3] @ ((np.array(shape) - 1) / 2) + [0.3, -0.2, 0.1]
        lattice = artefacts.ObliqueAcquisition(self.image.shape, self.affine, shape, target)
        np.testing.assert_array_equal(artefacts._box_counts(lattice.linear), [3, 3, 12])
        gm = self.fractions[..., 1]
        exact = self.integral(lattice, gm, 4)
        averaged = lattice.average(gm)
        boundary = (exact > 0.02) & (exact < 0.98)
        # 3 slices along every axis were 0.034 and 0.007 off
        self.assertLess(float(np.abs(averaged - exact).max()), 0.003)                # 0.0010
        self.assertLess(float(np.sqrt(np.mean((averaged - exact)[boundary] ** 2))), 0.001)   # 0.0003

    def test_the_average_is_exact_on_any_grid_along_the_axes(self):
        # 0.8 x 1 x 1.3 mm: no whole number of fine voxels, where counting at points is furthest off
        straight = artefacts.Acquisition(self.image.shape, self.affine, (0.8, 1.0, 1.3))
        oblique = artefacts.ObliqueAcquisition(self.image.shape, self.affine, straight.shape, straight.affine)
        np.testing.assert_allclose(oblique.average(self.fractions), straight.average(self.fractions), atol=1e-5)

    def test_the_image_is_the_mixture_of_the_averaged_fractions_away_from_the_boundaries(self):
        image = self.oblique.acquire(self.image)
        averaged = self.oblique.average(self.fractions)
        white = averaged[..., 2] > 0.9999
        self.assertGreater(int(white.sum()), 500)
        np.testing.assert_allclose(image[white], 400.0, rtol=0.06)       # it rings, as an acquisition does
        self.assertAlmostEqual(float(image[white].mean()), 400.0, delta=2.0)
        # nothing of the image where the new grid lies beyond it
        self.assertEqual(float(np.abs(image[self.oblique.coverage() == 0]).max()), 0.0)

    def test_the_coverage_and_a_map_that_is_continued(self):
        coverage = self.oblique.coverage(ramp=4.0)
        self.assertEqual((float(coverage.min()), float(coverage.max())), (0.0, 1.0))
        centre = tuple(n // 2 for n in self.shape)
        self.assertEqual(float(coverage[centre]), 1.0)
        self.assertEqual(float(coverage[0, 0, 0]), 0.0)
        field = np.broadcast_to(np.linspace(0.8, 1.2, 64, dtype=np.float32)[:, None, None], self.image.shape)
        cut, continued = self.oblique.sample(field), self.oblique.sample(field, extend=True)
        self.assertEqual(float(cut[0, 0, 0]), 0.0)
        self.assertGreaterEqual(float(continued.min()), 0.8 - 1e-6)
        np.testing.assert_allclose(cut[coverage == 1], continued[coverage == 1], atol=1e-6)

    def test_the_coverage_falls_only_where_something_lies_beyond(self):
        # a new grid that ends inside the fine one at the lower faces and beyond it at the upper ones
        target = np.diag([1.0, 1.0, 1.0, 1.0])
        target[:3, 3] = self.affine[:3, 3] + 4.0                    # starts 8 fine voxels in
        lattice = artefacts.ObliqueAcquisition(self.image.shape, self.affine, (40, 40, 40), target)
        coverage = lattice.coverage(ramp=6.0)
        # 8 fine voxels from the lower faces: no ramp there, nothing lies below them
        self.assertEqual(float(coverage[0, 20, 20]), 1.0)
        # the upper faces: new voxels lie beyond them
        self.assertEqual(float(coverage[-1, 20, 20]), 0.0)
        position = 8 + 2 * np.arange(40)                             # fine coordinate of the new voxels
        expected = np.clip((63.5 - position) / 6.0, 0, 1)            # from the outer face of the last fine voxel
        np.testing.assert_allclose(coverage[:, 20, 20], expected, atol=1e-6)

    def test_a_grid_that_ends_at_the_faces_of_the_fine_one_is_covered(self):
        face = self.affine[:3, 3] - 0.25                             # the outer face of the first fine voxel, in mm
        # voxels of 0.4 mm over the same field of view: the outermost centres lie between a face and a fine centre
        fine = np.diag([0.4, 0.4, 0.4, 1.0])
        fine[:3, 3] = face + 0.2
        # the fine grid itself, its affine off by rounding
        same = self.affine.copy()
        same[:3, :3] = same[:3, :3] @ rotation_about_x(1e-9)
        same[:3, 3] -= 1e-9
        for target, shape in ((fine, (80, 80, 80)), (same, self.image.shape)):
            lattice = artefacts.ObliqueAcquisition(self.image.shape, self.affine, shape, target)
            self.assertEqual(float(lattice.coverage().min()), 1.0)
            image = lattice.acquire(self.image + 50)                 # nothing of it is 0 by itself
            self.assertGreater(float(image.min()), 0.0)

    def test_a_volume_of_another_size_is_refused(self):
        with self.assertRaises(ValueError):
            self.oblique.acquire(np.zeros((10, 10, 10), np.float32))
        with self.assertRaises(ValueError):
            self.oblique.average(np.zeros((10, 10, 10, 3), np.float32))


if __name__ == "__main__":
    unittest.main()
