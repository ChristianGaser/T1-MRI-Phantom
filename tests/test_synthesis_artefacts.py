"""Synthesis (:mod:`t1sim.synthesis`) and the artefacts of the image domain (:mod:`t1sim.artefacts`)."""

import sys
import unittest
from pathlib import Path

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


class TestContrastAndNoise(unittest.TestCase):
    def test_the_contrast_change_keeps_the_range_and_the_order(self):
        image = np.linspace(10, 500, 50, dtype=np.float32).reshape(50, 1, 1)
        low = artefacts.change_contrast(image, 0.5)
        self.assertAlmostEqual(float(low.min()), 10, places=3)
        self.assertAlmostEqual(float(low.max()), 500, places=3)
        self.assertTrue(np.all(np.diff(low[:, 0, 0]) > 0))
        self.assertGreater(float(low[25, 0, 0]), float(image[25, 0, 0]))     # 0.5 lifts the middle
        np.testing.assert_array_equal(artefacts.change_contrast(image, 1), image)

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

    def test_a_contrast_change_maps_a_reference_intensity_with_the_image(self):
        image = np.linspace(0, 500, 51, dtype=np.float32).reshape(51, 1, 1)
        out, white = artefacts.change_contrast(image, 0.5, reference=400.0)
        self.assertAlmostEqual(white, float(out[40, 0, 0]), places=3)
        self.assertGreater(white, 400.0)
        same, white = artefacts.change_contrast(image, 1.0, reference=400.0)
        self.assertEqual(white, 400.0)

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
        self.assertLess(float(np.abs(acq.acquire(fast)).max()), 1e-3)
        got = acq.acquire(slow)
        expected = np.cos(2 * np.pi * (acq.first[0] + 2 * np.arange(acq.shape[0])) / 20)
        # the first and the last voxel see the border of the volume
        np.testing.assert_allclose(got[1:-1, 5, 5], expected[1:-1], atol=1e-3)

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


if __name__ == "__main__":
    unittest.main()
