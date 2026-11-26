"""The shared volume operations (:mod:`t1sim._volume`): sizes in mm, as SPM and CAT12 define them."""

import sys
import unittest
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import _volume as vol  # noqa: E402


class TestSmoothing(unittest.TestCase):
    def test_the_kernel_sums_to_one_and_has_the_width_asked_for(self):
        kernel = vol._smoothing_kernel(6.0)
        self.assertAlmostEqual(kernel.sum(), 1.0, places=12)
        x = np.arange(kernel.size) - kernel.size // 2
        # a Gaussian of that FWHM, widened by the triangle of the interpolation (variance 1/6)
        self.assertAlmostEqual(np.sum(kernel * x ** 2), (6.0 / 2.3548) ** 2 + 1 / 6, delta=0.02)

    def test_the_width_is_given_in_mm(self):
        point = np.zeros((41, 41, 41), np.float32)
        point[20, 20, 20] = 1
        out = vol.smooth(point, 4.0, (0.5, 1.0, 2.0))
        widths = [np.sum(out.sum(tuple(a for a in range(3) if a != axis)) > 1e-4) for axis in range(3)]
        self.assertGreater(widths[0], widths[1])
        self.assertGreater(widths[1], widths[2])
        self.assertAlmostEqual(float(out.sum()), 1.0, places=4)


class TestDistances(unittest.TestCase):
    def setUp(self):
        self.mask = np.zeros((30, 30, 30), bool)
        self.mask[10:20, 10:20, 10:20] = True

    def test_the_distance_is_euclidean_in_mm(self):
        d = vol.distance(self.mask, (0.5, 1.0, 2.0))
        self.assertEqual(d[15, 15, 15], 0)
        self.assertAlmostEqual(float(d[6, 15, 15]), 4 * 0.5)
        self.assertAlmostEqual(float(d[15, 15, 22]), 3 * 2.0)
        self.assertAlmostEqual(float(d[7, 7, 15]), np.hypot(3 * 0.5, 3 * 1.0), places=5)

    def test_an_empty_mask_is_infinitely_far(self):
        self.assertTrue(np.isinf(vol.distance(np.zeros((4, 4, 4), bool))).all())

    def test_dilation_and_erosion_are_balls_in_mm(self):
        grown = vol.dilate(self.mask, 2.0, 1.0)
        self.assertTrue(grown[8, 15, 15] and not grown[7, 15, 15])
        self.assertFalse(grown[8, 8, 15])                       # the corner is further than 2 mm
        shrunk = vol.erode(self.mask, 2.0, 1.0)
        self.assertTrue(shrunk[12, 15, 15] and not shrunk[11, 15, 15])
        np.testing.assert_array_equal(vol.erode(self.mask, 0, 1.0), self.mask)

    def test_closing_fills_a_gap_and_keeps_the_outline(self):
        gap = self.mask.copy()
        gap[14:16] = False
        closed = vol.closing(gap, 1.5, 1.0)
        self.assertTrue(closed[15, 15, 15])
        self.assertFalse(closed[9, 15, 15])
        # also when the mask touches the border of the volume
        edge = np.zeros((10, 10, 10), bool)
        edge[:4] = True
        np.testing.assert_array_equal(vol.closing(edge, 2.0, 1.0), edge)

    def test_opening_removes_what_is_thin(self):
        thin = self.mask.copy()
        thin[20:26, 15, 15] = True
        opened = vol.opening(thin, 1.5, 1.0)
        self.assertFalse(opened[23, 15, 15])
        self.assertTrue(opened[15, 15, 15])

    def test_the_nearest_voxel_of_the_mask(self):
        index = vol.nearest_index(self.mask)
        self.assertEqual([int(i[2, 15, 15]) for i in index], [10, 15, 15])


class TestTheRest(unittest.TestCase):
    def test_the_largest_component_only(self):
        mask = np.zeros((12, 12, 12), bool)
        mask[1:6, 1:6, 1:6] = True
        mask[8:10, 8:10, 8:10] = True
        mask[6, 6, 6] = True                                    # touches the cube by a corner only
        out = vol.largest_component(mask)
        self.assertEqual(int(out.sum()), 125)

    def test_grey_closing_fills_a_dark_gap(self):
        volume = np.ones((15, 15, 15), np.float32)
        volume[7] = 0.2
        np.testing.assert_allclose(vol.grey_close(volume, 2), 1.0)
        dark = np.zeros((15, 15, 15), np.float32)
        dark[7, 7, 7] = 1                                       # a bright voxel alone stays alone
        self.assertEqual(float(vol.grey_close(dark, 2).sum()), 1.0)

    def test_the_bounding_box_with_a_margin(self):
        mask = np.zeros((20, 20, 20), bool)
        mask[5:8, 0:3, 17:20] = True
        box = vol.bounding_box(mask, 2)
        self.assertEqual([(s.start, s.stop) for s in box], [(3, 10), (0, 5), (15, 20)])
        with self.assertRaises(ValueError):
            vol.bounding_box(np.zeros((3, 3, 3), bool))

    def test_the_percentile_ignores_what_is_not_finite(self):
        values = np.array([1.0, 2.0, 3.0, 4.0, np.nan, np.inf])
        self.assertEqual(vol.percentile(values, 50), 2.5)

    def test_the_voxel_size_of_a_grid(self):
        affine = np.diag([0.5, 0.75, 2.0, 1.0])
        np.testing.assert_allclose(vol.grid_voxel_size(affine), [0.5, 0.75, 2.0])
        with self.assertRaises(ValueError):
            vol.voxel_size([1, 2])


if __name__ == "__main__":
    unittest.main()
