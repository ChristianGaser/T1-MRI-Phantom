"""Changes of the anatomy (:mod:`t1sim.anatomy`) on a ball of white matter in a shell of cortex."""

import sys
import unittest
from pathlib import Path

import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import anatomy, frontend  # noqa: E402

AFFINE = np.eye(4)                      # 1 mm voxels
N, WM_R, GM_R, CSF_R = 64, 12.0, 15.0, 24.0


def _phantom():
    """WM ball, 3 mm of cortex, CSF, nothing: label, fractions, an image and the radius."""
    axis = np.arange(N) - (N - 1) / 2
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # linear partial volume across one voxel at every boundary
    label = (np.clip(WM_R + 0.5 - r, 0, 1) + np.clip(GM_R + 0.5 - r, 0, 1)
             + np.clip(CSF_R + 0.5 - r, 0, 1)).astype(np.float32)
    fractions = frontend.fractions_from_label(label)
    image = (fractions @ np.array([100, 300, 400], np.float32)).astype(np.float32)
    return label, fractions, image, r, x


def _gm_along(fractions, direction):
    """mm of grey matter along the x axis through the centre, on one side."""
    c = N // 2
    line = fractions[c:, c, c, 1] if direction > 0 else fractions[:c, c, c, 1]
    return float(line.sum())


class TestAtrophyByScaling(unittest.TestCase):
    def test_only_voxels_with_csf_change_and_the_sum_stays(self):
        _, fractions, _, r, x = _phantom()
        roi = x > 0
        out = anatomy.atrophy_by_scaling(fractions, roi, 2.0)
        np.testing.assert_allclose(out.sum(-1), fractions.sum(-1), atol=1e-6)
        no_csf = roi & (fractions[..., 0] == 0)
        np.testing.assert_array_equal(out[no_csf], fractions[no_csf])
        np.testing.assert_array_equal(out[~roi], fractions[~roi])
        mixed = roi & (fractions[..., 0] > 0) & (fractions[..., 1] > 0)
        self.assertTrue(np.all(out[mixed][:, 1] < fractions[mixed][:, 1]))


class TestAtrophyByDeformation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.label, cls.fractions, cls.image, cls.r, cls.x = _phantom()
        cls.roi = (cls.x > 0) & (cls.r < CSF_R)
        cls.results = {mode: anatomy.atrophy(cls.fractions, cls.image, AFFINE, cls.roi, -0.10, mode=mode)
                       for mode in ("thinning", "volume")}

    def test_the_tissue_of_the_region_loses_what_was_asked_for(self):
        for mode, res in self.results.items():
            with self.subTest(mode=mode):
                self.assertAlmostEqual(res["achieved"], -0.10, delta=0.004)
                self.assertGreater(float(res["jacobian"].min()), 0)

    def test_grey_matter_is_lost_and_csf_takes_its_place(self):
        before = self.fractions.sum((0, 1, 2), dtype=np.float64)
        for mode, res in self.results.items():
            with self.subTest(mode=mode):
                after = res["fractions"].sum((0, 1, 2), dtype=np.float64)
                self.assertLess(after[1], before[1] - 100)               # mm^3 of GM
                self.assertGreater(after[0], before[0] + 50)
                np.testing.assert_allclose(res["fractions"].sum(-1).max(), 1.0, atol=1e-4)

    def test_the_cortex_of_the_region_gets_thinner(self):
        # A closed shell around a ball can only lose volume by getting
        # thinner, so both kinds thin it here; in a folded cortex a change of
        # volume also shortens it.
        before = _gm_along(self.fractions, +1)
        for mode, res in self.results.items():
            with self.subTest(mode=mode):
                change = _gm_along(res["fractions"], +1) / before - 1
                self.assertLess(change, -0.05)
                self.assertGreater(change, -0.16)

    def test_the_other_side_keeps_its_cortex(self):
        before = _gm_along(self.fractions, -1)
        for mode, res in self.results.items():
            with self.subTest(mode=mode):
                self.assertAlmostEqual(_gm_along(res["fractions"], -1), before, delta=0.03 * before)

    def test_thinning_leaves_the_white_matter_where_it_is(self):
        wm_before = float(self.fractions[..., 2].sum(dtype=np.float64))
        wm_after = float(self.results["thinning"]["fractions"][..., 2].sum(dtype=np.float64))
        self.assertAlmostEqual(wm_after, wm_before, delta=0.01 * wm_before)

    def test_the_image_and_other_volumes_follow_the_same_deformation(self):
        res = self.results["thinning"]
        again = anatomy.warp(self.fractions, res["displacement"])
        np.testing.assert_allclose(again, res["fractions"], atol=1e-5)
        np.testing.assert_allclose(anatomy.warp(self.image, res["displacement"]), res["corrected"], atol=1e-3)
        np.testing.assert_allclose(res["label"], frontend.label_from_fractions(res["fractions"]), atol=1e-5)

    def test_a_positive_change_thickens(self):
        res = anatomy.atrophy(self.fractions, self.image, AFFINE, self.roi, 0.10, mode="thinning")
        self.assertGreater(_gm_along(res["fractions"], +1), _gm_along(self.fractions, +1) * 1.05)

    def test_what_cannot_be_done_is_refused(self):
        for kwargs, message in ((dict(change=0), "not 0"), (dict(change=-0.8), "between"),
                                (dict(change=-0.1, mode="other"), "mode"),
                                (dict(change=-0.1, tissues=("CSF",)), "tissues")):
            with self.subTest(kwargs=kwargs):
                with self.assertRaisesRegex(ValueError, message):
                    anatomy.atrophy(self.fractions, self.image, AFFINE, self.roi, **kwargs)
        with self.assertRaisesRegex(ValueError, "empty"):
            anatomy.atrophy(self.fractions, self.image, AFFINE, np.zeros_like(self.roi), -0.1)


class TestThickness(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.label, cls.fractions, _, cls.r, _ = _phantom()
        cls.regions = np.where(cls.r < CSF_R + 1, 100, 0).astype(np.uint8)      # all of it one cortical parcel

    def test_the_cortex_gets_the_thickness_asked_for(self):
        for mm in (2.0, 4.0):
            with self.subTest(mm=mm):
                fractions, label = anatomy.thickness(self.fractions, self.label, self.regions, AFFINE, mm)
                for direction in (+1, -1):
                    self.assertAlmostEqual(_gm_along(fractions, direction), mm, delta=0.35)
                np.testing.assert_allclose(fractions.sum(-1)[self.r < CSF_R - 2], 1.0, atol=1e-5)
                np.testing.assert_allclose(label, frontend.label_from_fractions(fractions), atol=1e-5)

    def test_the_white_matter_stays(self):
        fractions, _ = anatomy.thickness(self.fractions, self.label, self.regions, AFFINE, 2.0)
        before = float(self.fractions[..., 2].sum(dtype=np.float64))
        self.assertAlmostEqual(float(fractions[..., 2].sum(dtype=np.float64)), before, delta=0.04 * before)

    def test_a_structure_that_is_no_cortex_keeps_its_labels(self):
        regions = self.regions.copy()
        half = np.arange(N)[:, None, None] * np.ones((N, N, N)) > N / 2 + 6
        regions[half & (self.r < CSF_R + 1)] = 36                              # "caudate"
        fractions, _ = anatomy.thickness(self.fractions, self.label, regions, AFFINE, 5.0)
        self.assertAlmostEqual(_gm_along(fractions, +1), _gm_along(self.fractions, +1), delta=0.3)
        self.assertGreater(_gm_along(fractions, -1), 4.0)

    def test_one_value_or_three(self):
        with self.assertRaises(ValueError):
            anatomy.thickness(self.fractions, self.label, self.regions, AFFINE, (2.0, 3.0))


if __name__ == "__main__":
    unittest.main()
