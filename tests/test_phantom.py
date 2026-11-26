"""The thickness phantom (:mod:`t1sim.phantom`): is the thickness the requested one?

Small phantoms at 1 mm, so that the whole file runs in a few seconds.  Tested
is what the phantom promises: the two cortical boundaries are the requested
thickness apart, for both folding patterns; a flat boundary gets the linear
partial volume ramp of one voxel; the volumes are those of the geometry; the
sidecars record it; and the evaluation measures all of this from the label
image alone.

That label and T1w image equal those of ``matlab/thickness_phantom.m`` voxel
by voxel needs MATLAB and is not tested here.
"""

import functools
import importlib.util
import json
import os
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import nibabel as nib
import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import t1sim  # noqa: E402
from t1sim import phantom  # noqa: E402

#: A folded sphere that 40 voxels of 1 mm hold, with folds that 1 mm resolves:
#: the phantom for everything that is about accuracy.
SMALL = {"vx": 1.0, "dim": 40, "radius": 11.0, "amplitude": 1.5, "csf": 1.5}
FOLDS = {"cartesian": {"wavelength": 10.0}, "spherical": {"fold": "spherical", "folds": (3, 4)}}

#: The smallest phantom that still has all tissues, for everything else.
TINY = {"vx": 1.0, "dim": 24, "radius": 6.0, "amplitude": 1.0, "wavelength": 8.0, "csf": 1.0}

#: A sphere without folds.  The odd dimension puts a voxel on its centre, so
#: the three axes through the centre run through voxel centres.
SPHERE = {"vx": 1.0, "dim": 31, "radius": 7.2, "amplitude": 0.0, "csf": 1.6}


@functools.cache
def tiny(thickness=2.0):
    """The tiny phantom of one thickness; its arrays are shared, do not change them."""
    return phantom.thickness_phantom(thickness, **TINY)


def radius_of(volume_ml: float) -> float:
    """Radius in mm of the sphere of that volume."""
    return (3 * volume_ml * 1000 / (4 * np.pi)) ** (1 / 3)


class TestTheThickness(unittest.TestCase):
    """The thickness recovered from the label image alone is the requested one."""

    thickness = (1.5, 2.25, 3.0)

    @classmethod
    def setUpClass(cls):
        cls.results, cls.infos = {}, {}
        for fold, pattern in FOLDS.items():
            labels, _, infos = phantom.thickness_phantom(cls.thickness, **SMALL, **pattern)
            cls.infos[fold] = infos
            cls.results[fold] = phantom.evaluate_phantom(labels, infos, pbt=False)

    def check(self, fold):
        for known, result in zip(self.thickness, self.results[fold]):
            outer = result["geom"]["outer"]
            self.assertEqual(result["thickness"], known)
            self.assertGreater(outer["n"], 10000)
            # a twentieth of a voxel
            self.assertLess(abs(outer["bias"]), 0.05, (fold, known))
            self.assertLess(abs(outer["median"] - known), 0.05, (fold, known))
            self.assertLess(outer["sd"], 0.12, (fold, known))
            self.assertLess(outer["p5"], outer["median"])
            self.assertLess(outer["median"], outer["p95"])
            self.assertAlmostEqual(outer["mean"] - known, outer["bias"])

    def test_with_cartesian_folds(self):
        self.check("cartesian")

    def test_with_spherical_folds(self):
        self.check("spherical")

    def test_shallow_folds_bury_next_to_no_sulcus(self):
        for results in self.results.values():
            for result in results:
                self.assertLess(result["geom"]["buried"], 0.03)
                self.assertLess(abs(result["geom"]["inner"]["bias"]), 0.05)

    def test_the_label_volumes_are_those_the_phantom_reports(self):
        for fold, results in self.results.items():
            for result, info in zip(results, self.infos[fold]):
                self.assertIs(result["volume"]["ref"], info["volume_ref"])
                for tissue in ("csf", "gm", "wm"):
                    self.assertAlmostEqual(
                        result["volume"]["label"][tissue] / info["volume_pve"][tissue], 1, places=6)

    def test_a_thick_cortex_buries_a_narrow_sulcus(self):
        labels, _, infos = phantom.thickness_phantom(
            [1.0, 3.0], vx=1.0, dim=36, radius=8.0, amplitude=2.0, wavelength=7.0, csf=1.5)
        thin, thick = phantom.evaluate_phantom(labels, infos, pbt=False)
        self.assertLess(thin["geom"]["buried"], 0.1)
        self.assertGreater(thick["geom"]["buried"], 0.3)
        # outwards the distance stays the thickness, back to the WM it no longer is
        self.assertLess(abs(thick["geom"]["outer"]["bias"]), 0.1)
        self.assertGreater(thick["geom"]["inner"]["bias"], 0.3)


class TestASphereWithoutFolds(unittest.TestCase):
    """Along an axis through the centre the boundaries are flat and their position is known."""

    thickness = 2.3

    def profile_error(self, pve_steps, pve_range):
        label, _, info = phantom.thickness_phantom(
            self.thickness, pve_steps=pve_steps, pve_range=pve_range, **SPHERE)
        # The object ends at the last fine sample inside the analytic radius
        # and the WM/GM boundary lies wm_offset beyond it; the other two
        # boundaries follow at the thickness and the width of the CSF layer.
        vxf = info["vx"] / info["supersample"]
        analytic = info["radius"] - info["wm_offset"] + 0.5 * vxf
        boundary = np.floor(analytic / vxf) * vxf + info["wm_offset"]
        self.assertLessEqual(abs(boundary - info["radius"]), 0.5 * vxf)

        centre = SPHERE["dim"] // 2
        distance = np.abs(np.arange(SPHERE["dim"]) - centre) * info["vx"]

        def ramp(position):
            return np.clip(0.5 + (position - distance) / (pve_range * info["vx"]), 0, 1)

        ideal = (ramp(boundary) + ramp(boundary + self.thickness)
                 + ramp(boundary + self.thickness + info["csf"]))
        profiles = (label[:, centre, centre], label[centre, :, centre], label[centre, centre, :])
        self.assertTrue(np.any(profiles[0] % 1 > 0), "no partial volume on the axis")
        return max(np.abs(profile - ideal).max() for profile in profiles)

    def test_a_flat_boundary_gets_the_linear_ramp_of_one_voxel(self):
        # the midpoint rule is exact up to half a step of the average
        self.assertLessEqual(self.profile_error(15, 1.0), 0.5 / 15 + 1e-6)
        self.assertLessEqual(self.profile_error(101, 1.0), 0.5 / 101 + 1e-6)

    def test_the_ramp_is_as_wide_as_pve_range(self):
        self.assertLessEqual(self.profile_error(25, 0.5), 0.5 / 25 + 1e-6)
        self.assertLessEqual(self.profile_error(40, 0.8), 0.5 / 40 + 1e-6)

    def test_offsets_beyond_wm_offset_are_a_warning(self):
        # the distance map is 0 inside the object and cannot be shifted there
        with self.assertWarnsRegex(UserWarning, "no voxel will be pure WM"):
            label, _, _ = phantom.thickness_phantom(self.thickness, pve_range=1.5, **SPHERE)
        self.assertLess(float(label.max()), 3.0)

    def test_the_volumes_are_those_of_shells_of_the_requested_thickness(self):
        _, _, info = phantom.thickness_phantom(self.thickness, **SPHERE)
        for key, tolerance in (("volume_ref", 0.05), ("volume_pve", 0.1)):
            volume = info[key]
            wm = radius_of(volume["wm"])
            gm = radius_of(volume["wm"] + volume["gm"])
            csf = radius_of(volume["wm"] + volume["gm"] + volume["csf"])
            self.assertLess(abs(wm - SPHERE["radius"]), 1 / 3, key)
            self.assertLess(abs(gm - wm - self.thickness), tolerance, key)
            self.assertLess(abs(csf - gm - SPHERE["csf"]), tolerance, key)
        for tissue in ("csf", "gm", "wm"):
            self.assertLess(abs(info["volume_pve"][tissue] / info["volume_ref"][tissue] - 1), 0.03)

    def test_the_t1w_image_mixes_the_tissue_intensities(self):
        intensities = (0.3, 0.7, 1.1)
        label, t1w, _ = phantom.thickness_phantom(self.thickness, intensities=intensities, **SPHERE)
        self.assertEqual(label.dtype, np.float32)
        self.assertEqual(t1w.dtype, np.float32)
        for value, intensity in zip((0, 1, 2, 3), (0.0,) + intensities):
            self.assertTrue(np.any(label == value))
            np.testing.assert_allclose(t1w[label == value], intensity, rtol=1e-6)
        # neighbouring tissues mix linearly, so the fractions are the triangular weights
        mixture = sum(intensity * np.clip(1 - np.abs(label - value), 0, None)
                      for value, intensity in zip((1, 2, 3), intensities))
        np.testing.assert_allclose(t1w, mixture, atol=1e-6)
        self.assertEqual(float(label.max()), 3.0)
        self.assertEqual(float(label.min()), 0.0)


class TestTheConstruction(unittest.TestCase):
    def test_it_is_deterministic(self):
        first = phantom.thickness_phantom(2.0, **SMALL, **FOLDS["spherical"])
        second = phantom.thickness_phantom(2.0, **SMALL, **FOLDS["spherical"])
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[1], second[1])
        self.assertEqual(first[2], second[2])

    def test_a_series_shares_the_distance_map_and_changes_nothing(self):
        labels, images, infos = phantom.thickness_phantom([1.5, 2.0], **TINY)
        self.assertEqual((len(labels), len(images), len(infos)), (2, 2, 2))
        for i, thickness in enumerate((1.5, 2.0)):
            label, image, info = tiny(thickness)
            np.testing.assert_array_equal(labels[i], label)
            np.testing.assert_array_equal(images[i], image)
            self.assertEqual(infos[i], info)
        # the WM does not depend on the thickness, the GM grows with it
        self.assertEqual(infos[0]["volume_ref"]["wm"], infos[1]["volume_ref"]["wm"])
        self.assertLess(infos[0]["volume_ref"]["gm"], infos[1]["volume_ref"]["gm"])

    def test_one_value_gives_one_phantom_and_a_sequence_gives_lists(self):
        label, _, info = tiny()
        self.assertEqual(label.shape, (24, 24, 24))
        self.assertIsInstance(info, dict)
        labels, _, infos = phantom.thickness_phantom([2.0], **TINY)
        self.assertIsInstance(labels, list)
        np.testing.assert_array_equal(labels[0], label)
        self.assertEqual(infos[0], info)

    def test_it_is_centred_in_a_volume_that_is_no_cube(self):
        label, t1w, info = phantom.thickness_phantom(2.0, **dict(SMALL, dim=(36, 40, 44)))
        self.assertEqual(label.shape, (36, 40, 44))
        self.assertEqual(t1w.shape, (36, 40, 44))
        self.assertEqual(info["dim"], [36, 40, 44])
        for axis in range(3):
            np.testing.assert_array_equal(label, np.flip(label, axis))
            # the CSF layer is not cut open at the border
            self.assertFalse(np.take(label, 0, axis).any())
        # x and y are interchangeable in the cartesian pattern, the volume is not
        np.testing.assert_array_equal(label[:, 2:-2, :], label[:, 2:-2, :].transpose(1, 0, 2))

    def test_the_info_holds_the_geometry(self):
        _, _, info = phantom.thickness_phantom(
            2.0, **dict(TINY, wavelength=9.0), pve_steps=9, pve_range=0.8, intensities=(1, 2, 3))
        expected = {"thickness": 2.0, "vx": 1.0, "dim": [24, 24, 24], "radius": 6.0,
                    "amplitude": 1.0, "fold": "cartesian", "wavelength": 9.0, "folds": [6, 6],
                    "csf": 1.0, "wm_offset": 0.5, "supersample": 3, "pve_steps": 9,
                    "pve_range": 0.8, "intensities": [1.0, 2.0, 3.0],
                    "files": {"label": None, "t1w": None}}
        self.assertEqual({key: info[key] for key in expected}, expected)
        self.assertEqual(set(info) - set(expected), {"volume_ref", "volume_pve"})
        self.assertEqual(json.loads(json.dumps(info)), info)

    def test_the_distance_is_the_euclidean_distance_to_the_nearest_sample(self):
        index = np.indices((9, 10, 11)).reshape(3, -1).T
        obj = ((index * [7, 11, 13]).sum(axis=1) % 17 == 0).reshape(9, 10, 11)
        inside = index[obj.ravel()]
        brute = np.sqrt(((index[:, None, :] - inside[None, :, :]) ** 2).sum(axis=-1)).min(axis=1)
        dist = phantom._distance_to(obj, 0.4)
        self.assertEqual(dist.dtype, np.float32)
        np.testing.assert_allclose(dist.ravel(), 0.4 * brute, rtol=1e-6)
        self.assertTrue(np.isinf(phantom._distance_to(np.zeros((3, 3, 3), bool), 1.0)).all())


class TestTheFiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.folder.cleanup)
        cls.outdir = os.path.join(cls.folder.name, "new", "folder")
        cls.label, cls.t1w, cls.info = phantom.thickness_phantom(
            1.25, vx=0.5, dim=(28, 30, 32), radius=3.5, amplitude=0.5, wavelength=4.0, csf=0.6,
            outdir=cls.outdir, name="test")
        cls.base = os.path.join(cls.outdir, "test_res-0p5mm_desc-Thickness13mm")
        cls.result = phantom.evaluate_phantom(cls.label, cls.info, pbt=False)

    def sidecar(self, suffix):
        with open(self.base + suffix + ".json", encoding="utf-8") as handle:
            return json.load(handle)

    def test_two_values_with_one_name_are_not_written_over_each_other(self):
        with self.assertRaisesRegex(ValueError, "same file name"):
            phantom.thickness_phantom([1.25, 1.3], dim=40, radius=5, amplitude=1, wavelength=6, csf=1,
                                      outdir=self.outdir)
        # without files there is nothing to collide
        labels, _, _ = phantom.thickness_phantom([1.25, 1.3], dim=40, radius=5, amplitude=1, wavelength=6, csf=1)
        self.assertEqual(len(labels), 2)

    def test_the_names_follow_mri_simulate(self):
        self.assertEqual(self.info["files"], {"label": self.base + "_dseg.nii.gz",
                                              "t1w": self.base + "_T1w.nii.gz"})
        self.assertEqual(sorted(os.listdir(self.outdir)), sorted(
            os.path.basename(self.base) + end
            for end in ("_dseg.nii.gz", "_dseg.json", "_T1w.nii.gz", "_T1w.json")))

    def test_the_images_come_back_unchanged(self):
        for key, data in (("label", self.label), ("t1w", self.t1w)):
            image = nib.load(self.info["files"][key])
            self.assertEqual(image.get_data_dtype(), np.float32)
            np.testing.assert_array_equal(np.asarray(image.dataobj), data)
            np.testing.assert_allclose(image.header.get_zooms(), 0.5)
            # the world origin is the centre of the volume
            centre = (np.array(data.shape) - 1) / 2
            np.testing.assert_allclose(image.affine @ np.append(centre, 1), [0, 0, 0, 1], atol=1e-6)

    def test_the_label_sidecar_records_the_geometry(self):
        meta = self.sidecar("_dseg")
        self.assertIn("partial volume label", meta.pop("Description"))
        self.assertEqual(meta, {
            "GeneratedBy": [{"Name": "t1sim", "Version": t1sim.__version__}],
            "Manual": False,
            "Labels": {"CSF": 1, "GM": 2, "WM": 3},
            "PhantomParameters": {
                "Thickness": 1.25, "CSFLayer": 0.6, "Radius": 3.5, "Amplitude": 0.5,
                "Fold": "cartesian", "FoldWavelength": 4.0, "WMOffset": 0.25,
                "VoxelSize": [0.5, 0.5, 0.5], "Supersample": 3, "PVESteps": 15, "PVERange": 1.0},
        })

    def test_the_t1w_sidecar_names_its_source_and_the_intensities(self):
        meta = self.sidecar("_T1w")
        self.assertEqual(meta["Sources"], ["bids::test_res-0p5mm_desc-Thickness13mm_dseg.nii.gz"])
        self.assertEqual(meta["GeneratedBy"], self.sidecar("_dseg")["GeneratedBy"])
        parameters = self.sidecar("_dseg")["PhantomParameters"]
        self.assertEqual(meta["PhantomParameters"],
                         dict(parameters, Intensities=[30.0, 100.0, 150.0]))
        self.assertNotIn("Labels", meta)

    def test_spherical_folds_are_recorded_by_their_number(self):
        label, t1w, info = phantom.thickness_phantom(
            2.0, **dict(TINY, fold="spherical", folds=(3, 4)))
        files = phantom.write_phantom(label, t1w, info, os.path.join(self.folder.name, "spherical"))
        self.assertEqual(os.path.basename(files["label"]),
                         "thicknessPhantom_res-1mm_desc-Thickness20mm_dseg.nii.gz")
        with open(files["label"].replace(".nii.gz", ".json"), encoding="utf-8") as handle:
            parameters = json.load(handle)["PhantomParameters"]
        self.assertEqual(parameters["Fold"], "spherical")
        self.assertEqual(parameters["FoldNumber"], [3, 4])
        self.assertNotIn("FoldWavelength", parameters)

    def test_the_evaluation_reads_thickness_and_voxel_size_from_the_files(self):
        from_file = phantom.evaluate_phantom(self.info["files"]["label"], pbt=False)
        self.assertEqual(from_file["name"], os.path.basename(self.base) + "_dseg.nii.gz")
        self.assertEqual(self.result["name"], "label image")
        self.assertEqual((from_file["thickness"], from_file["vx"]), (1.25, 0.5))
        self.assertEqual(from_file["geom"], self.result["geom"])
        self.assertEqual(from_file["volume"]["label"], self.result["volume"]["label"])
        self.assertIsNone(from_file["volume"]["ref"])
        self.assertLess(abs(from_file["geom"]["outer"]["bias"]), 0.05 * 0.5)

    def test_an_explicit_thickness_wins_over_the_sidecar(self):
        result = phantom.evaluate_phantom(self.info["files"]["label"], thickness=1.0, pbt=False)
        self.assertEqual(result["thickness"], 1.0)
        self.assertAlmostEqual(result["geom"]["outer"]["bias"],
                               self.result["geom"]["outer"]["bias"] + 0.25)


class TestTheEvaluation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.label, _, cls.info = tiny()
        cls.result = phantom.evaluate_phantom(cls.label, cls.info, pbt=False)

    def test_the_statistics_follow_the_matlab_conventions(self):
        stats = phantom._stats([4.0, 1.0, np.nan, 3.0, 2.0, np.inf], 2.0)
        self.assertEqual(stats["n"], 4)
        self.assertEqual(stats["mean"], 2.5)
        self.assertEqual(stats["bias"], 0.5)
        # the sample of rank round(q n), not an interpolated quantile
        self.assertEqual((stats["median"], stats["p5"], stats["p95"]), (2.0, 1.0, 4.0))
        self.assertAlmostEqual(stats["sd"], np.sqrt(5 / 3))
        self.assertAlmostEqual(stats["rmse"], np.sqrt(1.5))
        self.assertEqual(phantom._stats([3.0], 2.0)["sd"], 0.0)
        empty = phantom._stats([np.nan], 2.0)
        self.assertEqual(empty["n"], 0)
        self.assertTrue(all(np.isnan(empty[key]) for key in empty if key != "n"))

    def test_supersampling_is_trilinear_and_keeps_the_voxel_centres(self):
        x, y, z = np.meshgrid(np.arange(5.0), np.arange(6.0), np.arange(7.0), indexing="ij")
        volume = (10 + 2 * x - y + 0.5 * z).astype(np.float32)
        for ss in (1, 3, 5):
            fine = phantom._upsample(volume, ss)
            half = (ss - 1) // 2
            self.assertEqual(fine.shape, (5 * ss, 6 * ss, 7 * ss))
            np.testing.assert_array_equal(fine[half::ss, half::ss, half::ss], volume)
            # a linear image stays linear between the outermost voxel centres
            inner = fine[tuple(slice(half, n - half) for n in fine.shape)]
            i, j, k = np.meshgrid(*(np.arange(n) / ss for n in inner.shape), indexing="ij")
            np.testing.assert_allclose(inner, 10 + 2 * i - j + 0.5 * k, atol=1e-5)
            # and beyond them is background
            self.assertEqual(np.count_nonzero(fine), inner.size)

    def test_a_finer_measurement_grid_is_more_accurate(self):
        coarse = phantom.evaluate_phantom(self.label, self.info, supersample=1, pbt=False)
        even = phantom.evaluate_phantom(self.label, self.info, supersample=2, pbt=False)
        self.assertLess(abs(self.result["geom"]["outer"]["bias"]),
                        abs(coarse["geom"]["outer"]["bias"]))
        self.assertGreater(self.result["geom"]["outer"]["n"], 4 * coarse["geom"]["outer"]["n"])
        # an even factor is increased by one
        self.assertEqual(even["geom"], self.result["geom"])

    def test_a_thickness_map_is_compared_with_the_known_thickness(self):
        gm = (self.label > 1.5) & (self.label < 2.5)
        thickness_map = np.full(self.label.shape, 2.1, np.float32)
        undefined = gm & (np.indices(self.label.shape)[0] < 12)
        thickness_map[undefined] = 0
        position = np.where(gm & ~undefined & (self.label > 1.9), 0.5, 0.0).astype(np.float32)

        result = phantom.evaluate_phantom(self.label, self.info, thickness_map=thickness_map)
        self.assertEqual(result["geom"], self.result["geom"])
        self.assertEqual(result["measure"]["name"], "thickness map")
        self.assertIsNone(result["measure"]["central"])
        self.assertGreater(np.count_nonzero(undefined), 100)
        self.assertEqual(result["measure"]["gm"]["n"], np.count_nonzero(gm & ~undefined))
        self.assertAlmostEqual(result["measure"]["gm"]["bias"], 0.1, places=5)
        self.assertAlmostEqual(result["measure"]["gm"]["sd"], 0.0, places=5)

        result = phantom.evaluate_phantom(self.label, self.info, thickness_map=thickness_map,
                                          position_map=position)
        self.assertGreater(np.count_nonzero(position), 100)
        self.assertEqual(result["measure"]["central"]["n"], np.count_nonzero(position))
        self.assertAlmostEqual(result["measure"]["central"]["rmse"], 0.1, places=5)

    def test_without_a_measure_only_the_phantom_is_checked(self):
        self.assertIsNone(self.result["measure"])
        self.assertEqual(self.result["volume"]["ref"], self.info["volume_ref"])

    def test_a_missing_cat_surf_is_a_warning_and_no_error(self):
        with mock.patch.dict(sys.modules, {"cat_surf": None}), \
                self.assertWarnsRegex(UserWarning, "cat_surf"):
            result = phantom.evaluate_phantom(self.label, self.info)
        self.assertIsNone(result["measure"])
        self.assertEqual(result["geom"], self.result["geom"])

    @unittest.skipUnless(importlib.util.find_spec("cat_surf"), "cat_surf is not installed")
    def test_the_projection_based_thickness_of_cat_surf(self):
        plain = phantom.evaluate_phantom(self.label, self.info,
                                         pbt_options={"correct_thickness": 0.0})["measure"]
        shifted = phantom.evaluate_phantom(self.label, self.info,
                                           pbt_options={"correct_thickness": 0.4})["measure"]
        self.assertEqual(plain["name"], "cat_surf.vol_thickness_pbt")
        self.assertGreater(plain["central"]["n"], 0)
        self.assertGreater(plain["gm"]["n"], plain["central"]["n"])
        self.assertLess(abs(plain["gm"]["bias"]), 0.5)
        # the options reach the measure
        self.assertAlmostEqual(shifted["gm"]["bias"] - plain["gm"]["bias"], 0.4, places=4)

    def test_a_series_is_a_list_of_results(self):
        labels, infos = [tiny(1.5)[0], tiny(2.0)[0]], [tiny(1.5)[2], tiny(2.0)[2]]
        results = phantom.evaluate_phantom(labels, infos, pbt=False)
        self.assertEqual([r["name"] for r in results], ["label image 1", "label image 2"])
        self.assertEqual([r["thickness"] for r in results], [1.5, 2.0])
        self.assertEqual(results[1]["geom"], self.result["geom"])
        # without the info: one thickness per phantom and the voxel size
        again = phantom.evaluate_phantom(labels, thickness=[1.5, 2.0], vx=1.0, pbt=False)
        self.assertEqual([r["geom"] for r in again], [r["geom"] for r in results])
        self.assertIsNone(again[0]["volume"]["ref"])

    def test_the_report_names_what_was_measured(self):
        result = phantom.evaluate_phantom(self.label, self.info,
                                          thickness_map=np.full(self.label.shape, 2.1, np.float32))
        report = phantom._report(result)
        for text in ("Known thickness 2.000 mm", "geometry: WM surf -> GM/CSF surf",
                     "thickness map, all GM", "buried sulcus", "reference", "error"):
            self.assertIn(text, report)

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is not installed")
    def test_the_figure_is_written(self):
        labels, infos = [tiny(1.5)[0], tiny(2.0)[0]], [tiny(1.5)[2], tiny(2.0)[2]]
        maps = [np.full(labels[0].shape, 1.6, np.float32), None]
        with tempfile.TemporaryDirectory() as folder:
            fname = os.path.join(folder, "qc.png")
            phantom.evaluate_phantom(labels, infos, thickness_map=maps, pbt=False, fig=fname)
            with open(fname, "rb") as handle:
                self.assertEqual(handle.read(8), b"\x89PNG\r\n\x1a\n")


class TestBadParameters(unittest.TestCase):
    def assert_refused(self, message, **kwargs):
        with self.assertRaisesRegex(ValueError, message):
            phantom.thickness_phantom(**dict(TINY, **kwargs))

    def test_the_thickness_has_to_be_positive(self):
        for thickness in (0, -1.0, [2.0, 0.0], [], float("nan")):
            self.assert_refused("thickness", thickness=thickness)

    def test_an_unknown_fold_is_refused(self):
        self.assert_refused("unknown fold 'radial'", fold="radial")
        label, _, info = phantom.thickness_phantom(2.0, **dict(TINY, fold="Cartesian"))
        self.assertEqual(info["fold"], "cartesian")
        np.testing.assert_array_equal(label, tiny()[0])

    def test_the_phantom_has_to_fit_into_the_field_of_view(self):
        self.assert_refused("field of view only offers 12 mm", thickness=4.0)
        self.assert_refused("field of view", dim=(24, 24, 18))
        # the largest thickness of a series counts
        self.assert_refused("field of view", thickness=[2.0, 4.0])

    def test_an_empty_sphere_is_refused(self):
        self.assert_refused("empty", radius=0.2, amplitude=0.0)

    def test_parameters_out_of_range_are_refused(self):
        self.assert_refused("dim", dim=(24, 24))
        self.assert_refused("dim", dim=0)
        self.assert_refused("vx", vx=0.0)
        self.assert_refused("csf", csf=-1.0)
        self.assert_refused("wavelength", wavelength=0.0)
        self.assert_refused("pve_steps", pve_steps=0)
        self.assert_refused("pve_steps", pve_steps=7.5)
        self.assert_refused("intensities", intensities=(1, 2))
        self.assert_refused("folds", fold="spherical", folds=(6,))

    def test_an_even_supersampling_is_increased_with_a_warning(self):
        with self.assertWarnsRegex(UserWarning, "increased to 3"):
            label, _, info = phantom.thickness_phantom(2.0, **TINY, supersample=2)
        self.assertEqual(info["supersample"], 3)
        np.testing.assert_array_equal(label, tiny()[0])

    def test_an_offset_below_the_sampling_of_the_distance_map_is_a_warning(self):
        with self.assertWarnsRegex(UserWarning, "scalloped"):
            phantom.thickness_phantom(2.0, **TINY, wm_offset=0.3, pve_range=0.5)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _, _, info = phantom.thickness_phantom(2.0, **TINY, wm_offset=float("nan"))
        self.assertEqual(info["wm_offset"], 0.5)

    def test_the_evaluation_needs_thickness_and_voxel_size(self):
        label = np.zeros((8, 8, 8), np.float32)
        with self.assertRaisesRegex(ValueError, "vx is required"):
            phantom.evaluate_phantom(label, thickness=2.0, pbt=False)
        with self.assertRaisesRegex(ValueError, "thickness of label image is unknown"):
            phantom.evaluate_phantom(label, vx=1.0, pbt=False)
        with self.assertRaisesRegex(ValueError, "3-D"):
            phantom.evaluate_phantom(label[0], thickness=2.0, vx=1.0, pbt=False)

    def test_a_series_needs_one_entry_per_phantom(self):
        labels = [np.zeros((8, 8, 8), np.float32)] * 2
        with self.assertRaisesRegex(ValueError, "thickness has 3 entries but 2 phantoms"):
            phantom.evaluate_phantom(labels, thickness=[1, 2, 3], vx=1.0, pbt=False)
        with self.assertRaisesRegex(ValueError, "info has to be a sequence"):
            phantom.evaluate_phantom(labels, {"thickness": 2.0, "vx": 1.0}, pbt=False)

    def test_a_thickness_map_has_to_be_on_the_grid_of_the_label(self):
        label, _, info = tiny()
        with self.assertRaisesRegex(ValueError, "shape of the label image"):
            phantom.evaluate_phantom(label, info, thickness_map=np.ones((24, 24, 23)))

    def test_a_label_without_tissue_gives_no_measurement(self):
        result = phantom.evaluate_phantom(np.zeros((8, 8, 8), np.float32), thickness=2.0, vx=1.0,
                                          pbt=False)
        self.assertEqual(result["geom"]["outer"]["n"], 0)
        self.assertTrue(np.isnan(result["geom"]["outer"]["bias"]))
        self.assertTrue(np.isnan(result["geom"]["buried"]))

    def test_images_that_do_not_belong_together_are_not_written(self):
        label, t1w, info = tiny()
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "one shape"):
                phantom.write_phantom(label, t1w[:-1], info, folder)
            self.assertEqual(os.listdir(folder), [])


if __name__ == "__main__":
    unittest.main()
