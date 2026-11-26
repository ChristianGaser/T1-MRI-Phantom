"""Names and sidecars (:mod:`t1sim.bids`) and the pipeline (:mod:`t1sim.simulate`) on a small made-up brain."""

import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

import nibabel as nib
import numpy as np

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import bids, frontend, simulate  # noqa: E402


class TestNames(unittest.TestCase):
    def test_labels_are_alphanumeric(self):
        self.assertEqual(bids.label("snr2.5 Rf-20_A"), "snr2p5Rf20A")
        self.assertEqual(bids.number(40.0), "40")
        self.assertEqual(bids.number(0.75), "0p75")

    def test_the_entities_of_the_input_are_kept(self):
        self.assertEqual(bids.prefix("/d/sub-01_ses-2_acq-mprage_desc-old_T1w.nii.gz"), "sub-01_ses-2_acq-mprage")
        self.assertEqual(bids.prefix("colin27_t1_tal_hires.nii"), "colin27_t1_tal_hires")

    def test_the_resolution_entity(self):
        self.assertEqual(bids.resolution_entity([0.497, 0.4976, 0.497]), "res-0p5mm")
        self.assertEqual(bids.resolution_entity([1, 1, 1.5]), "res-1x1x1p5mm")

    def test_the_outputs_of_a_bids_dataset_go_to_its_derivatives(self):
        folder, pipeline, source = bids.output_folder("/data/study/sub-01/ses-1/anat/sub-01_ses-1_T1w.nii.gz")
        self.assertEqual(pipeline, Path("/data/study/derivatives") / f"t1sim-{bids.GENERATED_BY['Version']}")
        self.assertEqual(folder, pipeline / "sub-01" / "ses-1" / "anat")
        self.assertEqual(source, "bids::sub-01/ses-1/anat/sub-01_ses-1_T1w.nii.gz")

    def test_other_scans_get_a_derivatives_folder_next_to_theirs(self):
        folder, pipeline, source = bids.output_folder("/data/scans/colin.nii")
        self.assertEqual(folder, Path("/data/derivatives") / f"t1sim-{bids.GENERATED_BY['Version']}")
        self.assertEqual((pipeline, source), (folder, ""))
        self.assertEqual(bids.output_folder("/data/scans/colin.nii", derivatives=False)[0], Path("/data/scans"))


class TestFiles(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.folder = Path(self._tmp.name)

    def test_a_fraction_is_stored_in_8_bit_and_read_back_as_a_fraction(self):
        volume = np.linspace(0, 1, 60, dtype=np.float32).reshape(3, 4, 5)
        path = bids.save(volume, np.diag([0.5, 0.5, 0.5, 1]), self.folder / "a" / "x_probseg.nii.gz",
                         scale=1 / 255, dtype=np.uint8)
        image = nib.load(path)
        self.assertEqual(image.get_data_dtype(), np.uint8)
        np.testing.assert_allclose(image.get_fdata(), volume, atol=0.5 / 255 + 1e-6)
        np.testing.assert_allclose(image.affine, np.diag([0.5, 0.5, 0.5, 1]))

    def test_the_sidecar_takes_numpy_values(self):
        path = bids.write_sidecar(self.folder / "x_T1w.nii.gz", {"a": np.float32(1.5), "b": np.arange(3), "c": {"d": np.bool_(True)}})
        self.assertEqual(path.name, "x_T1w.json")
        self.assertEqual(json.loads(path.read_text()), {"a": 1.5, "b": [0, 1, 2], "c": {"d": True}})

    def test_the_dataset_description_is_written_once(self):
        bids.write_dataset_description(self.folder / "t1sim-x")
        text = json.loads((self.folder / "t1sim-x" / "dataset_description.json").read_text())
        self.assertEqual((text["DatasetType"], text["GeneratedBy"][0]["Name"]), ("derivative", "t1sim"))
        (self.folder / "t1sim-x" / "dataset_description.json").write_text("{}")
        bids.write_dataset_description(self.folder / "t1sim-x")
        self.assertEqual((self.folder / "t1sim-x" / "dataset_description.json").read_text(), "{}")


def _seg(shape=(40, 44, 36)):
    """What the front end returns, for a ball of WM in cortex in CSF inside a head."""
    axes = [np.arange(n) - (n - 1) / 2 for n in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    label = (np.clip(8.5 - r, 0, 1) + np.clip(11.5 - r, 0, 1) + np.clip(14.5 - r, 0, 1)).astype(np.float32)
    means = {"CSF": 100.0, "GM": 300.0, "WM": 400.0}
    corrected = np.where(r < 17, 250.0, 5.0).astype(np.float32)          # a "scalp" around the brain
    grid = np.stack(np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij"), -1)
    return {"label": label, "normalised": label.copy(), "fractions": frontend.fractions_from_label(label),
            "corrected": corrected, "mask": label > 0, "means": means,
            "grid": nib.Nifti1Image(np.zeros(shape, np.uint8), np.diag([0.5, 0.5, 0.5, 1.0])),
            "iy": grid, "template_mat": np.diag([0.5, 0.5, 0.5, 1.0])}


class TestThePipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seg = _seg()
        cls.truth = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=cls.seg, clean=False)

    def test_the_ideal_image_is_the_fractions_times_the_tissue_intensities(self):
        truth = self.truth
        pure_wm = truth["label"] == 3
        np.testing.assert_allclose(truth["ideal"][pure_wm], 400.0, rtol=1e-5)
        np.testing.assert_allclose(truth["ideal"][truth["label"] == 0], self.seg["corrected"][truth["label"] == 0])
        self.assertIs(truth["truth_fractions"], truth["fractions"])       # nothing cleaned

    def test_an_acquisition_has_the_noise_and_the_names_asked_for(self):
        acq = simulate.acquire(self.truth, snr=20, bias_percent=0)
        clean = simulate.acquire(self.truth, snr=0, bias_percent=0)
        self.assertEqual(acq["tags"], "snr20")
        self.assertEqual(clean["tags"], "")
        wm = self.truth["label"] == 3
        self.assertAlmostEqual(float((acq["image"] - clean["image"])[wm].std()), 20.0, delta=2.0)
        self.assertEqual(acq["parameters"]["NoiseSD"], 20.0)

    def test_every_file_gets_noise_of_its_own_that_is_the_same_in_every_run(self):
        a = simulate.acquire(self.truth, snr=20, bias_percent=0)
        b = simulate.acquire(self.truth, snr=20, bias_percent=0)
        other = simulate.ground_truth("/data/sub-08/anat/sub-08_T1w.nii.gz", seg=self.seg, clean=False)
        c = simulate.acquire(other, snr=20, bias_percent=0)
        np.testing.assert_array_equal(a["image"], b["image"])
        self.assertFalse(np.array_equal(a["image"], c["image"]))
        fixed = simulate.acquire(other, snr=20, bias_percent=0, seed=5)
        np.testing.assert_array_equal(fixed["image"], simulate.acquire(self.truth, snr=20, bias_percent=0, seed=5)["image"])
        self.assertTrue(fixed["tags"].endswith("Seed5"))

    def test_the_tags_of_the_acquisition(self):
        acq = simulate.acquire(self.truth, snr=30, bias_percent=-20, bias_type=(3, 1), contrast=0.5)
        self.assertEqual(acq["tags"], "snr30RfNeg20T3ConLow")
        self.assertEqual(acq["parameters"]["BiasFieldType"], "T3")
        acq = simulate.acquire(self.truth, snr=0, noise_percent=3, bias_percent=0, contrast=1.3)
        self.assertEqual(acq["tags"], "pn3Con1p3")

    def test_the_snr_refers_to_white_matter_after_a_contrast_change(self):
        acq = simulate.acquire(self.truth, snr=25, bias_percent=0, contrast=0.5)
        clean = simulate.acquire(self.truth, snr=0, bias_percent=0, contrast=0.5)
        wm = float(np.median(clean["image"][self.truth["label"] == 3]))
        self.assertAlmostEqual(acq["parameters"]["NoiseSD"], wm / 25, delta=0.02)

    def test_another_voxel_size_averages_the_ground_truth(self):
        acq = simulate.acquire(self.truth, snr=0, bias_percent=0, resolution=1.0)
        self.assertEqual(acq["image"].shape, (20, 22, 18))
        self.assertEqual(acq["label"].shape, (20, 22, 18))
        np.testing.assert_allclose(np.sqrt((acq["affine"][:3, :3] ** 2).sum(0)), 1.0)
        fine = self.truth["truth_fractions"].sum((0, 1, 2), dtype=np.float64) * 0.125
        np.testing.assert_allclose(acq["fractions"].sum((0, 1, 2), dtype=np.float64), fine, rtol=1e-4)

    def test_back_on_the_grid_the_ground_truth_keeps_its_voxels(self):
        acq = simulate.acquire(self.truth, snr=0, bias_percent=0, resolution=1.0, upsample=True)
        self.assertEqual(acq["image"].shape, self.truth["ideal"].shape)
        self.assertIs(acq["label"], self.truth["label"])
        self.assertEqual(acq["tags"], "Acq1mm")
        sharp = simulate.acquire(self.truth, snr=0, bias_percent=0)
        edge = lambda image: float(np.abs(np.diff(image, axis=0)).max())
        self.assertLess(edge(acq["image"]), edge(sharp["image"]))          # the resolution is that of 1 mm

    def test_the_files_pair_the_way_the_training_scripts_pair_them(self):
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(self.truth, simulate.acquire(self.truth, snr=40), out_dir=folder)
            names = sorted(p.name for p in Path(folder).iterdir())
            base = "sub-07_space-MNI152NLin2009cAsym_res-0p5mm"
            self.assertEqual(files["T1w"].name, f"{base}_desc-snr40Rf30T2_T1w.nii.gz")
            self.assertEqual(files["dseg"].name, f"{base}_dseg.nii.gz")
            self.assertIn(f"{base}_label-GM_probseg.nii.gz", names)
            self.assertIn(f"{base}_desc-snr40Rf30T2_T1w.json", names)
            label = nib.load(files["dseg"]).get_fdata()
            np.testing.assert_allclose(label, self.truth["label"], atol=1 / 765)
            sidecar = json.loads(Path(str(files["T1w"]).replace(".nii.gz", ".json")).read_text())
            self.assertEqual(sidecar["SimulationParameters"]["SNR"], 40.0)
            self.assertEqual(sidecar["GeneratedBy"][0]["Name"], "t1sim")
            # deepmriprep-train: base before desc, the label's tags end the image's
            split = lambda n: re.match(r"^(.*?)(?:_desc-([A-Za-z0-9]+))?_([A-Za-z0-9-]+)$", n[:-7]).groups()
            image, truth = split(files["T1w"].name), split(files["dseg"].name)
            self.assertEqual(image[0], truth[0])
            self.assertTrue(image[1].endswith(truth[1] or ""))

    def test_the_deepmriprep_grid_is_cut_from_the_working_grid_only(self):
        with self.assertRaisesRegex(ValueError, "339 x 411 x 339"):
            simulate.acquire(self.truth, grid="deepmriprep")
        with self.assertRaises(ValueError):
            simulate.acquire(self.truth, grid="other")

    def test_thickness_and_atrophy_exclude_each_other_and_options_are_checked(self):
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            simulate.ground_truth("x", seg=self.seg, thickness=2.5, atrophy={"regions": [28], "change": -0.1})
        with self.assertRaisesRegex(ValueError, "unknown cleanup"):
            simulate.ground_truth("x", seg=self.seg, clean={"vessel": 1})
        with self.assertRaisesRegex(ValueError, "regions"):
            simulate.ground_truth("x", seg=self.seg, clean=False, atrophy={"change": -0.1})

    def test_the_crop_of_the_deepmriprep_grid(self):
        shape = tuple(len(range(*c.indices(n))) for c, n in zip(simulate.DEEPMRIPREP_CROP, (339, 411, 339)))
        self.assertEqual(shape, (336, 384, 336))
        moved = simulate._crop_affine(np.diag([0.5, 0.5, 0.5, 1.0]), simulate.DEEPMRIPREP_CROP)
        np.testing.assert_allclose(moved[:3, 3], [0.5, 7.5, 0.0])


if __name__ == "__main__":
    unittest.main()
