"""Names and sidecars (:mod:`t1sim.bids`) and the pipeline (:mod:`t1sim.simulate`) on a small made-up brain."""

import json
import re
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

from t1sim import artefacts, bids, frontend, simulate  # noqa: E402

try:
    import torch  # noqa: F401
    HAVE_TORCH = True
except ImportError:
    HAVE_TORCH = False


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
        self.assertEqual(pipeline, Path("/data/derivatives") / f"t1sim-{bids.GENERATED_BY['Version']}")
        self.assertEqual((folder, source), (pipeline / "scans", ""))
        self.assertEqual(bids.output_folder("/data/scans/colin.nii", derivatives=False)[0], Path("/data/scans"))

    def test_scans_of_one_name_in_two_folders_do_not_share_their_outputs(self):
        first = bids.output_folder("/data/subj01/T1.nii")[0]
        second = bids.output_folder("/data/subj02/T1.nii")[0]
        self.assertNotEqual(first, second)
        self.assertEqual(first.parent, second.parent)

    def test_the_thickness_in_a_name(self):
        self.assertEqual(bids.thickness_tag(2.5), "Thickness25mm")
        self.assertEqual(bids.thickness_tag(np.array([0.8])), "Thickness08mm")
        self.assertEqual(bids.thickness_tag(1.25), "Thickness13mm")       # rounded as MATLAB rounds
        self.assertEqual(bids.thickness_tag([2.5, 1.5, 2.0]), "Thickness15to25mm")


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
    fractions = frontend.fractions_from_label(label)
    # as in the front end, the label map is that of the fractions: at the rim of the brain a voxel is CSF or nothing
    return {"label": frontend.label_from_fractions(fractions), "normalised": label.copy(), "fractions": fractions,
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
        self.assertIs(truth["image_fractions"], truth["fractions"])       # nothing cleaned

    def test_an_acquisition_has_the_noise_and_the_names_asked_for(self):
        acq = simulate.acquire(self.truth, snr=20, bias_percent=0, grid="working")
        clean = simulate.acquire(self.truth, snr=0, bias_percent=0, grid="working")
        self.assertEqual(acq["tags"], "snr20")
        self.assertEqual(clean["tags"], "")
        wm = self.truth["label"] == 3
        self.assertAlmostEqual(float((acq["image"] - clean["image"])[wm].std()), 20.0, delta=2.0)
        self.assertEqual(acq["parameters"]["NoiseSD"], 20.0)

    def test_every_file_gets_noise_of_its_own_that_is_the_same_in_every_run(self):
        a = simulate.acquire(self.truth, snr=20, bias_percent=0, grid="working")
        b = simulate.acquire(self.truth, snr=20, bias_percent=0, grid="working")
        other = simulate.ground_truth("/data/sub-08/anat/sub-08_T1w.nii.gz", seg=self.seg, clean=False)
        c = simulate.acquire(other, snr=20, bias_percent=0, grid="working")
        np.testing.assert_array_equal(a["image"], b["image"])
        self.assertFalse(np.array_equal(a["image"], c["image"]))
        fixed = simulate.acquire(other, snr=20, bias_percent=0, seed=5, grid="working")
        np.testing.assert_array_equal(fixed["image"], simulate.acquire(self.truth, snr=20, bias_percent=0, seed=5, grid="working")["image"])
        self.assertTrue(fixed["tags"].endswith("Seed5"))
        # a seed of 0 is a seed that was given: its file is not that of the default
        self.assertEqual(simulate.acquire(self.truth, snr=20, bias_percent=0, seed=0, grid="working")["tags"], "snr20Seed0")

    def test_scans_of_one_name_in_two_folders_get_noise_and_a_field_of_their_own(self):
        first = simulate.ground_truth("/data/subj01/T1.nii", seg=self.seg, clean=False)
        second = simulate.ground_truth("/data/subj02/T1.nii", seg=self.seg, clean=False)
        a = simulate.acquire(first, snr=20, bias_type=(2, None), grid="working")
        b = simulate.acquire(second, snr=20, bias_type=(2, None), grid="working")
        self.assertNotEqual(a["parameters"]["Seed"], b["parameters"]["Seed"])
        self.assertNotEqual(a["parameters"]["BiasFieldSeed"], b["parameters"]["BiasFieldSeed"])
        self.assertFalse(np.array_equal(a["bias"], b["bias"]))

    def test_the_ground_truth_is_not_changed_by_an_acquisition(self):
        before = {key: self.truth[key].copy() for key in ("fractions", "label", "ideal", "corrected")}
        simulate.acquire(self.truth, snr=20, contrast=0.5, resolution=1.0, grid="working")
        for key, value in before.items():
            np.testing.assert_array_equal(self.truth[key], value)

    def test_the_tags_of_the_acquisition(self):
        acq = simulate.acquire(self.truth, snr=30, bias_percent=-20, bias_type=(3, 0), contrast=0.5, grid="working")
        self.assertEqual(acq["tags"], "snr30RfNeg20T3ConLow")
        self.assertEqual(acq["parameters"]["BiasFieldType"], "T3")
        # another seed of the field is another file
        other = simulate.acquire(self.truth, snr=30, bias_percent=-20, bias_type=(3, 7), contrast=0.5, grid="working")
        self.assertEqual(other["tags"], "snr30RfNeg20T3S7ConLow")
        self.assertEqual((other["parameters"]["BiasFieldType"], other["parameters"]["BiasFieldSeed"]), ("T3", 7))
        acq = simulate.acquire(self.truth, snr=0, noise_percent=3, bias_percent=0, contrast=1.3, grid="working")
        self.assertEqual(acq["tags"], "pn3Con1p3")
        with self.assertRaisesRegex(ValueError, "positive"):
            simulate.acquire(self.truth, contrast=0, grid="working")

    def test_the_snr_refers_to_white_matter_after_a_contrast_change(self):
        acq = simulate.acquire(self.truth, snr=25, bias_percent=0, contrast=0.5, grid="working")
        clean = simulate.acquire(self.truth, snr=0, bias_percent=0, contrast=0.5, grid="working")
        wm = float(np.median(clean["image"][self.truth["label"] == 3]))
        self.assertAlmostEqual(wm, 400.0, delta=0.01)                   # white matter keeps its intensity
        self.assertAlmostEqual(acq["parameters"]["NoiseSD"], wm / 25, delta=0.02)

    def test_a_contrast_change_keeps_a_mixed_voxel_the_mixture_of_its_tissues(self):
        truth = self.truth
        brain = truth["fractions"].sum(-1) > 0.999
        for exponent in (0.5, 1.5):
            with self.subTest(exponent=exponent):
                image = simulate.acquire(truth, snr=0, bias_percent=0, contrast=exponent, grid="working")
                means = [artefacts.change_contrast(truth["means"][t], exponent, 400.0) for t in frontend.TISSUES]
                np.testing.assert_allclose(image["image"][brain], (truth["fractions"] @ np.float32(means))[brain],
                                           rtol=1e-4)
                self.assertEqual(image["parameters"]["TissueIntensities"]["GM"], round(means[1], 3))
                # outside the brain the scan is mapped: 250 of the "scalp"
                scalp = (truth["label"] == 0) & (truth["corrected"] == 250.0)
                np.testing.assert_allclose(image["image"][scalp],
                                           artefacts.change_contrast(250.0, exponent, 400.0), rtol=1e-4)

    def test_a_contrast_change_leaves_the_bias_field_what_was_asked_for(self):
        plain = simulate.acquire(self.truth, snr=0, bias_percent=0, contrast=0.5, grid="working")["image"]
        biased = simulate.acquire(self.truth, snr=0, bias_percent=30, contrast=0.5, grid="working")
        np.testing.assert_allclose(biased["image"], plain * biased["bias"], rtol=1e-4, atol=1e-3)

    def test_another_voxel_size_averages_the_ground_truth(self):
        acq = simulate.acquire(self.truth, snr=0, bias_percent=0, resolution=1.0, grid="working")
        self.assertEqual(acq["image"].shape, (20, 22, 18))
        self.assertEqual(acq["label"].shape, (20, 22, 18))
        np.testing.assert_allclose(np.sqrt((acq["affine"][:3, :3] ** 2).sum(0)), 1.0)
        fine = self.truth["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.125
        np.testing.assert_allclose(acq["fractions"].sum((0, 1, 2), dtype=np.float64), fine, rtol=1e-4)

    def test_back_on_the_grid_the_ground_truth_keeps_its_voxels(self):
        acq = simulate.acquire(self.truth, snr=0, bias_percent=0, resolution=1.0, upsample=True, grid="working")
        self.assertEqual(acq["image"].shape, self.truth["ideal"].shape)
        self.assertIs(acq["label"], self.truth["label"])
        self.assertEqual(acq["tags"], "Acq1mm")
        sharp = simulate.acquire(self.truth, snr=0, bias_percent=0, grid="working")
        edge = lambda image: float(np.abs(np.diff(image, axis=0)).max())
        self.assertLess(edge(acq["image"]), edge(sharp["image"]))          # the resolution is that of 1 mm

    def test_the_files_pair_the_way_the_training_scripts_pair_them(self):
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(self.truth, simulate.acquire(self.truth, snr=40, grid="working"), out_dir=folder)
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

    def test_a_ground_truth_of_other_settings_is_not_replaced(self):
        acquired = simulate.acquire(self.truth, snr=40, grid="working")
        with tempfile.TemporaryDirectory() as folder:
            simulate.write(self.truth, acquired, out_dir=folder)
            simulate.write(self.truth, simulate.acquire(self.truth, snr=20, grid="working"), out_dir=folder)   # the same one
            # the same names, but made with a thickness: the images in the folder would lose their ground truth
            other = dict(self.truth, settings={**self.truth["settings"], "Thickness": [2.5]})
            with self.assertRaisesRegex(FileExistsError, "Thickness"):
                simulate.write(other, acquired, out_dir=folder)
            # and so they would with a ground truth on another grid
            coarse = simulate.acquire(self.truth, snr=40, resolution=1.0, grid="working")
            coarse["affine"] = acquired["affine"]                   # as if the name were the same
            with self.assertRaisesRegex(FileExistsError, "Shape"):
                simulate.write(self.truth, coarse, out_dir=folder)
            # a ground truth that was written without the record of its settings is checked by its grid
            dseg = next(Path(folder).glob("*_dseg.json"))
            dseg.write_text("{}")
            simulate.write(other, acquired, out_dir=folder)
            dseg.write_text("{}")
            with self.assertRaisesRegex(FileExistsError, "another grid"):
                simulate.write(self.truth, coarse, out_dir=folder)

    def test_the_bias_field_is_written_under_the_name_of_its_image(self):
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(self.truth, simulate.acquire(self.truth, snr=40, grid="working"), out_dir=folder, save_bias=True)
            self.assertEqual(files["bias"].name,
                             "sub-07_space-MNI152NLin2009cAsym_res-0p5mm_desc-snr40Rf30T2Biasfield_T1w.nii.gz")

    def test_the_bias_type_of_the_command_line(self):
        self.assertEqual(simulate._bias_type_argument(["2", "0"]), (2.0, 0))
        self.assertEqual(simulate._bias_type_argument(["3"]), (3.0, 0))
        self.assertEqual(simulate._bias_type_argument(["2", "scan"]), (2.0, None))
        self.assertEqual(simulate._bias_type_argument(["A"]), "A")
        self.assertEqual(simulate._bias_type_argument(["/fields/rf.nii"]), "/fields/rf.nii")
        with self.assertRaises(SystemExit):
            simulate._bias_type_argument(["A", "0"])

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


def _seg_with_scan():
    """The made-up brain of :func:`_seg`, with the scan it would come from.

    The scan has voxels of 0.7 mm, is stored with its first axis running to
    the left, is turned against the working grid, shows the head 10% smaller
    (the registration to the template enlarged it) and reaches further than
    the working grid.  It is half as bright as the bias-corrected scan.
    """
    from scipy import ndimage

    seg = _seg()
    seg["corrected"] = seg["corrected"].copy()
    seg["corrected"][12:28, 14:30, :8] = 250.0          # a neck that leaves the working grid at its lower end
    a, b = np.deg2rad(8.0), np.deg2rad(-5.0)
    turn = (np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
            @ np.array([[np.cos(b), -np.sin(b), 0], [np.sin(b), np.cos(b), 0], [0, 0, 1]]))
    world = np.eye(4)                                   # millimetres of the scan -> of the template
    world[:3, :3] = 1.1 * turn
    world[:3, 3] = [0.6, -0.4, 0.3]
    shape = (40, 42, 38)
    centre = np.linalg.solve(world, seg["grid"].affine @ np.append((np.array(seg["label"].shape) - 1) / 2, 1))[:3]
    affine = np.diag([-0.7, 0.7, 0.7, 1.0])
    affine[:3, 3] = centre - affine[:3, :3] @ ((np.array(shape) - 1) / 2)
    to_grid = np.linalg.solve(seg["grid"].affine, world @ affine)
    data = 0.5 * ndimage.affine_transform(seg["corrected"], to_grid[:3, :3], offset=to_grid[:3, 3],
                                          output_shape=shape, order=1, mode="nearest")
    seg["scan"] = nib.Nifti1Image(data.astype(np.float32), affine)
    seg["scan"].header.set_sform(affine, code=1)
    seg["world"] = world
    return seg


class TestTheGridOfTheScan(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seg = _seg_with_scan()
        cls.scan = cls.seg["scan"]
        cls.truth = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=cls.seg, clean=False)
        cls.plain = simulate.acquire(cls.truth, grid="scan", snr=0, bias_percent=0)

    def test_the_image_and_its_ground_truth_lie_on_the_scan(self):
        acq = self.plain
        self.assertEqual(acq["image"].shape, self.scan.shape)
        self.assertEqual(acq["fractions"].shape, self.scan.shape + (3,))
        np.testing.assert_allclose(acq["affine"], self.scan.affine)
        np.testing.assert_allclose(acq["parameters"]["VoxelSize"], [0.7, 0.7, 0.7], atol=1e-4)
        self.assertEqual((acq["parameters"]["Grid"], acq["parameters"]["Space"]), ("scan", "scan"))
        np.testing.assert_allclose(acq["label"], frontend.label_from_fractions(acq["fractions"]), atol=1e-6)
        # the same object as the scan, voxel by voxel: bright where the scan is bright
        scan = np.asarray(self.scan.dataobj)
        self.assertGreater(float(np.corrcoef(acq["image"].ravel(), scan.ravel())[0, 1]), 0.75)

    def test_the_head_has_the_size_of_the_subject(self):
        on_grid = self.truth["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.5 ** 3
        on_scan = self.plain["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.7 ** 3
        # the registration enlarged the head by 1.1 along every axis
        np.testing.assert_allclose(on_scan, on_grid / 1.1 ** 3, rtol=2e-3)
        self.assertLessEqual(float(self.plain["fractions"].sum(-1).max()), 1 + 1e-5)

    def test_pure_white_matter_has_its_intensity(self):
        white = self.plain["fractions"][..., 2] > 0.9999
        self.assertGreater(int(white.sum()), 200)
        self.assertAlmostEqual(float(self.plain["image"][white].mean()), 400.0, delta=4.0)

    def test_beyond_the_working_grid_the_image_is_the_scan(self):
        from t1sim import artefacts

        truth = self.truth
        lattice = artefacts.ObliqueAcquisition(truth["ideal"].shape, truth["affine"], self.scan.shape,
                                               truth["world"] @ self.scan.affine)
        beyond = lattice.coverage() == 0
        self.assertGreater(float(beyond.mean()), 0.2)
        scan = np.asarray(self.scan.dataobj)
        bright = beyond & (scan > 50)
        self.assertGreater(int(bright.sum()), 500)
        # brought to the intensities of the corrected scan, which is twice as bright here
        np.testing.assert_allclose(self.plain["image"][bright], 2 * scan[bright], rtol=0.05)
        # and with the bias field of the simulation, continued beyond the grid it was made on
        biased = simulate.acquire(truth, grid="scan", snr=0, bias_percent=30)
        self.assertEqual(biased["bias"].shape, self.scan.shape)
        np.testing.assert_allclose(biased["image"][bright], (self.plain["image"] * biased["bias"])[bright], rtol=1e-3)
        self.assertGreater(float(biased["bias"][beyond].min()), 0.7)

    def test_the_noise_is_white_on_the_scan_and_the_ground_truth_is_made_once(self):
        noisy = simulate.acquire(self.truth, grid="scan", snr=20, bias_percent=0)
        white = self.plain["fractions"][..., 2] > 0.9999
        self.assertAlmostEqual(float((noisy["image"] - self.plain["image"])[white].std()), 20.0, delta=2.0)
        again = simulate.acquire(self.truth, grid="scan", snr=40, bias_percent=0)
        self.assertIs(again["fractions"], noisy["fractions"])
        self.assertIs(again["label"], noisy["label"])
        np.testing.assert_array_equal(noisy["fractions"], self.plain["fractions"])

    def test_the_jacobian_of_an_atrophy_comes_along_and_is_one_where_nothing_was_deformed(self):
        # a copy of the ground truth with a Jacobian: what is kept for the scan must not be that of the original
        shrunk = dict(self.truth, jacobian=np.full(self.truth["label"].shape, 0.9, np.float32))
        acq = simulate.acquire(shrunk, grid="scan", snr=0, bias_percent=0)
        self.assertEqual(acq["jacobian"].shape, self.scan.shape)
        brain = acq["fractions"].sum(-1) > 0.999
        np.testing.assert_allclose(acq["jacobian"][brain], 0.9, atol=1e-4)
        np.testing.assert_allclose(acq["jacobian"][0, 0, 0], 1.0, atol=1e-6)        # beyond the working grid
        self.assertIsNone(simulate.acquire(self.truth, grid="scan", snr=0, bias_percent=0)["jacobian"])

    def test_another_voxel_size_stays_on_the_axes_of_the_scan(self):
        coarse = simulate.acquire(self.truth, grid="scan", snr=0, bias_percent=0, resolution=1.4)
        self.assertEqual(coarse["image"].shape, (20, 21, 19))
        np.testing.assert_allclose(coarse["parameters"]["VoxelSize"], [1.4, 1.4, 1.4], atol=1e-4)
        np.testing.assert_allclose(coarse["affine"][:3, :3], 2 * self.scan.affine[:3, :3], atol=1e-6)
        fine = self.plain["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.7 ** 3
        np.testing.assert_allclose(coarse["fractions"].sum((0, 1, 2), dtype=np.float64) * 1.4 ** 3, fine, rtol=5e-3)

    def test_the_names_carry_no_template_space(self):
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(self.truth, simulate.acquire(self.truth, grid="scan", snr=40), out_dir=folder)
            self.assertEqual(files["T1w"].name, "sub-07_desc-snr40Rf30T2_T1w.nii.gz")
            self.assertEqual(files["dseg"].name, "sub-07_dseg.nii.gz")
            self.assertEqual(files["probseg"]["GM"].name, "sub-07_label-GM_probseg.nii.gz")
            image = nib.load(files["T1w"])
            self.assertEqual(image.shape, self.scan.shape)
            np.testing.assert_allclose(image.affine, self.scan.affine, atol=1e-5)
            self.assertEqual(int(image.header["sform_code"]), 1)           # the scanner, not the template
            coarse = simulate.acquire(self.truth, grid="scan", snr=40, resolution=1.4)
            self.assertEqual(simulate.write(self.truth, coarse, out_dir=folder)["T1w"].name,
                             "sub-07_res-1p4mm_desc-snr40Rf30T2_T1w.nii.gz")
            # in template space the code stays that of the template
            files = simulate.write(self.truth, simulate.acquire(self.truth, snr=40, grid="working"), out_dir=folder)
            self.assertEqual(int(nib.load(files["T1w"]).header["sform_code"]), 4)
            self.assertIn("_space-MNI152NLin2009cAsym_res-0p5mm_", files["T1w"].name)

    def test_the_grid_of_the_scan_is_the_default(self):
        acq = simulate.acquire(self.truth, snr=0, bias_percent=0)
        self.assertEqual(acq["parameters"]["Grid"], "scan")
        np.testing.assert_array_equal(acq["image"], self.plain["image"])

    def test_nothing_is_below_zero_in_what_is_written(self):
        # the cut of k-space and the splines swing below 0 next to the bright scalp, and so does ringing
        for options in (dict(), dict(ringing=2), dict(resolution=1.4), dict(grid="working", resolution=1.0, upsample=True)):
            with self.subTest(options=options):
                acq = simulate.acquire(self.truth, **{"snr": 0, "bias_percent": 30, **options})
                self.assertGreaterEqual(float(acq["image"].min()), 0.0)
                self.assertGreaterEqual(float(acq["fractions"].min()), 0.0)
                self.assertLessEqual(float(acq["fractions"].sum(-1).max()), 1 + 1e-5)

    def test_what_cannot_be_done_on_the_scan_is_refused(self):
        with self.assertRaisesRegex(ValueError, "upsample"):
            simulate.acquire(self.truth, resolution=1.4, upsample=True)
        without = simulate.ground_truth("x", seg=_seg(), clean=False)
        with self.assertRaisesRegex(ValueError, "grid='working'"):
            simulate.acquire(without)


class TestTheCleanedGroundTruth(unittest.TestCase):
    """The pipeline with a cleanup: what is detected is given, how it is used is tested."""

    def setUp(self):
        self.seg = _seg()
        label = self.seg["label"]
        cortex = np.argwhere(label == 2)
        self.vessel = np.zeros(label.shape, bool)
        self.vessel[tuple(cortex[len(cortex) // 2])] = True
        self.rim = np.zeros(label.shape, bool)
        self.rim[tuple(cortex[len(cortex) // 3])] = True
        masks = {"vessels": self.vessel, "dura": np.zeros(label.shape, bool), "pve": self.rim}
        for target, value in (("atlas.to_grid", np.full(label.shape, 100, np.uint8)),      # one cortical parcel
                              ("cleanup.detect", masks)):
            patcher = mock.patch(f"t1sim.simulate.{target}", return_value=value if target == "atlas.to_grid" else None,
                                 side_effect=None if target == "atlas.to_grid" else self.detect)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.masks = masks

    def detect(self, *args, **kwargs):
        return {name: mask.copy() for name, mask in self.masks.items()}

    def test_the_image_shows_what_the_ground_truth_rejects(self):
        truth = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=self.seg)
        self.assertIsNot(truth["fractions"], truth["image_fractions"])
        # the vessel is CSF in the ground truth and has the intensity of GM in the image
        np.testing.assert_array_equal(truth["fractions"][self.vessel], [[1, 0, 0]])
        self.assertEqual(float(truth["label"][self.vessel][0]), 1.0)
        self.assertAlmostEqual(float(truth["ideal"][self.vessel][0]), 300.0, places=3)
        # the rim is CSF and WM with the intensity of the image: (400 - 300) / (400 - 100) of CSF
        np.testing.assert_allclose(truth["fractions"][self.rim], [[1 / 3, 0, 2 / 3]], atol=1e-6)
        self.assertAlmostEqual(float(truth["label"][self.rim][0]), 3 - 2 / 3, places=5)
        np.testing.assert_allclose(truth["label"], frontend.label_from_fractions(truth["fractions"]), atol=1e-6)
        self.assertEqual(truth["parameters"]["LabelCleanup"]["VolumeMl"]["vessels"], round(0.125 / 1000, 2))
        self.assertTrue(truth["cleaned"])

    def test_the_files_pair_as_deepmriprep_train_pairs_them(self):
        truth = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=self.seg)
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(truth, simulate.acquire(truth, snr=40, grid="working"), out_dir=folder)
            self.assertEqual(files["dseg"].name, "sub-07_space-MNI152NLin2009cAsym_res-0p5mm_desc-Clean_dseg.nii.gz")
            # pair_simulations: the desc of the label without "Clean" ends the desc of the image
            image_desc = re.search(r"_desc-([A-Za-z0-9]+)_T1w", files["T1w"].name).group(1)
            label_desc = re.search(r"_desc-([A-Za-z0-9]+)_dseg", files["dseg"].name).group(1)
            self.assertTrue(label_desc.endswith("Clean"))
            self.assertTrue(image_desc.endswith(label_desc[:-5]))
            self.assertEqual(files["T1w"].name.split("_desc-")[0], files["dseg"].name.split("_desc-")[0])
            np.testing.assert_allclose(nib.load(files["probseg"]["GM"]).get_fdata(), truth["fractions"][..., 1],
                                       atol=0.5 / 255 + 1e-6)
            sidecar = json.loads(Path(str(files["dseg"]).replace(".nii.gz", ".json")).read_text())
            self.assertEqual(sidecar["GroundTruthSettings"]["Cleanup"], {"vessels": 0.5, "dura": 4.0, "pve": True})
            self.assertEqual(sidecar["GroundTruthSettings"]["Grid"], "working")

    def test_a_rebuilt_cortex_is_not_cleaned(self):
        # everything the detection would take lies in the cortex a thickness rebuilds
        self.masks["dura"] = self.seg["label"] == 2
        truth = simulate.ground_truth("x", seg=self.seg, thickness=2.5)
        # 2.5 mm are 5 voxels of cortex along a ray, in the ground truth as in the image
        centre = tuple(n // 2 for n in truth["label"].shape)
        for name in ("fractions", "image_fractions"):
            along = truth[name][centre[0]:, centre[1], centre[2], 1]
            self.assertAlmostEqual(float(along.sum()) * 0.5, 2.5, delta=0.3, msg=name)
        self.assertEqual(truth["tags"], "Thickness25mm")
        self.assertEqual(truth["parameters"]["LabelCleanup"]["VolumeMl"]["dura"], 0.0)

    def test_a_structure_that_a_thickness_keeps_is_still_cleaned(self):
        # one half of the brain is no cortex: it keeps its labels, and with them what the cleanup finds there
        shape = self.seg["label"].shape
        half = np.arange(shape[0])[:, None, None] * np.ones(shape) > shape[0] / 2 + 3
        regions = np.where(half, 36, 100).astype(np.uint8)
        self.masks["dura"] = self.seg["label"] == 2
        with mock.patch("t1sim.simulate.atlas.to_grid", return_value=regions):
            truth = simulate.ground_truth("x", seg=self.seg, thickness=2.5)
        centre = tuple(n // 2 for n in shape)
        kept = (slice(centre[0] + 6, None), centre[1], centre[2], 1)
        rebuilt = truth["fractions"][:centre[0], centre[1], centre[2], 1]
        # the image keeps the 1.5 mm of cortex there; in the ground truth all of it that was pure became CSF
        self.assertAlmostEqual(float(truth["image_fractions"][kept].sum()) * 0.5, 1.5, delta=0.2)
        self.assertLess(float(truth["fractions"][kept].sum()) * 0.5, 0.7)
        self.assertAlmostEqual(float(rebuilt.sum()) * 0.5, 2.5, delta=0.3)


@unittest.skipUnless(HAVE_TORCH, "torch is needed for the deformations")
class TestAtrophyInThePipeline(unittest.TestCase):
    """Two neighbouring regions of an atlas that is given as a file."""

    @classmethod
    def setUpClass(cls):
        cls.seg = _seg()
        shape = cls.seg["label"].shape
        x = np.arange(shape[0])[:, None, None] * np.ones(shape)
        regions = np.where(cls.seg["label"] > 0, np.where(x < shape[0] / 2, 1, 2), 0).astype(np.uint8)
        cls._tmp = tempfile.TemporaryDirectory()
        cls.atlas = str(Path(cls._tmp.name) / "halves.nii.gz")
        nib.save(nib.Nifti1Image(regions, cls.seg["template_mat"]), cls.atlas)
        cls.regions = regions
        cls.both = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=cls.seg, clean=False,
                                         atrophy={"atlas": cls.atlas, "regions": [1, 2], "change": -0.1})
        cls.one = simulate.ground_truth("/data/sub-07/anat/sub-07_T1w.nii.gz", seg=cls.seg, clean=False,
                                        atrophy={"atlas": cls.atlas, "regions": [1], "change": -0.1})

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_every_region_reaches_what_was_asked_for(self):
        done = self.both["parameters"]["Atrophy"]
        self.assertEqual([step["Region"] for step in done], [1, 2])
        for step in done:
            self.assertAlmostEqual(step["Achieved"], -0.1, delta=0.005)
            self.assertLess(step["ChangeInRegion"], -0.05)
            self.assertLess(step["TissueChangeMl"]["GM"], 0)
            self.assertGreater(step["TissueChangeMl"]["CSF"], 0)

    def test_the_sidecar_holds_the_change_of_the_fractions_that_are_written(self):
        before = self.seg["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.125 / 1000
        after = self.both["fractions"].sum((0, 1, 2), dtype=np.float64) * 0.125 / 1000
        steps = self.both["parameters"]["Atrophy"]
        for index, tissue in enumerate(frontend.TISSUES):
            total = sum(step["TissueChangeMl"][tissue] for step in steps)
            self.assertAlmostEqual(total, after[index] - before[index], delta=0.002)

    def test_the_names(self):
        self.assertEqual(self.one["tags"], "halvesRoi1Thin10")
        self.assertEqual(self.both["tags"], "halvesMulti")
        thicker = simulate.ground_truth("x", seg=self.seg, clean=False,
                                        atrophy={"atlas": self.atlas, "regions": [1], "change": 0.1})
        self.assertEqual(thicker["tags"], "halvesRoi1ThinPlus10")       # not the name of a loss of 10%

    def test_the_jacobian_of_two_regions_is_chained(self):
        jacobian = self.both["jacobian"]
        self.assertEqual(jacobian.shape, self.seg["label"].shape)
        first = self.one["jacobian"]
        # far from the second region the first deformation is all there is
        far = np.zeros(jacobian.shape, bool)
        far[:6] = True
        np.testing.assert_allclose(jacobian[far], first[far], atol=0.02)
        # the grey matter that was there, times the Jacobian, is the grey matter that is left
        gm = self.seg["fractions"][..., 1].astype(np.float64)
        left = float(self.both["fractions"][..., 1].sum(dtype=np.float64))
        self.assertAlmostEqual(float((gm * jacobian).sum()) / left, 1.0, delta=0.03)

    def test_the_jacobian_is_written_with_the_space_it_lies_in(self):
        acquired = simulate.acquire(self.one, snr=40, grid="working")
        with tempfile.TemporaryDirectory() as folder:
            files = simulate.write(self.one, acquired, out_dir=folder)
            self.assertTrue(files["jacobian"].name.endswith("_desc-halvesRoi1Thin10_jacobian.nii.gz"))
            np.testing.assert_allclose(nib.load(files["jacobian"]).get_fdata(), self.one["jacobian"], atol=1e-6)
            sidecar = json.loads(Path(str(files["jacobian"]).replace(".nii.gz", ".json")).read_text())
            self.assertIn("BEFORE the atrophy", sidecar["Description"])
            self.assertEqual(sidecar["SimulationParameters"]["Atrophy"][0]["Region"], 1)
            self.assertNotIn("jacobian", simulate.write(self.one, acquired, out_dir=folder, save_jacobian=False))


if __name__ == "__main__":
    unittest.main()
