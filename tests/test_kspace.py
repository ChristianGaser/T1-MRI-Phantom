"""The k-space artefacts (:mod:`t1sim.kspace`): movement and ringing.

Properties of the model, on small volumes: what the options resolve to, how
the pose time course is built, what a translation, a rotation and the two
ringing types do to an image, and that arrays and tensors come back as what
they were.  That the images are those of ``matlab/mri_simulate.m`` was
checked against MATLAB itself; of that comparison a handful of values of
SPM's sinc interpolation are kept here.
"""

import json
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import torch

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import kspace  # noqa: E402

IDENTITY = np.eye(4)


def oblique_affine():
    """Sagittal storage (voxel axes along the world's y, z and x), anisotropic, tilted by 10 degrees."""
    c, s = np.cos(np.deg2rad(10)), np.sin(np.deg2rad(10))
    tilt = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    permute = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    affine = np.eye(4)
    affine[:3, :3] = tilt @ permute @ np.diag([0.9, 1.2, 1.1])
    affine[:3, 3] = [-12, -14, -9]
    return affine


def phantom(shape=(24, 28, 20), dtype=np.float32):
    """A smooth object with sharp edges: an ellipsoid with a bright core and a dark box."""
    x, y, z = np.meshgrid(*[np.arange(n) - (n - 1) / 2 for n in shape], indexing="ij")
    r2 = (x / (0.36 * shape[0])) ** 2 + (y / (0.38 * shape[1])) ** 2 + (z / (0.34 * shape[2])) ** 2
    box = (abs(x + 3) < 2.5) & (abs(y - 2) < 3.5) & (abs(z) < 3)
    return (0.7 * (r2 < 1) + 0.3 * np.exp(-3 * r2) - 0.3 * box).astype(dtype)


def band_limited(shape, shift=(0.0, 0.0, 0.0)):
    """A positive object of a few low spatial frequencies, moved by ``shift`` voxels."""
    grid = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in shape], indexing="ij")
    waves = [((1, 0, 0), 0.5, 0.3), ((0, 2, 1), 0.3, 1.1), ((2, -1, 3), 0.2, 2.0), ((-3, 4, -2), 0.25, 0.7)]
    volume = np.full(shape, 2.0)
    for k, amplitude, phase in waves:
        angle = sum(2 * np.pi * k[a] * (grid[a] - shift[a]) / shape[a] for a in range(3))
        volume += amplitude * np.cos(angle + phase)
    return volume


def centre_of_mass(volume):
    volume = np.asarray(volume, dtype=np.float64)
    grid = np.meshgrid(*[np.arange(n) for n in volume.shape], indexing="ij")
    return np.array([(volume * g).sum() / volume.sum() for g in grid])


def gain(after, before):
    """What an artefact multiplied every point of k-space by, where there was signal."""
    after, before = np.fft.fftn(after), np.fft.fftn(before)
    return np.where(np.abs(before) > 1e-9, after / np.where(before == 0, 1, before), np.nan)


def distance(shape, axis):
    """Distance of every point of k-space from its centre along an axis, 1 at Nyquist."""
    index = np.arange(shape[axis])
    k = np.minimum(index, shape[axis] - index) / (shape[axis] / 2)
    return np.broadcast_to(k.reshape([-1 if a == axis else 1 for a in range(3)]), shape)


class TestTheMotionOptions(unittest.TestCase):
    shape = (40, 48, 36)

    def test_the_severity_table_of_the_documentation(self):
        for severity, events, amplitude in [(1, 1, 1.0), (2, 3, 2.8), (3, 5, 5.2)]:
            motion = kspace.resolve_motion(severity, IDENTITY, self.shape)
            self.assertEqual(motion["events"], events)
            self.assertAlmostEqual(motion["translation"], amplitude, delta=0.05)
            self.assertAlmostEqual(motion["rotation"], amplitude, delta=0.05)

    def test_events_and_amplitudes_grow_with_the_severity_to_the_power_of_one_and_a_half(self):
        for severity in (0.4, 1.5, 2.5, 4):
            motion = kspace.resolve_motion(severity, IDENTITY, self.shape)
            self.assertEqual(motion["translation"], severity ** 1.5)
            self.assertEqual(motion["rotation"], severity ** 1.5)
            self.assertEqual(motion["events"], max(1, int(np.floor(severity ** 1.5 + 0.5))))

    def test_the_defaults(self):
        motion = kspace.resolve_motion(1, IDENTITY, self.shape)
        self.assertEqual(
            {key: motion[key] for key in ("blocks", "continuous", "pe", "ordering", "centre", "lines")},
            {"blocks": 32, "continuous": 0.2, "pe": 1, "ordering": "linear", "centre": True, "lines": 48})

    def test_no_motion_is_none(self):
        for options in (None, 0, -1, {"severity": 0}):
            self.assertIsNone(kspace.resolve_motion(options, IDENTITY, self.shape))

    def test_a_dict_without_a_severity_is_mild_motion(self):
        self.assertEqual(kspace.resolve_motion({}, IDENTITY, self.shape),
                         kspace.resolve_motion(1, IDENTITY, self.shape))

    def test_a_dict_overrides_single_values_and_the_severity_gives_the_rest(self):
        motion = kspace.resolve_motion(
            {"severity": 2, "events": 6, "translation": 4, "ordering": "Centric", "centre": 0},
            IDENTITY, self.shape)
        self.assertEqual((motion["events"], motion["translation"], motion["rotation"]), (6, 4.0, 2 ** 1.5))
        self.assertEqual((motion["ordering"], motion["centre"]), ("centric", False))

    def test_there_are_no_more_blocks_than_lines(self):
        motion = kspace.resolve_motion({"severity": 1, "pe": 2}, IDENTITY, (40, 48, 20))
        self.assertEqual((motion["blocks"], motion["lines"]), (20, 20))

    def test_a_world_axis_is_the_closest_voxel_axis_of_an_oblique_permuted_grid(self):
        affine = oblique_affine()
        for world, voxel in (("x", 2), ("y", 0), ("z", 1), ("Y", 0)):
            motion = kspace.resolve_motion({"pe": world}, affine, self.shape)
            self.assertEqual((motion["pe"], motion["lines"]), (voxel, self.shape[voxel]))
            self.assertEqual(kspace.resolve_ringing({"pe": world}, affine)["axes"], [voxel])

    def test_a_voxel_axis_is_taken_as_it_is(self):
        for axis in (0, 1, 2):
            self.assertEqual(kspace.resolve_motion({"pe": axis}, oblique_affine(), self.shape)["pe"], axis)

    def test_what_makes_no_sense_is_an_error(self):
        wrong = [
            ({"ordering": "spiral"}, "ordering"), ({"blocks": 1}, "blocks"), ({"events": 0}, "events"),
            ({"translation": -1}, "negative"), ({"translation": 0, "rotation": 0}, "nothing would move"),
            ({"continuous": -0.1}, "continuous"), ({"pe": 3}, "pe"), ({"pe": "w"}, "pe"),
            ({"center": 1}, "center"), ({"severity": "mild"}, "severity"), ("mild", "scalar"),
            ([1, 2], "scalar"),
        ]
        for options, message in wrong:
            with self.assertRaisesRegex(ValueError, message):
                kspace.resolve_motion(options, IDENTITY, self.shape)


class TestTheRingingOptions(unittest.TestCase):
    def test_a_scalar_is_the_strength_of_a_notch_along_y(self):
        self.assertEqual(kspace.resolve_ringing(2, IDENTITY),
                         {"strength": 2.0, "type": "notch", "axes": [1], "k0": 0.6, "gain": -1.0})

    def test_gibbs_keeps_less_of_k_space_the_stronger_it_is(self):
        for strength, fraction in [(1, 0.55), (2, 0.40), (3, 0.25), (5, 0.15)]:
            ringing = kspace.resolve_ringing({"strength": strength, "type": "Gibbs"}, IDENTITY)
            self.assertAlmostEqual(ringing["fraction"], fraction)
            self.assertNotIn("gain", ringing)

    def test_all_axes(self):
        self.assertEqual(kspace.resolve_ringing({"pe": "all"}, IDENTITY)["axes"], [0, 1, 2])

    def test_no_ringing_is_none(self):
        for options in (None, 0, {"strength": 0}):
            self.assertIsNone(kspace.resolve_ringing(options, IDENTITY))

    def test_what_makes_no_sense_is_an_error(self):
        wrong = [(-1, "negative"), ({"type": "hann"}, "type"), ({"k0": 1}, "k0"), ({"k0": 0}, "k0"),
                 ({"pe": 4}, "pe"), ({"fraction": 0.5}, "fraction")]
        for options, message in wrong:
            with self.assertRaisesRegex(ValueError, message):
                kspace.resolve_ringing(options, IDENTITY)


class TestTheAcquisition(unittest.TestCase):
    def test_linear_runs_from_minus_kmax_to_kmax(self):
        line, block, centre = kspace._acquisition(8, 4, "linear")
        np.testing.assert_array_equal(line, [4, 5, 6, 7, 0, 1, 2, 3])      # k = -4 ... 3 in the FFT
        np.testing.assert_array_equal(block, [0, 0, 1, 1, 2, 2, 3, 3])
        self.assertEqual(centre, 2)

    def test_centric_starts_at_the_centre_and_alternates(self):
        line, _, centre = kspace._acquisition(7, 3, "centric")
        np.testing.assert_array_equal(line, [0, 1, 6, 2, 5, 3, 4])         # k = 0, 1, -1, 2, -2, 3, -3
        self.assertEqual(centre, 0)

    def test_every_line_is_acquired_once(self):
        for lines, blocks, ordering in [(48, 32, "linear"), (47, 32, "centric"), (31, 31, "linear")]:
            line, block, _ = kspace._acquisition(lines, blocks, ordering)
            np.testing.assert_array_equal(np.sort(line), np.arange(lines))
            np.testing.assert_array_equal(np.unique(block), np.arange(blocks))

    def test_the_blocks_end_where_matlab_rounds_them(self):
        # round(linspace(0, 10, 5)) is 0 3 5 8 10 in MATLAB, which rounds 2.5 up
        _, block, _ = kspace._acquisition(10, 4, "linear")
        np.testing.assert_array_equal(np.bincount(block), [3, 2, 3, 2])


class TestThePoseTimeCourse(unittest.TestCase):
    shape = (40, 48, 36)

    def motion(self, **options):
        return kspace.resolve_motion({"severity": 2, **options}, IDENTITY, self.shape)

    def test_one_row_per_block_and_one_block_per_event(self):
        pose, events = kspace.sample_motion(self.motion(), 3)
        self.assertEqual(pose.shape, (32, 6))
        self.assertEqual(len(events), 3)
        self.assertTrue(all(0 <= at < 32 for at in events))

    def test_a_seed_gives_the_same_motion_and_another_seed_another(self):
        first, events = kspace.sample_motion(self.motion(), 7)
        again, events_again = kspace.sample_motion(self.motion(), 7)
        other, _ = kspace.sample_motion(self.motion(), 8)
        np.testing.assert_array_equal(first, again)
        self.assertEqual(events, events_again)
        self.assertGreater(np.abs(first - other).max(), 0.1)

    def test_the_first_event_hits_the_centre_of_k_space(self):
        for ordering, centre in (("linear", 16), ("centric", 0)):
            for seed in range(5):
                _, events = kspace.sample_motion(self.motion(ordering=ordering), seed)
                self.assertEqual(events[0], centre)

    def test_without_centre_the_events_are_anywhere(self):
        first = {kspace.sample_motion(self.motion(centre=False), seed)[1][0] for seed in range(20)}
        self.assertGreater(len(first), 5)

    def test_an_event_is_an_excursion_and_a_residual_offset_of_30_percent(self):
        motion = self.motion(events=1, continuous=0, translation=3, rotation=2)
        for seed in range(10):
            pose, (at,) = kspace.sample_motion(motion, seed)
            np.testing.assert_array_equal(pose[:at], 0)
            np.testing.assert_allclose(pose[at + 1:], np.tile(0.3 * pose[at], (31 - at, 1)))
            # one amplitude for both, between 30% and 100% of the maximum
            amplitude = np.linalg.norm(pose[at, :3]) / 3
            self.assertAlmostEqual(np.linalg.norm(pose[at, 3:]) / 2, amplitude)
            self.assertTrue(0.3 <= amplitude <= 1.0)

    def test_pitch_dominates_the_nods(self):
        motion = self.motion(events=1, continuous=0)
        nods = np.array([pose[at] for pose, (at,) in
                         (kspace.sample_motion(motion, seed) for seed in range(300))])
        size = np.abs(nods).mean(axis=0)
        self.assertLess(size[0], 0.6 * min(size[1], size[2]))       # little left-right translation
        self.assertGreater(size[3], 1.6 * max(size[4], size[5]))    # mostly rotation around x

    def test_the_continuous_part_is_a_bounded_translation(self):
        events_only, _ = kspace.sample_motion(self.motion(continuous=0), 5)
        with_drift, _ = kspace.sample_motion(self.motion(continuous=0.2), 5)
        drift = with_drift - events_only
        np.testing.assert_array_equal(drift[:, 3:], 0)
        self.assertTrue((np.ptp(drift[:, :3], axis=0) > 0.1).all())
        # drift of 0.5 and oscillations of 0.3 and 0.2 of the amplitude
        self.assertLessEqual(np.abs(drift).max(), 0.2 * 2 ** 1.5)
        self.assertGreater(len(np.unique(with_drift[:, :3], axis=0)), 30)

    def test_a_generator_is_used_as_it_is(self):
        pose, _ = kspace.sample_motion(self.motion(), np.random.default_rng(1))
        again, _ = kspace.sample_motion(self.motion(), np.random.default_rng(1))
        seeded, _ = kspace.sample_motion(self.motion(), 1)
        np.testing.assert_array_equal(pose, again)
        self.assertGreater(np.abs(pose - seeded).max(), 0.1)        # the seed has a stream of its own


class TestTheMotion(unittest.TestCase):
    def setUp(self):
        self.volume = phantom()
        self.motion = kspace.resolve_motion({"severity": 2, "blocks": 14}, IDENTITY, self.volume.shape)

    def step(self):
        """The head moves by (-2, 3, 0) mm after the fourth block and stays there."""
        pose = np.zeros((14, 6))
        pose[4:, :2] = [-2.0, 3.0]
        return pose

    def test_a_head_that_does_not_move_gives_the_image(self):
        out, realised = kspace.apply_motion(self.volume, np.zeros((14, 6)), self.motion, IDENTITY)
        np.testing.assert_array_equal(out, self.volume)
        self.assertIsNot(out, self.volume)
        np.testing.assert_array_equal(realised["Pose"], np.zeros((14, 6)))

    def test_one_pose_for_all_blocks_is_no_motion(self):
        pose = np.tile([1.5, -2.0, 0.7, 2.0, -1.0, 3.0], (14, 1))
        out, realised = kspace.apply_motion(self.volume, pose, self.motion, IDENTITY)
        np.testing.assert_allclose(out, self.volume, atol=1e-5)
        np.testing.assert_allclose(realised["Pose"], 0, atol=1e-12)

    def test_a_translation_of_all_blocks_is_the_shifted_image(self):
        shape, affine = (20, 24, 18), oblique_affine()
        motion = kspace.resolve_motion({"severity": 1, "blocks": 6}, affine, shape)
        millimetres = np.array([1.3, -2.45, 0.6])
        voxels = np.linalg.solve(affine[:3, :3], millimetres)
        pose = np.tile(np.concatenate([millimetres, np.zeros(3)]), (6, 1))
        for dtype, tolerance in ((np.float64, 1e-12), (np.float32, 1e-5)):
            out, _ = kspace.apply_motion(band_limited(shape).astype(dtype), pose, motion, affine,
                                         subtract_mean=False, recentre=False)
            np.testing.assert_allclose(out, band_limited(shape, voxels), atol=tolerance)

    def test_translations_only_move_the_energy_around(self):
        volume = phantom(dtype=np.float64)
        pose = np.zeros((14, 6))
        pose[:, :3] = np.random.default_rng(0).normal(size=(14, 3))
        out, _ = kspace.apply_motion(volume, pose, self.motion, oblique_affine())
        self.assertGreater(np.abs(out - volume).max(), 0.1)
        self.assertAlmostEqual((out ** 2).sum(), (volume ** 2).sum(), places=8)

    def test_the_mean_pose_is_subtracted_and_the_object_stays_in_place(self):
        at = centre_of_mass(self.volume)
        kept, _ = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY, recentre=False)
        moved, _ = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY, subtract_mean=False,
                                       recentre=False)
        self.assertLess(np.linalg.norm(centre_of_mass(kept) - at), 0.25)
        self.assertGreater(np.linalg.norm(centre_of_mass(moved) - at), 2.5)

    def test_the_subtracted_mean_is_weighted_by_the_energy_of_the_blocks(self):
        _, realised = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY, recentre=False)
        applied = np.array(realised["Pose"])
        # the blocks before the step hold the outer lines, with 1% of the energy
        energy = np.abs(np.fft.fft(self.volume, axis=1)) ** 2
        line, block, _ = kspace._acquisition(28, 14, "linear")
        early = energy[:, line[block < 4]].sum() / energy.sum()
        self.assertLess(early, 0.02)
        np.testing.assert_allclose(applied[4:, :2], np.tile([-2.0 * early, 3.0 * early], (10, 1)), atol=1e-6)
        np.testing.assert_allclose(applied[4:] - applied[:4].mean(axis=0), self.step()[4:], atol=1e-12)

    def test_a_pose_all_blocks_share_changes_nothing(self):
        out, realised = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY)
        shared = self.step() + [5.0, -7.0, 2.0, 1.0, -2.0, 3.0]
        same, same_realised = kspace.apply_motion(self.volume, shared, self.motion, IDENTITY)
        np.testing.assert_allclose(same, out, atol=1e-5)
        np.testing.assert_allclose(same_realised["Pose"], realised["Pose"], atol=1e-12)

    def test_blocks_with_one_rotation_share_one_resampling(self):
        pose = np.zeros((14, 6))
        pose[:, 0] = np.linspace(-1, 1, 14)             # every block has its own translation
        pose[3:9, 3:] = [2.0, -1.0, 0.5]
        pose[9:, 3:] = [0.6, -0.3, 0.15]
        with mock.patch.object(kspace, "_resample", wraps=kspace._resample) as resample:
            kspace.apply_motion(self.volume, pose, self.motion, IDENTITY, subtract_mean=False)
        self.assertEqual(resample.call_count, 2)        # and none for the blocks that are not rotated

    def test_the_realised_motion_is_what_the_sidecar_gets(self):
        pose, events = kspace.sample_motion(self.motion, 4)
        _, realised = kspace.apply_motion(self.volume, pose, self.motion, IDENTITY, events=events)
        self.assertEqual(json.loads(json.dumps(realised)), realised)
        self.assertEqual(
            {key: value for key, value in realised.items() if key not in ("Pose", "Recentred")},
            {"Severity": 2.0, "Events": 3, "EventBlocks": events, "MaxTranslation": 2 ** 1.5,
             "MaxRotation": 2 ** 1.5, "Blocks": 14, "PhaseEncodingAxis": 1, "Ordering": "linear",
             "CentreBlock": 7, "Continuous": 0.2})
        self.assertEqual(np.shape(realised["Recentred"]), (3,))
        self.assertEqual(np.shape(realised["Pose"]), (14, 6))
        self.assertEqual(events[0], 7)

    def test_a_seed_gives_the_same_image_and_another_seed_another(self):
        first, realised = kspace.add_motion(self.volume, self.motion, IDENTITY, 11)
        again, realised_again = kspace.add_motion(self.volume, self.motion, IDENTITY, 11)
        other, _ = kspace.add_motion(self.volume, self.motion, IDENTITY, 12)
        np.testing.assert_array_equal(first, again)
        self.assertEqual(realised, realised_again)
        self.assertGreater(np.abs(first - other).max(), 0.05)
        self.assertGreater(np.abs(first - self.volume).max(), 0.05)
        self.assertGreaterEqual(first.min(), 0)         # a magnitude image

    def test_add_motion_is_sample_and_apply(self):
        pose, events = kspace.sample_motion(self.motion, 11)
        expected, realised = kspace.apply_motion(self.volume, pose, self.motion, IDENTITY, events=events)
        out, realised_out = kspace.add_motion(self.volume, self.motion, IDENTITY, 11)
        np.testing.assert_array_equal(out, expected)
        self.assertEqual(realised_out, realised)

    def test_without_motion_add_motion_adds_the_ringing_alone(self):
        ringing = kspace.resolve_ringing(2, IDENTITY)
        out, realised = kspace.add_motion(self.volume, None, IDENTITY, 11, ringing)
        np.testing.assert_array_equal(out, kspace.add_ringing(self.volume, ringing))
        self.assertIsNone(realised)

    def test_the_applied_pose_reproduces_the_image(self):
        out, realised = kspace.add_motion(self.volume, self.motion, IDENTITY, 11)
        replayed, _ = kspace.apply_motion(self.volume, realised["Pose"], self.motion, IDENTITY,
                                          subtract_mean=False, recentre=False)
        np.testing.assert_allclose(replayed, out, atol=2e-5)

    def test_the_image_is_brought_back_onto_the_object(self):
        # A step in the pose that the mean does not take out: the magnitude
        # image lies elsewhere than the complex one.
        reference = np.fft.rfftn(self.volume).conj()

        def offset(image):
            import torch
            return kspace._displacement(torch.as_tensor(np.ascontiguousarray(image)),
                                        torch.as_tensor(reference), coarse=True)

        left, plain = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY, recentre=False)
        centred, realised = kspace.apply_motion(self.volume, self.step(), self.motion, IDENTITY)
        self.assertGreater(np.abs(offset(left)).max(), 0.02)
        self.assertLess(np.abs(offset(centred)).max(), 2e-3)
        self.assertIsNone(plain["Recentred"])
        # the pose is the movement relative to the image that is returned
        np.testing.assert_allclose(np.array(realised["Pose"])[:, :3],
                                   np.array(plain["Pose"])[:, :3] - realised["Recentred"], atol=1e-12)
        np.testing.assert_array_equal(np.array(realised["Pose"])[:, 3:], np.array(plain["Pose"])[:, 3:])

    def test_a_whole_translation_is_taken_out_again(self):
        pose = np.tile([3.0, -2.0, 1.0, 0, 0, 0], (14, 1)).astype(float)
        out, realised = kspace.apply_motion(self.volume, pose, self.motion, IDENTITY, subtract_mean=False)
        np.testing.assert_allclose(out, self.volume, atol=1e-4)
        np.testing.assert_allclose(realised["Recentred"], [3.0, -2.0, 1.0], atol=1e-3)
        np.testing.assert_allclose(realised["Pose"], 0, atol=1e-3)

    def test_recentring_leaves_a_head_that_did_not_move_alone(self):
        out, realised = kspace.apply_motion(self.volume, np.zeros((14, 6)), self.motion, IDENTITY)
        np.testing.assert_array_equal(out, self.volume)
        self.assertEqual(realised["Recentred"], [0.0, 0.0, 0.0])

    def test_the_removed_shift_is_given_in_mm_of_the_world(self):
        affine = oblique_affine()
        volume = phantom(dtype=np.float64)
        motion = kspace.resolve_motion({"severity": 2, "blocks": 14}, affine, volume.shape)
        millimetres = np.array([1.2, -0.9, 0.5])
        pose = np.tile(np.concatenate([millimetres, np.zeros(3)]), (14, 1))
        _, realised = kspace.apply_motion(volume, pose, motion, affine, subtract_mean=False)
        np.testing.assert_allclose(realised["Recentred"], millimetres, atol=5e-3)

    def test_a_volume_or_a_pose_of_another_size_is_refused(self):
        with self.assertRaisesRegex(ValueError, "lines"):
            kspace.apply_motion(phantom((24, 30, 20)), np.zeros((14, 6)), self.motion, IDENTITY)
        with self.assertRaisesRegex(ValueError, "pose"):
            kspace.apply_motion(self.volume, np.zeros((32, 6)), self.motion, IDENTITY)
        with self.assertRaisesRegex(ValueError, "3-D"):
            kspace.add_ringing(np.zeros((4, 4)), None)

    def test_options_that_were_not_resolved_are_refused(self):
        with self.assertRaisesRegex(TypeError, "resolve_motion"):
            kspace.add_motion(self.volume, 2, IDENTITY, 11)
        with self.assertRaisesRegex(TypeError, "resolve_motion"):
            kspace.sample_motion({"severity": 2}, 11)
        with self.assertRaisesRegex(TypeError, "resolve_ringing"):
            kspace.add_ringing(self.volume, 2)
        with self.assertRaisesRegex(TypeError, "resolve_ringing"):
            kspace.add_motion(self.volume, self.motion, IDENTITY, 11, {"strength": 2})


class TestTheRotation(unittest.TestCase):
    def rotate(self, volume, degrees, affine):
        return kspace._resample(torch.as_tensor(volume),
                                kspace._rotation(degrees, affine, volume.shape)).numpy()

    def test_the_angles_are_those_of_spm_matrix(self):
        # a quarter turn lands on the grid, where the kernel is a single point
        volume = np.random.default_rng(0).random((9, 9, 9))
        for axis, plane in enumerate([(1, 2), (0, 2), (0, 1)]):
            degrees = [0.0, 0.0, 0.0]
            degrees[axis] = 90.0
            np.testing.assert_array_equal(self.rotate(volume, degrees, IDENTITY),
                                          np.rot90(volume, -1, axes=plane))

    def test_it_is_a_rotation_of_the_world_around_the_centre_of_the_volume(self):
        shape, affine = (32, 36, 30), oblique_affine()
        index = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1)
        world = index @ affine[:3, :3].T + affine[:3, 3]
        centre = affine[:3, :3] @ ((np.array(shape) - 1) / 2) + affine[:3, 3]
        blob = centre + [3.0, -2.5, 2.0]

        def gaussian(at):
            return np.exp(-((world - at) ** 2).sum(axis=-1) / (2 * 2.5 ** 2))

        # 20 degrees around the world's x: (y, z) -> (y cos + z sin, -y sin + z cos)
        c, s = np.cos(np.deg2rad(20)), np.sin(np.deg2rad(20))
        turned = centre + [3.0, -2.5 * c + 2.0 * s, 2.5 * s + 2.0 * c]
        out = self.rotate(gaussian(blob), [20.0, 0.0, 0.0], affine)
        # the blob is where it belongs to a micrometre; its values carry the
        # ripple of the six-point kernel, which sharpens between the voxels
        at = (out[..., None] * world).sum(axis=(0, 1, 2)) / out.sum()
        np.testing.assert_allclose(at, turned, atol=1e-3)
        self.assertGreater(np.linalg.norm(turned - blob), 1.0)
        np.testing.assert_allclose(out, gaussian(turned), atol=0.015)
        self.assertGreater(np.abs(gaussian(blob) - gaussian(turned)).max(), 0.25)

    def test_no_rotation_is_the_volume(self):
        volume = phantom(dtype=np.float64)
        np.testing.assert_array_equal(self.rotate(volume, [0.0, 0.0, 0.0], oblique_affine()), volume)

    def test_what_comes_from_outside_the_field_of_view_is_zero(self):
        out = self.rotate(np.ones((12, 12, 12)), [0.0, 0.0, 45.0], IDENTITY)
        self.assertEqual(out[0, 0, 6], 0)
        np.testing.assert_allclose(out[3:9, 3:9], 1, atol=1e-12)    # the weights sum to 1
        self.assertTrue(np.isfinite(out).all())

    def test_the_interpolation_is_the_sinc_of_spm(self):
        # rotate_volume of mri_simulate.m (spm_slice_vol with a hold of -5) for this volume
        i, j, k = np.meshgrid(np.arange(7), np.arange(6), np.arange(5), indexing="ij")
        volume = (np.sin(0.9 * i + 0.3 * j) + np.cos(0.5 * j - 0.7 * k) + 0.01 * i * k).astype(np.float32)
        affine = np.array([[0, 0, 1.5, -3], [-0.8, 0, 0, 2], [0, 1.2, 0, -1], [0, 0, 0, 1]])
        out = self.rotate(volume.astype(np.float64), [4.0, -3.0, 5.0], affine)
        np.testing.assert_allclose(
            out[:, 2, 2], [1.57283306, 1.9756273, 1.68007636, 0.872780442, 0.146738887,
                           0.0186910722, 0.635314345], atol=2e-7)
        np.testing.assert_allclose(
            out[3, :, 1], [1.33225834, 1.12004471, 0.686592937, 0.0976702049, -0.557103038, 0],
            atol=2e-7)
        np.testing.assert_array_equal(out[0, 0, :], 0)
        np.testing.assert_allclose(self.rotate(volume, [4.0, -3.0, 5.0], affine), out, atol=2e-6)


class TestTheRinging(unittest.TestCase):
    def test_gibbs_overshoots_a_step_by_about_9_percent_whatever_the_fraction(self):
        step = np.zeros((4, 96, 6))
        step[:, 24:72] = 1
        overshoot = []
        for strength in (1, 2, 3):
            out = kspace.add_ringing(step, kspace.resolve_ringing({"strength": strength, "type": "gibbs"},
                                                                  IDENTITY))
            # between the voxels as well: the peak need not fall onto one
            spectrum = np.fft.fftshift(np.fft.fft(out[0, :, 0]))
            fine = np.fft.ifft(np.fft.ifftshift(np.pad(spectrum, 96 * 8))).real * 17
            overshoot.append(fine.max() - 1)
            self.assertAlmostEqual(fine.min(), -overshoot[-1], places=6)
            np.testing.assert_allclose(out, out[:1, :, :1] + 0 * out, atol=1e-12)   # along y alone
        # 8.9% for a continuous step, a little more for a sampled one
        self.assertTrue(all(0.085 < value < 0.11 for value in overshoot), overshoot)

    def test_gibbs_truncates_the_phase_encoding_axis_alone(self):
        volume = np.random.default_rng(1).random((16, 20, 12))
        for axis in (0, 1, 2):
            ringing = kspace.resolve_ringing({"strength": 2, "type": "gibbs", "pe": axis}, IDENTITY)
            kept = distance(volume.shape, axis) <= ringing["fraction"]
            np.testing.assert_allclose(gain(kspace.add_ringing(volume, ringing), volume), kept, atol=1e-9)
            self.assertAlmostEqual(kept.mean(), 0.4, delta=0.07)

    def test_gibbs_along_all_axes_truncates_each_of_them(self):
        volume = np.random.default_rng(1).random((16, 20, 12))
        ringing = kspace.resolve_ringing({"strength": 1, "type": "gibbs", "pe": "all"}, IDENTITY)
        kept = np.all([distance(volume.shape, axis) <= ringing["fraction"] for axis in range(3)], axis=0)
        np.testing.assert_allclose(gain(kspace.add_ringing(volume, ringing), volume), kept, atol=1e-9)

    def test_the_gain_of_the_notch_is_one_minus_the_strength_in_a_symmetric_band_at_k0(self):
        volume = np.random.default_rng(2).random((12, 40, 10))
        for strength, k0 in [(1, 0.6), (2.5, 0.6), (3, 0.35)]:
            ringing = kspace.resolve_ringing({"strength": strength, "k0": k0}, IDENTITY)
            out = kspace.add_ringing(volume, ringing)
            k = distance(volume.shape, 1)
            band = (k > 0.95 * k0) & (k < 1.05 * k0)
            # a band on one side only would end up as the mean of the two gains
            np.testing.assert_allclose(gain(out, volume), np.where(band, 1 - strength, 1), atol=1e-9)
            rows = np.flatnonzero(band[0, :, 0])
            np.testing.assert_array_equal(np.sort(40 - rows), rows)
            self.assertIn(round(k0 * 20), rows)
            self.assertTrue(np.isrealobj(out))

    def test_the_notch_along_all_axes_is_a_spherical_shell(self):
        volume = np.random.default_rng(3).random((20, 24, 16))
        out = kspace.add_ringing(volume, kspace.resolve_ringing({"strength": 2, "pe": "all"}, IDENTITY))
        radius = np.sqrt(sum(distance(volume.shape, axis) ** 2 for axis in range(3)))
        shell = (radius > 0.57) & (radius < 0.63)
        self.assertGreater(shell.sum(), 100)
        np.testing.assert_allclose(gain(out, volume), np.where(shell, -1, 1), atol=1e-9)

    def test_ringing_keeps_the_mean_and_makes_ripples(self):
        volume = phantom(dtype=np.float64)
        for options in (2, {"strength": 2, "type": "gibbs"}):
            out = kspace.add_ringing(volume, kspace.resolve_ringing(options, IDENTITY))
            self.assertAlmostEqual(out.mean(), volume.mean(), places=12)
            self.assertGreater(np.abs(out - volume).max(), 0.05)
            self.assertLess(out.min(), -0.01)           # the real part, not the magnitude

    def test_no_ringing_is_a_copy_of_the_volume(self):
        volume = phantom()
        out = kspace.add_ringing(volume, None)
        np.testing.assert_array_equal(out, volume)
        self.assertIsNot(out, volume)

    def test_ringing_and_motion_share_one_k_space(self):
        volume = phantom(dtype=np.float64)
        motion = kspace.resolve_motion({"severity": 2, "blocks": 14}, IDENTITY, volume.shape)
        pose, _ = kspace.sample_motion(motion, 6)
        moved, _ = kspace.apply_motion(volume, pose, motion, IDENTITY)
        for options in (2, {"strength": 2, "type": "gibbs"}):
            ringing = kspace.resolve_ringing(options, IDENTITY)
            # the gain on the complex k-space of the motion, before the magnitude ...
            together, _ = kspace.apply_motion(volume, pose, motion, IDENTITY, ringing)
            self.assertGreater(np.abs(together - moved).max(), 0.05)
            self.assertGreaterEqual(together.min(), 0)
            # ... which is not the ringing of the magnitude image
            after = kspace.add_ringing(moved, ringing)
            self.assertGreater(np.abs(together - after).max(), 0.02)
            self.assertGreater(np.linalg.norm(together - after), 0.05 * np.linalg.norm(together - moved))


class TestWhatGoesInComesOut(unittest.TestCase):
    def setUp(self):
        self.volume = phantom()
        self.motion = kspace.resolve_motion({"severity": 1, "blocks": 14}, IDENTITY, self.volume.shape)
        self.ringing = kspace.resolve_ringing(1, IDENTITY)

    def test_an_array_gives_an_array_of_its_precision(self):
        for dtype in (np.float32, np.float64):
            out, _ = kspace.add_motion(self.volume.astype(dtype), self.motion, IDENTITY, 1, self.ringing)
            self.assertIsInstance(out, np.ndarray)
            self.assertEqual((out.dtype, out.shape), (dtype, self.volume.shape))
            self.assertEqual(kspace.add_ringing(self.volume.astype(dtype), self.ringing).dtype, dtype)

    def test_single_and_double_precision_agree(self):
        single, _ = kspace.add_motion(self.volume, self.motion, IDENTITY, 1, self.ringing)
        double, _ = kspace.add_motion(self.volume.astype(np.float64), self.motion, IDENTITY, 1, self.ringing)
        np.testing.assert_allclose(single, double, atol=1e-5)

    def test_integers_are_computed_and_returned_in_single_precision(self):
        counts = np.round(100 * self.volume).astype(np.int16)
        out = kspace.add_ringing(counts, self.ringing)
        self.assertEqual(out.dtype, np.float32)
        np.testing.assert_allclose(out, kspace.add_ringing(counts.astype(np.float32), self.ringing))

    def test_a_read_only_array_is_not_touched(self):
        volume = self.volume.copy()
        volume.flags.writeable = False
        out, _ = kspace.add_motion(volume, self.motion, IDENTITY, 1)
        np.testing.assert_array_equal(volume, self.volume)
        np.testing.assert_array_equal(out, kspace.add_motion(self.volume, self.motion, IDENTITY, 1)[0])

    def test_a_tensor_gives_a_tensor_of_its_precision(self):
        expected, _ = kspace.add_motion(self.volume, self.motion, IDENTITY, 1, self.ringing)
        for dtype in (torch.float32, torch.float64):
            tensor = torch.from_numpy(self.volume).to(dtype)
            out, _ = kspace.add_motion(tensor, self.motion, IDENTITY, 1, self.ringing)
            self.assertIsInstance(out, torch.Tensor)
            self.assertEqual((out.dtype, tuple(out.shape)), (dtype, self.volume.shape))
            np.testing.assert_allclose(out.numpy(), expected, atol=1e-5)
            self.assertIsInstance(kspace.add_ringing(tensor, self.ringing), torch.Tensor)
        np.testing.assert_array_equal(tensor.numpy(), self.volume.astype(np.float64))

    @unittest.skipUnless(torch.backends.mps.is_available(), "no MPS device")
    def test_a_tensor_on_mps_stays_there_and_gives_the_image_of_the_cpu(self):
        tensor = torch.from_numpy(self.volume).to("mps")
        motion = dict(self.motion)
        expected, _ = kspace.add_motion(self.volume, motion, IDENTITY, 1, self.ringing)
        out, _ = kspace.add_motion(tensor, motion, IDENTITY, 1, self.ringing)
        self.assertEqual((out.device.type, out.dtype), ("mps", torch.float32))
        np.testing.assert_allclose(out.cpu().numpy(), expected, atol=1e-5)
        ringing = kspace.add_ringing(tensor, self.ringing)
        self.assertEqual(ringing.device.type, "mps")
        np.testing.assert_allclose(ringing.cpu().numpy(), kspace.add_ringing(self.volume, self.ringing),
                                   atol=1e-5)


if __name__ == "__main__":
    unittest.main()
