# Simulating images

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim.simulate` runs the whole path from a real scan to simulated images with their ground truth. It has three steps, which can also be called on their own.

- [In one call](#in-one-call)
- [The three steps](#the-three-steps)
- [Options of the anatomy](#options-of-the-anatomy)
- [Options of the acquisition](#options-of-the-acquisition)
- [What is written](#what-is-written)
- [Command line](#command-line)
- [How long it takes](#how-long-it-takes)

## In one call

```python
from t1sim import simulate

result = simulate.simulate("sub-01_T1w.nii.gz", snr=30, bias_percent=20)
result["files"]["T1w"]     # the simulated image
result["files"]["dseg"]    # its ground truth label map
result["files"]["probseg"] # {"CSF": ..., "GM": ..., "WM": ...}: the tissue fractions
```

## The three steps

```python
truth = simulate.ground_truth("sub-01_T1w.nii.gz")          # once per anatomy, about 2.5 minutes
for snr in (20, 40):
    acquired = simulate.acquire(truth, snr=snr)             # seconds each
    simulate.write(truth, acquired)
```

| Step | Does | Cost |
| --- | --- | --- |
| `ground_truth(scan, ...)` | Runs the [front end](frontend.md), changes the [anatomy](anatomy.md) if asked, [cleans](cleanup.md) the ground truth and renders the image without artefacts | Minutes, once |
| `acquire(truth, ...)` | Makes one image: [bias field, voxel size, contrast, movement, ringing, noise](artefacts.md) | Seconds |
| `write(truth, acquired, ...)` | Stores the image and the ground truth with names and sidecars after BIDS | Seconds |

Any number of acquisitions can share one ground truth. The image is rendered from the fractions before the cleanup and the ground truth is the fractions after it, so a vessel stays bright in the image and is CSF in the ground truth.

## Options of the anatomy

Arguments of `ground_truth` (and of `simulate`):

| Argument | Meaning (default) |
| --- | --- |
| `thickness` | Give the cortex this thickness in mm; one value, or three for the mainly occipital parcels, the rest and the mainly frontal parcels (`0`: keep the cortex) |
| `atrophy` | A dict, or a list of dicts, with `regions` (numbers of the atlas) and `change` (relative volume change of their grey matter, e.g. `-0.1`), and optionally `atlas` (`"Hammers"`), `mode` (`"thinning"` or `"volume"`), `tissues`, `reach` (`None`) |
| `clean` | Clean the ground truth: `True`, `False`, or a dict with `vessels`, `dura`, `pve` (`True`) |
| `seg` | A result of `frontend.segment(scan, deformation=True)` to start from |

`thickness` and `atrophy` cannot be combined. See [Changing the anatomy](anatomy.md).

## Options of the acquisition

Arguments of `acquire` (and of `simulate`):

| Argument | Meaning (default) |
| --- | --- |
| `snr` | Rician noise for this SNR in white matter; `0` uses `noise_percent` (`40`) |
| `noise_percent` | Gaussian noise in percent of the white matter intensity, with `snr=0` (`0`) |
| `bias_percent` | Peak-to-peak amplitude of the bias field in percent; `0` for none, negative to invert it (`30`) |
| `bias_type` | `(strength, seed)` of a simulated field, or `"A"`, `"B"`, `"C"` for a predefined one. A seed of `None` gives every scan a field of its own (`(2, 0)`) |
| `contrast` | Exponent of the contrast change (`1`) |
| `resolution` | Voxel size of the acquisition in mm, one value or three (`None`: that of the grid) |
| `upsample` | Bring the acquired image back onto the 0.5 mm grid; the ground truth then keeps its voxels (`False`) |
| `motion` | Movement artefacts: a severity, or a dict (`0`) |
| `ringing` | Ringing: a strength, or a dict with `strength`, `type` (`"notch"` or `"gibbs"`), `pe`, `k0` (`0`) |
| `seed` | Seed of noise and movement. `None`: derived from the name of the scan and the settings, so every output file has noise of its own that is the same in every run (`None`) |
| `grid` | `"working"` (339 x 411 x 339), or `"deepmriprep"` for the 336 x 384 x 336 voxels its networks are trained on (`"working"`) |

All outputs lie in the space of the working grid: affinely registered to MNI152NLin2009cAsym.

## What is written

For `sub-01_T1w.nii.gz` with `snr=30, bias_percent=20`:

```
sub-01_space-MNI152NLin2009cAsym_res-0p5mm_desc-snr30Rf20T2_T1w.nii.gz     the simulated image
sub-01_space-MNI152NLin2009cAsym_res-0p5mm_desc-snr30Rf20T2_T1w.json       its parameters
sub-01_space-MNI152NLin2009cAsym_res-0p5mm_desc-Clean_dseg.nii.gz          label map of the ground truth
sub-01_space-MNI152NLin2009cAsym_res-0p5mm_desc-Clean_dseg.json
sub-01_space-MNI152NLin2009cAsym_res-0p5mm_desc-Clean_label-GM_probseg.nii.gz   and _label-WM_, _label-CSF_
```

- **Where.** In `derivatives/t1sim-<version>/` of the dataset: at the root of a BIDS dataset with the folders of subject and session kept, next to the folder of the scan otherwise. `out_dir` names another folder, `derivatives=False` writes next to the scan.
- **Names.** The scheme is that of `matlab/mri_simulate.m`. The `desc` of the image holds the tags of the acquisition followed by those of the anatomy; the `desc` of the ground truth holds those of the anatomy and `Clean`. Images that differ in the acquisition alone share one ground truth.
- **Ground truth.** The fractions (`probseg`, 8 bit in steps of 1/255) are the ground truth. The label map (`dseg`, 16 bit in steps of 1/765) is the view of them that is continuous: 2.5 is half GM and half WM.
- **Sidecar.** `SimulationParameters` holds every setting, the seed, the noise SD that was used, the tissue intensities, the volumes the cleanup relabelled, the volume change an atrophy reached, and for movement the pose of every block of k-space.
- **Atrophy.** The Jacobian determinant of the deformation is written as `..._jacobian.nii.gz`.

Tags of the acquisition: `snr30` or `pn3`, `Rf20T2` / `Rf20A` / `RfNeg20T2`, `ConLow` / `ConHigh` / `Con1p3`, `Acq1mm` (acquired at 1 mm and brought back), `Motion2`, `Ringing2` / `Ringing2Gibbs`, `Seed5` (only for a seed that was given). Tags of the anatomy: `Thickness25mm`, `Thickness15to25mm`, `HammersRoi28Thin10`, `HammersRoi28Vol10`, `HammersMulti`.

## Command line

```bash
t1sim sub-01_T1w.nii.gz --snr 20 40 --bias 20                       # two images of one ground truth
t1sim sub-01_T1w.nii.gz --resolution 1 1 1.5 --motion 2             # a 1 x 1 x 1.5 mm scan with movement
t1sim sub-01_T1w.nii.gz --atrophy Hammers 28 -0.1                   # 10% less grey matter in a region
t1sim sub-01_T1w.nii.gz --thickness 2.5                             # a cortex of 2.5 mm
t1sim sub-*_T1w.nii.gz --grid deepmriprep --resolution 1 --upsample --bias-type 2 scan   # training data
```

`t1sim --help` lists all options. Several values of `--snr` or `--seed` give several images of one ground truth. Without the installed command: `python -m t1sim.simulate`.

## How long it takes

Measured on an M2 MacBook Air, one 0.75 mm scan, 0.5 mm working grid:

| Step | Time |
| --- | --- |
| Front end with the deformation | 91 s |
| Cleanup | 35 to 40 s |
| Atrophy of one region | 20 to 30 s |
| Cortical thickness | 7 to 8 minutes |
| One acquisition on the 0.5 mm grid | 8 s |
| One acquisition at 1 mm, brought back to 0.5 mm | 18 s |
| Movement at 1 mm, severity 1 / 2 / 3 | 10 / 15 / 23 s |
