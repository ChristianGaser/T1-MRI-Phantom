# Artefacts of the acquisition

Part of the [documentation of t1sim](../README.md#documentation).

`simulate.acquire` turns the image without artefacts into a simulated scan. Every step takes the image and explicit parameters; the ground truth is never touched. The functions are in `t1sim.artefacts` (bias field, voxel size, contrast, noise) and `t1sim.kspace` (movement, ringing).

- [Order of the steps](#order-of-the-steps)
- [Bias field](#bias-field)
- [Voxel size](#voxel-size)
- [Contrast](#contrast)
- [Movement and ringing](#movement-and-ringing)
- [Noise](#noise)
- [What was checked](#what-was-checked)

![One detail under every artefact](t1sim-examples.png)

## Order of the steps

1. **Bias field.** It belongs to the object and the coil, so it comes first and everything after it acts on the biased image.
2. **Voxel size.** The image is acquired on the grid of the acquisition.
3. **Contrast.**
4. **Movement and ringing**, on the grid of the acquisition, because they are made on the acquisition matrix.
5. **Noise**, last: the noise of the scanner is no part of the moving object, and it is white on the grid of the acquisition.
6. Optionally the image is brought back onto the 0.5 mm grid (`upsample=True`).

## Bias field

A smooth field with a mean of 1 that the image is multiplied with.

- **Simulated** (`bias_type=(strength, seed)`): a random field of `2**strength` voxels per axis, low-pass filtered, brought to the size of the image and smoothed with 15 mm. A strength of 3 to 4 gives the finer field of 7 T. `bias_percent` is its peak-to-peak amplitude; a negative value swaps bright and dark.
- **Predefined** (`bias_type="A"`, `"B"`, `"C"`): measured fields in template space (`matlab/rf100_*.nii`), brought onto the grid through the deformation of the segmentation.

For one seed the simulated field is the field `matlab/mri_simulate.m` makes.

## Voxel size

`resolution` is the voxel size of the acquisition. This is not a port of the MATLAB version, which interpolates the image at the new voxel centres and the label linearly.

- **The image is acquired.** Its k-space is cut at the Nyquist frequency of the new voxel size along every axis that gets coarser, which is what a smaller acquisition matrix does, Gibbs ringing included. The band-limited image is then read at the new voxel centres.
- **The ground truth is averaged.** Every new voxel holds the mean of the fractions over the box it covers, which is the partial volume of that voxel. Volumes are kept exactly.
- **The grid** has as many whole voxels as fit into the field of view, centred, so the voxel size is exactly the one asked for.

With `upsample=True` the acquired image is brought back onto the 0.5 mm grid by cubic B-spline interpolation, and the ground truth keeps its 0.5 mm voxels. This is a scan of that voxel size after resampling, the form in which most real data reaches a network that works at 0.5 mm.

## Contrast

`contrast` is the exponent of a power law: the image is scaled to 0..1, raised to it and scaled back. 0.5 lowers the contrast between the tissues, 1.5 raises it.

## Movement and ringing

Both are made in k-space; the model is that of `matlab/mri_simulate.m` and is described in full in [matlab/docs/artefacts.md](../matlab/docs/artefacts.md).

- **Movement** (`motion`): k-space is assembled from blocks of phase-encoding lines that are taken from rigidly moved copies of the image. Translations are exact phase ramps; rotations resample the volume with the windowed sinc SPM uses. Events are nods with a residual offset, with a slow drift and tremor between them. A severity of 1, 2, 3 is mild, moderate, severe.
- **Ringing** (`ringing`): `"notch"` damps or inverts a narrow band of k-space, which lays pronounced ripples over the image without blurring it; `"gibbs"` truncates k-space, which is what a finite acquisition matrix does. Together with movement the two share one k-space.

Voxel axes count from 0 here (`pe=1` is the second axis), and the sidecar holds the pose of every block.

One step is added to the MATLAB model: the image is brought back onto the object. `mri_simulate.m` subtracts the mean pose, weighted by the k-space energy of the blocks, so that the image stays where its ground truth is. That is where the complex image lies; the magnitude image lies elsewhere. `t1sim` estimates the displacement that is left from the cross-correlation with the unmoved image and removes it with a phase ramp on the assembled k-space, which is exact. The sidecar holds what was removed (`Recentred`), and its `Pose` is the movement relative to the image that is written.

## Noise

`snr` gives Rician noise: Gaussian noise of one standard deviation in the real and in the imaginary part of the signal, `sqrt((s + n1)**2 + n2**2)`. The standard deviation is the intensity of white matter in the image divided by `snr`.

- The intensity of white matter is taken after a contrast change, so the SNR is that of the image that is written.
- The bias field has a mean of 1 and is not part of it: where the field is dark, the local SNR is lower, as in a scan.
- The noise has no upper limit.
- `noise_percent` gives Gaussian noise instead, in percent of the white matter intensity.
- With `seed=None` every output file gets noise of its own.

With `upsample=True` the noise is white at the voxel size of the acquisition and correlated on the 0.5 mm grid: the SNR is that of the acquisition, and on the fine grid the standard deviation is 0.82 of it for a 1 mm acquisition, with a correlation of 0.69 between neighbours.

## What was checked

Noise, each case against the same acquisition without noise, in voxels of one pure tissue (0.75 mm MPRAGE, SNR 25):

| Case | Nominal SD | SD in WM | Mean shift in WM (Rician theory) | SNR in WM | Correlation of neighbours |
| --- | --- | --- | --- | --- | --- |
| No bias field | 20.04 | 20.02 | +0.42 (+0.40) | 25.0 | 0.001 |
| Bias field 45%, strength 4 | 20.04 | 20.02 | +0.42 (+0.40) | 24.9 | 0.001 |
| Contrast 0.5 | 26.16 | 26.14 | +0.55 (+0.52) | 25.0 | 0.001 |
| Contrast 1.5 | 15.20 | 15.18 | +0.32 (+0.30) | 25.0 | 0.001 |
| Acquired at 1 mm | 20.04 | 20.12 | +0.38 (+0.40) | 24.9 | 0.001 |
| 1 mm, back on 0.5 mm | 20.04 | 16.41 | +0.38 | 30.5 | 0.69 |
| 5% Gaussian | 25.05 | 25.04 | +0.03 (0) | 20.0 | 0.001 |

GM and CSF have the same SD as WM in every case (19.96 to 20.06 without contrast change), and the background has the Rayleigh floor (mean 26.1, SD 13.6 for an expected 25.1 and 13.1 above a clean background of 6.7).

Against `matlab/mri_simulate.m`:

| Step | Result |
| --- | --- |
| Simulated bias field | The same field for the same seed: largest difference 8e-8 for three settings |
| Movement, 17 cases, same pose, without the recentring | Largest difference 3e-7 to 5e-7, rotations included |
| Ringing, 11 cases | Largest difference 3e-7 to 7e-7 |
| Noise | Differs on purpose, see below |

Movement, the distance between the image and its ground truth (the translation that fits the unmoved image best):

| Severity | Head-sized object, without / with recentring | Real head at 1 mm, without / with |
| --- | --- | --- |
| 1 | 0.47 / 0.005 voxels | 0.40 / 0.11 mm |
| 2 | 1.37 / 0.03 voxels | 0.85 / 0.21 mm |
| 3 | 1.95 / 0.04 voxels | 1.56 / 0.46 mm |

The object was measured for four seeds each, with a fit that finds no better shift than none afterwards. The real head is one seed, measured with a simpler fit that reads 0.06 to 0.24 mm for ringing alone, which does not move the image at all; its values with recentring are at that floor for severity 1 and 2. A rotation of the apparent object is not removed.

Two properties of the MATLAB version were not carried over, because they are wrong for the image:

- **A ceiling on the noisy image.** `mri_simulate.m` limits the noisy image to the brightest voxel of the image before the noise. In its own output without a bias field (SNR 25) white matter has a noise SD of 11.8 instead of 19.5 and is 6.5 darker, because white matter is the brightest tissue there and the upper half of its noise is cut off. GM and CSF are not affected (19.5).
- **A limit at the 98th percentile.** It limits the synthesised image to its 98th percentile, which on a whole head at 0.5 mm is the intensity of white matter: more than 2% of the voxels are pure white matter. The scalp fat, 1.4 to 3 times as bright in the scans this was measured on, is cut down to it, and that is what makes white matter the brightest tissue. `t1sim` limits at the 99.99th percentile.
