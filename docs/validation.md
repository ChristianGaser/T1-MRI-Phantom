# Agreement with the MATLAB tool

Part of the [documentation of t1sim](../README.md#documentation).

`matlab/mri_simulate.m` makes its ground truth by the same steps as `t1sim`, with SPM's unified segmentation and CAT12's LAS in place of CAT-Surface's and T1Prep's. This page compares the two labels on one scan. It is a comparison of two implementations of one idea, not a test against a truth.

- [What was compared](#what-was-compared)
- [Result](#result)
- [The intensity of pure GM](#the-intensity-of-pure-gm)
- [Limits](#limits)

## What was compared

| Item | Value |
| --- | --- |
| Scan | One 0.75 mm MPRAGE, whole head |
| MATLAB label | `mri_simulate.m` 0.10.2, before its cleanup of vessels and dura |
| `t1sim` label | `frontend.segment`, on the working grid |
| Common grid | The grid of the scan: the `t1sim` label was resampled onto it once, by linear interpolation |

`mri_simulate.m` writes its label after the cleanup only. The label before it was recovered from the simulation without noise and bias field, which is that label rendered through three tissue intensities and can be inverted. It agrees with the written label to within 0.1 in 98.8% of the brain voxels (RMSE 0.004 there); the remaining 19.7 ml are what the cleanup relabelled.

| Measure | Meaning |
| --- | --- |
| Soft Dice | Overlap of the fraction maps of a tissue: twice the sum of the voxel-wise minimum, over the sum of both maps |
| Boundary shift | Mean displacement of a boundary in mm: the difference in the volume beyond it, divided by the area of the boundary. Positive: the boundary of `t1sim` lies further out |
| Volume | Sum of the fraction map of a tissue |

## Result

| Measure | `t1sim` against MATLAB | Interpolation alone |
| --- | --- | --- |
| Soft Dice CSF | 0.820 | 0.913 |
| Soft Dice GM | 0.908 | 0.949 |
| Soft Dice WM | 0.948 | 0.968 |
| Shift of the GM/WM boundary | +0.04 mm | -0.05 mm |
| Shift of the pial boundary | -0.01 mm | +0.01 mm |

| Volume | `t1sim` | MATLAB |
| --- | --- | --- |
| CSF | 340 ml | 348 ml |
| GM | 695 ml | 708 ml |
| WM | 587 ml | 577 ml |

"Interpolation alone" is the MATLAB label resampled onto the working grid and back, compared with itself. It shows what two linear interpolations cost and so bounds what the first column can reach; the `t1sim` label was interpolated once.

## The intensity of pure GM

The label depends on the intensity that LAS takes for pure GM. `mri_simulate.m` hands LAS the GM mean of SPM's mixture of Gaussians. T1Prep's LAS estimates the peak itself unless it is given one, and on this scan the two differ:

| Source | Pure GM, relative to WM |
| --- | --- |
| SPM's mixture (`mri_simulate.m`) | 0.760 |
| CAT-Surface's mixture (`t1sim`) | 0.758 |
| T1Prep's LAS, own estimate | 0.72 |

The two unified segmentations agree; the difference lies in LAS. Its effect on the label, each against the MATLAB label:

| Intensities handed to LAS | Soft Dice CSF / GM / WM | GM/WM boundary | Pial boundary | Volumes CSF / GM / WM |
| --- | --- | --- | --- | --- |
| None: LAS estimates all three | 0.804 / 0.895 / 0.935 | +0.18 mm | +0.12 mm | 302 / 692 / 627 ml |
| GM from the mixture (what `t1sim` does) | 0.820 / 0.908 / 0.948 | +0.04 mm | -0.01 mm | 340 / 695 / 587 ml |
| CSF, GM and WM from the mixture | 0.822 / 0.910 / 0.948 | +0.06 mm | +0.06 mm | 321 / 708 / 593 ml |

With the GM mean of the mixture, pure GM sits at 379.4 in the corrected scan, against a mixture mean of 380.0; with LAS's own estimate it sat at 358.4. Rendering the fractions with the tissue intensities of the mixture then differs from the corrected scan inside the brain by -0.6% of the WM intensity in the median (RMS 3.8%).

For the peak to hold, it has to be used twice in LAS: the global mapping runs through it, and the local field of the tissue is scaled to it. With the global mapping alone the GM/WM boundary was still 0.14 mm off, because the local field is fitted on the voxels a first label map calls GM and follows their median.

## Limits

- **One scan.** Nothing here says how the numbers vary between scans, scanners or field strengths.
- **No truth.** The MATLAB label is a second implementation. Where the two agree, both may be wrong in the same way.
- **Before the cleanup.** `t1sim` has no cleanup of vessels and dura yet, so the labels were compared without it.
- **Not reproducible from the repository.** The scan and the scripts of this comparison are not part of it.
