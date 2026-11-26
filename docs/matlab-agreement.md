# Agreement with the MATLAB tools

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim` was ported from the [MATLAB tools](../matlab/README.md). This page holds everything about the relation of the two: the ground truth they make of one scan, the steps that give the same result, and the steps where `t1sim` differs on purpose. It is for those who know or use the MATLAB tools; nothing here is needed to use `t1sim`.

The comparison of the ground truth is one of two implementations of one idea, not a test against a truth: `mri_simulate.m` uses SPM's unified segmentation and CAT12's LAS, `t1sim` CAT-Surface's and T1Prep's.

- [Summary](#summary)
- [The ground truth of one scan](#the-ground-truth-of-one-scan)
- [The intensity of pure GM](#the-intensity-of-pure-gm)
- [Limits of the comparison](#limits-of-the-comparison)
- [Geometric phantom](#geometric-phantom)
- [Where t1sim differs on purpose](#where-t1sim-differs-on-purpose)

## Summary

| Step | Compared how | Result |
| --- | --- | --- |
| Ground truth before the cleanup | One scan, on the grid of the scan | Soft Dice 0.82 / 0.91 / 0.95 for CSF / GM / WM, boundaries within 0.04 mm; [below](#the-ground-truth-of-one-scan) |
| Cleanup | The voxels both relabel as CSF, same scan | Dice 0.47 inside MATLAB's brain mask; [Cleanup](cleanup.md#measured) |
| Cortical thickness | Not against MATLAB; against the thickness asked for | 1.50 / 2.49 / 3.49 mm for 1.5 / 2.5 / 3.5 mm on a sphere; [Anatomy](anatomy.md#cortical-thickness) |
| Simulated bias field | Same seed, three settings | The same field: largest difference 8e-8 of a field with a mean of 1 |
| Movement | MATLAB's pose fed to `t1sim`, 17 cases, rotations included, without the recentring | Largest difference 3e-7 to 5e-7 of the brightest voxel |
| Ringing | 11 cases | Largest difference 3e-7 to 7e-7 of the brightest voxel |
| Geometric phantom | Voxel by voxel, 29 phantoms | Identical in 27; [below](#geometric-phantom) |
| Atrophy, voxel size, contrast, noise, seeds, output | - | Differ on purpose; [below](#where-t1sim-differs-on-purpose) |

## The ground truth of one scan

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

## Limits of the comparison

- **One scan.** Nothing here says how the numbers vary between scans, scanners or field strengths.
- **No truth.** The MATLAB label is a second implementation. Where the two agree, both may be wrong in the same way.
- **Before the cleanup.** The labels were compared before the cleanup of either tool; the two cleanups are compared in [Cleaning the ground truth](cleanup.md#measured).
- **Volumes are those of the scan.** Everything on this page was measured on the grid of the scan. What `t1sim` writes lies in the space of the template, where this brain is 1.26 times as large; see [The working grid](frontend.md#the-working-grid).
- **Not reproducible from the repository.** The scan and the scripts of this comparison are not part of it.

## Geometric phantom

| Compared | Result |
| --- | --- |
| Label map and image, 29 phantoms of 20 parameter sets | Identical in 27 |
| The other two | One partial volume step (1/15) in 96 of 262144 and 536 of 591360 voxels |
| Fine-grid tissue counts | Identical in all 29 |
| Evaluation, 11 phantoms | Every statistic within 2e-6 mm |
| Geometry columns of the [table of results](phantom.md#results) | Those of the MATLAB documentation, to the digits it gives |

The two exceptions have voxels whose distance ties exactly with a shifted level; CAT12's `cat_bwdist` is one float32 step off the correctly rounded distance in 15 to 20% of the samples, and that bit decides the tie. Fed with MATLAB's distance map, the port gives MATLAB's label exactly.

## Where t1sim differs on purpose

| Step | `mri_simulate.m` | `t1sim` | Why |
| --- | --- | --- | --- |
| Noise | The noisy image is limited to the brightest voxel of the image before the noise | No upper limit | The limit cuts off the upper half of the noise in the brightest tissue. In MATLAB's output without a bias field, white matter has a noise SD of 11.8 instead of 19.5 |
| Synthesis | Limited to the 98th percentile of the image | Limited to the 99.99th | On a whole head at 0.5 mm the 98th percentile is the intensity of white matter, so the scalp fat is cut down to it |
| SNR after a contrast change | Refers to white matter before the change | Refers to white matter in the image that is written | Otherwise the SNR that is asked for is not the SNR of the image |
| Voxel size | Sinc interpolation of the image, linear interpolation of the label | k-space cut at the new Nyquist frequency; fractions averaged over the box of a voxel | Interpolation at the new voxel centres aliases and samples the label at a point |
| Atrophy | CSF fraction scaled inside the region | A deformation: thinning, or a change of volume | See [Which kind of atrophy](anatomy.md#which-kind-of-atrophy) |
| Seeds | One noise pattern for all images by default | Every output file its own | The same pattern in every image can be learned by a network |
| Movement | The mean pose is subtracted | The image is also brought back onto the object | The magnitude image was left 0.4 to 1.6 mm from its ground truth on a real head |
| Output grid | The grid of the scan, or an affinely registered one | Always affinely registered at 0.5 mm | The segmentation of the front end works there |
| Vessels | CAT12's atlas partitioning | T1Prep's detection and voxels brighter than white matter | The partitioning is not available in Python |
