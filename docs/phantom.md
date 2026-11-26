# Geometric thickness phantom

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim.phantom` builds an object whose cortical thickness is exact by construction: a ball of white matter with regular folds, a band of grey matter of constant thickness around it, and CSF around that. It has no artefacts and no random part, so the error a thickness measure shows on it is the error of the measure. It is a port of `matlab/thickness_phantom.m` and `matlab/thickness_phantom_eval.m`.

- [Usage](#usage)
- [Why the thickness is exact](#why-the-thickness-is-exact)
- [Evaluation](#evaluation)
- [Results](#results)
- [How exact it is](#how-exact-it-is)
- [Agreement with the MATLAB tool](#agreement-with-the-matlab-tool)

## Usage

```python
from t1sim import phantom

# one phantom with 2.5 mm thickness at 0.5 mm, then evaluate it
label, t1w, info = phantom.thickness_phantom(2.5)
result = phantom.evaluate_phantom(label, info)

# a series; the phantoms share one distance map
labels, t1ws, infos = phantom.thickness_phantom([1.5, 2.0, 2.5, 3.0, 3.5], outdir="phantoms")
results = phantom.evaluate_phantom(labels, infos, fig="phantom_qc.png")
```

| Argument | Meaning (default) |
| --- | --- |
| `thickness` | Thickness in mm, one value or a sequence (`2.5`) |
| `dim`, `vx` | Size and voxel size of the volume (`(128, 128, 128)`, `0.5`) |
| `radius`, `amplitude` | Radius of the white matter ball and amplitude of its folds in mm (`22`, `2.5`) |
| `fold` | `"cartesian"` with a `wavelength` (`12`), or `"spherical"` with `folds` (`(6, 6)`), which get arbitrarily fine towards the poles |
| `csf` | Thickness of the CSF layer in mm (`2`) |
| `supersample` | Factor of the finer grid the object is built on (`3`) |
| `pve_steps`, `pve_range` | Number and range of the sub-voxel offsets that are averaged into partial volumes (`15`, `1`) |
| `intensities` | Intensity of CSF, GM and WM in the image (`(30, 100, 150)`) |
| `outdir` | Write `_dseg` and `_T1w` files with JSON sidecars there (`None`: write nothing) |

`thickness_phantom` returns the label map (CSF = 1, GM = 2, WM = 3 with partial volumes), the image that belongs to it, and the geometry as a dict.

## Why the thickness is exact

The cortex is not made by enlarging the radius of the folded ball. A radial offset of `t` gives a thickness of `t * cos(alpha)` along the normal, where `alpha` is the angle between the radius and the normal, so the cortex would be too thin on the flanks of the folds. Instead the Euclidean distance to the white matter is computed once, and both boundaries of the cortex are level sets of that one distance map: the white matter boundary at `wm_offset` and the pial boundary at `wm_offset + thickness`. Two level sets of a distance map are `thickness` apart everywhere, whatever the folds look like.

Partial volumes come from shifting both boundaries by sub-voxel offsets of the distance, making hard labels for each and averaging them. For a boundary that is locally flat this gives the exact linear ramp.

## Evaluation

`evaluate_phantom` checks two independent things:

1. **The phantom.** The distance between the two boundaries is recovered from the label map alone: it is supersampled, the distance from the white matter surface is computed, and that distance is read at the sub-voxel crossings of the pial boundary along the three axes. Tissue volumes are compared with the supersampled grid the phantom was built on. It also reports how much of the white matter surface faces a buried sulcus, a fold narrower than twice the thickness, in which any thickness measure overestimates.
2. **A thickness measure.** By default the projection-based thickness of cat-surf (`cat_surf.vol_thickness_pbt`); `thickness_map=` takes any other map, and `pbt=False` skips it.

## Results

Default phantom (radius 22 mm, amplitude 2.5 mm, wavelength 12 mm) at 0.5 mm. "Geometry" is the thickness recovered from the label map; the two last columns are cat-surf's projection-based thickness over the grey matter voxels.

| Asked for | Geometry, bias | SD | Buried sulci | PBT bias, default | PBT bias, `correct_thickness=0` |
| --- | --- | --- | --- | --- | --- |
| 1.5 mm | -0.014 | 0.036 | 0.0% | +0.256 | +0.006 |
| 2.0 mm | -0.017 | 0.037 | 0.3% | +0.238 | -0.012 |
| 2.5 mm | -0.018 | 0.036 | 4.2% | +0.249 | -0.001 |
| 3.0 mm | -0.016 | 0.039 | 6.9% | +0.230 | -0.020 |
| 3.5 mm | -0.018 | 0.035 | 10.3% | +0.267 | +0.017 |

- The geometry columns equal those of the MATLAB documentation to the digits it gives.
- cat-surf's projection-based thickness adds 0.25 mm by default (`correct_thickness`), a correction for the label maps of its own pipeline. Without it, it reads the phantom within 0.02 mm.
- Tissue volumes agree with the supersampled reference within 0.43% for CSF and WM and within 1.46% for GM.

## How exact it is

Checked on labels that were written analytically, without any distance map, these two effects are in the phantom and its evaluation, in the MATLAB version as well:

- **The readout of the evaluation is too small by about 0.4 samples of its grid**: -0.067 mm at 0.5 mm with a supersampling of 3, and -0.118 / -0.082 / -0.065 mm at 1 mm with 3 / 5 / 7. The half-sample correction it applies is right only for surfaces along the grid.
- **The cortex of the default phantom is about 0.03 to 0.04 mm thicker than asked for at 0.5 mm**, and 0.07 mm at 1 mm. For a ball without folds the tissue volumes give +0.039 mm; with `wm_offset=1.0` it is +0.010 mm.

The geometry bias of -0.017 mm in the table is therefore the two nearly cancelling, not a measure of the accuracy of the phantom. The thickness of the phantom is right to about 0.04 mm at 0.5 mm.

## Agreement with the MATLAB tool

| Compared | Result |
| --- | --- |
| Label map and image, 29 phantoms of 20 parameter sets | Identical in 27 |
| The other two | One partial volume step (1/15) in 96 of 262144 and 536 of 591360 voxels |
| Fine-grid tissue counts | Identical in all 29 |
| Evaluation, 11 phantoms | Every statistic within 2e-6 mm |

The two exceptions have voxels whose distance ties exactly with a shifted level; CAT12's `cat_bwdist` is one float32 step off the correctly rounded distance in 15 to 20% of the samples, and that bit decides the tie. Fed with MATLAB's distance map, the port gives MATLAB's label exactly.

The default phantom takes 6 to 8 s and its evaluation about 15 s.
