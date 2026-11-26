# Geometric thickness phantom

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim.phantom` builds an object whose cortical thickness is known by construction: a ball of white matter with regular folds, a band of grey matter of constant thickness around it, and CSF around that. It has no artefacts and no random part, so the error a thickness measure shows on it is the error of the measure, down to the [accuracy of the phantom itself](#how-exact-it-is), about 0.05 mm at 0.5 mm voxels. It is a port of `matlab/thickness_phantom.m` and `matlab/thickness_phantom_eval.m`.

- [Usage](#usage)
- [Why the thickness is exact](#why-the-thickness-is-exact)
- [Evaluation](#evaluation)
- [Results](#results)
- [How exact it is](#how-exact-it-is)
- [Adding artefacts](#adding-artefacts)

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
| `wm_offset` | Level of the distance map that is the white matter boundary, in mm (`None`: half a voxel). It has to stay above the sampling of the distance map; a larger value makes the phantom more exact, see [below](#how-exact-it-is) |
| `supersample` | Factor of the finer grid the object is built on (`3`) |
| `pve_steps`, `pve_range` | Number and range of the sub-voxel offsets that are averaged into partial volumes (`15`, `1`) |
| `intensities` | Intensity of CSF, GM and WM in the image (`(30, 100, 150)`) |
| `outdir` | Write `_dseg` and `_T1w` files with JSON sidecars there (`None`: write nothing) |

`thickness_phantom` returns the label map (CSF = 1, GM = 2, WM = 3 with partial volumes), the image that belongs to it, and the geometry as a dict, which also holds the tissue volumes of the supersampled grid (`volume_ref`) and of the label map (`volume_pve`). The files are named `thicknessPhantom_res-0p5mm_desc-Thickness25mm_dseg.nii.gz` and `..._T1w.nii.gz`; the thickness in the name is rounded to a tenth of a mm, and a series in which two values get the same name is refused.

## Why the thickness is exact

The cortex is not made by enlarging the radius of the folded ball. A radial offset of `t` gives a thickness of `t * cos(alpha)` along the normal, where `alpha` is the angle between the radius and the normal, so the cortex would be too thin on the flanks of the folds. Instead the Euclidean distance to the white matter is computed once, and both boundaries of the cortex are level sets of that one distance map: the white matter boundary at `wm_offset` and the pial boundary at `wm_offset + thickness`. Two level sets of a distance map are `thickness` apart everywhere, whatever the folds look like.

Partial volumes come from shifting both boundaries by sub-voxel offsets of the distance, making hard labels for each and averaging them. For a boundary that is locally flat this gives the exact linear ramp.

## Evaluation

`evaluate_phantom` checks two independent things:

1. **The phantom.** The distance between the two boundaries is recovered from the label map alone: it is supersampled, the distance from the white matter surface is computed, and that distance is read at the sub-voxel crossings of the pial boundary along the three axes. Tissue volumes are compared with the supersampled grid the phantom was built on. It also reports how much of the white matter surface faces a buried sulcus, a fold narrower than twice the thickness, in which any thickness measure overestimates.
2. **A thickness measure.** By default the projection-based thickness of cat-surf (`cat_surf.vol_thickness_pbt`); `thickness_map=` takes any other map, and `pbt=False` skips it. `pbt_options` are passed on to cat-surf: its default adds 0.25 mm to every thickness (`correct_thickness`), a correction for the label maps of its own pipeline that the phantom does not need.

```python
result = phantom.evaluate_phantom(label, info, pbt_options={"correct_thickness": 0.0})
result["geom"]["outer"]["bias"]     # the phantom: recovered thickness minus the known one, mm
result["geom"]["buried"]            # share of the white matter surface that faces a sulcus without CSF
result["measure"]["gm"]["bias"]     # the measure over all grey matter voxels, mm
result["volume"]                    # tissue volumes of the label map and of the reference grid, ml
```

Every entry of `geom` and `measure` holds `n`, `mean`, `median`, `sd`, `bias`, `rmse`, `p5` and `p95`.

## Results

Default phantom (radius 22 mm, amplitude 2.5 mm, wavelength 12 mm) at 0.5 mm. "Geometry" is the thickness recovered from the label map; the two last columns are cat-surf's projection-based thickness over the grey matter voxels.

| Asked for | Geometry, bias | SD | Buried sulci | PBT bias, default | PBT bias, `correct_thickness=0` |
| --- | --- | --- | --- | --- | --- |
| 1.5 mm | -0.014 | 0.036 | 0.0% | +0.256 | +0.006 |
| 2.0 mm | -0.017 | 0.037 | 0.3% | +0.238 | -0.012 |
| 2.5 mm | -0.018 | 0.036 | 4.2% | +0.249 | -0.001 |
| 3.0 mm | -0.016 | 0.039 | 6.9% | +0.230 | -0.020 |
| 3.5 mm | -0.018 | 0.035 | 10.3% | +0.267 | +0.017 |

- cat-surf's projection-based thickness adds 0.25 mm by default. Without that, its mean over the grey matter voxels is within 0.02 mm of the thickness asked for, which is less than the accuracy of the phantom itself (next section).
- Tissue volumes agree with the supersampled reference within 0.43% for CSF and WM and within 1.46% for GM.

## How exact it is

The two level sets of the distance map are exactly `thickness` apart. What limits the phantom is the sampling of that map. Two effects were measured, both on a ball without folds, where the right answer can be written down, and both are in the MATLAB version as well:

- **The cortex is thicker than asked for**: by 0.039 mm at 0.5 mm voxels and by 0.07 mm at 1 mm, going by the tissue volumes of the label map. With `wm_offset=1.0` it is 0.010 mm at 0.5 mm.
- **The readout of the evaluation is too small by about 0.4 samples of its grid**: -0.067 mm at 0.5 mm with a supersampling of 3, and -0.118 / -0.082 / -0.065 mm at 1 mm with 3 / 5 / 7, on a label map written down analytically without any distance map. The half-sample correction it applies is right only for surfaces along the grid.

Together they give -0.03 mm on the ball; the folded phantom reads -0.014 to -0.018 mm in the table above, so the folds change them by about 0.01 mm. The geometry bias of the table is therefore two errors that largely cancel and not the accuracy of the phantom: take the thickness of the default phantom as right to about 0.05 mm at 0.5 mm voxels, and a difference between two thickness measures below that as not resolved by it.

## Adding artefacts

The phantom has no noise and no bias field. `simulate.acquire` takes the ground truth of a scan and not a phantom, but the artefacts are plain functions of an image:

```python
from t1sim import artefacts

white = info["intensities"][2]
field = artefacts.simulated_bias_field(t1w.shape, info["vx"], percent=20, seed=1)
image = artefacts.add_noise(t1w * field, sigma=white / 30, seed=1)      # SNR 30 in white matter
```

The label map stays the ground truth of that image. Movement and ringing are in `t1sim.kspace`; see [Artefacts of the acquisition](artefacts.md).

The default phantom takes 6 to 8 s and its evaluation about 15 s. The label maps of 27 of 29 phantoms are identical to those of the MATLAB tool, voxel by voxel; see [Agreement with the MATLAB tools](matlab-agreement.md#geometric-phantom).
