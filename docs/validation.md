# Validating a method

Part of the [documentation of t1sim](../README.md#documentation).

A simulated image comes with the tissue fractions it was made from. This page says how to use the two to test a segmentation or a cortical thickness method, and how far the ground truth can be trusted. For the other use of `t1sim` see [Training data for deepmriprep](training-data.md).

- [What to compare](#what-to-compare)
- [Settings for a validation set](#settings-for-a-validation-set)
- [How exact the ground truth is](#how-exact-the-ground-truth-is)
- [What a result does not show](#what-a-result-does-not-show)

## What to compare

| Question | Ground truth | Where |
| --- | --- | --- |
| Tissue segmentation | The three fractions, `..._label-GM_probseg.nii.gz` and so on; or the label map `..._dseg.nii.gz` for a method that writes one | [What is written](simulation.md#what-is-written) |
| Tissue volumes | The sum of a fraction map times the voxel volume | below |
| Robustness against noise, bias field, contrast, voxel size, movement | The same ground truth under every acquisition: only the image changes | [Artefacts](artefacts.md) |
| Cortical thickness in a real anatomy | The thickness asked for with `thickness=` | [Cortical thickness](anatomy.md#cortical-thickness) |
| Cortical thickness itself, without an anatomy | The geometric phantom | [Phantom](phantom.md) |
| Sensitivity to a local loss of tissue | Two simulations of one scan, with and without `atrophy=`: the difference of their fractions, the Jacobian determinant, and the numbers of the sidecar | [Atrophy](anatomy.md#atrophy) |

Image and ground truth lie on one grid and share their affine, so they can be compared voxel by voxel. A method that works on the grid of its input needs no resampling; one that resamples (to a template, to 1 mm) has to bring its result back, or the ground truth to it, and that interpolation is then part of what is measured.

**Volumes** on the grid of the scan and on the grid of an anatomy made on the scan, the defaults, are those of the subject, and so are the numbers of the anatomy in the sidecar (`VolumeMl`, `TissueChangeMl`). In template space they are those of the template ([The grid of the anatomy](frontend.md#the-grid-of-the-anatomy)); there, too, nothing has to be converted to compare with what a method reports for the simulated image, because the method sees the same enlarged head. `AffineVolumeScaling` in the sidecar is the factor between the two.

```python
import nibabel as nib
import numpy as np

gm = nib.load("sub-01_desc-Clean_label-GM_probseg.nii.gz")
voxel_ml = abs(np.linalg.det(gm.affine[:3, :3])) / 1000
volume_ml = gm.get_fdata().sum() * voxel_ml          # of the simulated image
```

## Settings for a validation set

```python
from t1sim import simulate

truth = simulate.ground_truth("sub-01_T1w.nii.gz")
for snr in (15, 30, 60):
    for bias in (0, 30, 60):
        simulate.write(truth, simulate.acquire(truth, snr=snr, bias_percent=bias, seed=1), out_dir="validation")
```

| Point | Setting | Why |
| --- | --- | --- |
| Grid | `grid="scan"` (the default), or `grid="anatomy"` | On the grid of the scan the image has the voxels, the field of view and the size of a real scan of that subject, which is what a pipeline is built for; its ground truth is the mean of the fractions over those voxels. On the grid of the anatomy the ground truth is at 0.5 mm, as it was made, with nothing averaged. Both are in the space of the subject |
| Seed | A fixed `seed` | The images of a series then differ in what is varied and not in the noise pattern. For training data the opposite holds |
| Voxel size | `resolution=1` | Image and ground truth at 1 mm, along the axes of the scan, the ground truth as the average of the fractions over each voxel. On the grid of the anatomy `upsample=True` gives a 1 mm scan resampled to 0.5 mm instead |
| Cleanup | `clean=True` (the default) | Without it, vessels and dura count as grey matter in the ground truth |
| Bias field | The default `(2, 0)`, or the stored fields `"A"`, `"B"`, `"C"` | One field for all scans makes the scans comparable |
| Folder | One `out_dir` per kind of run | `write` refuses to replace a ground truth made with other settings, see [Names](simulation.md#names) |

Do not validate on scans a network was trained on, real or simulated.

## How exact the ground truth is

Exact, by construction:

- **The image without artefacts.** Inside the brain every voxel is the sum of the tissue intensities weighted by the fractions it is rendered from. Before the cleanup those fractions are the ground truth.
- **Contrast.** A change of the contrast changes the tissue intensities, and the image is mixed again; a voxel of two tissues stays their mixture.
- **Bias field and noise.** Neither moves anything. The bias field has a mean of 1 over the volume, and the noise has the standard deviation that the sidecar states (checked to 0.1%, see [Artefacts](artefacts.md#what-was-checked)).

Not exact, and by how much:

| Source | Size | Why |
| --- | --- | --- |
| Cleanup | 29 ml and 23 ml of the ground truth on two scans differ from what the image shows | On purpose: vessels, dura and ventricle rims look like grey matter. Use `clean=False` to test a method against the image alone |
| Voxel size | On a ball with sharp boundaries acquired at 1 mm, the image differed from the mixture of the averaged fractions by 6.5% of the GM/WM contrast (RMS) and by up to 20% next to a boundary; pure white matter rang by -2 to +4% | The image has the point spread function of an acquisition (a cut of k-space, Gibbs ringing included), the ground truth the average over the box of a voxel. In a scanner it is the same, and the fractions stay the right answer to "how much of this voxel is grey matter" |
| Movement | The image lies 0.1 to 0.5 mm from its ground truth for severity 1 to 3 | The offset of the magnitude image is estimated and removed, to that residue; a rotation of the apparent object is not removed |
| Cortical thickness | 0.01 mm along rays on a sphere; the partial volume at its boundaries is 1.15 mm wide on the 0.5 mm grid | The label map is smoothed before the cortex is rebuilt, so the ground truth is smoother than its grid |
| Atrophy | For a loss of 10%, three ways of counting what was reached give -10.0, -9.7 and -11.1% with the thinning, and -10.0, -6.1 and -7.6% with the change of volume | Tissue crosses the border of the region, and its neighbours are deformed a little as well; see [what `change` means](anatomy.md#what-change-means) |
| Geometric phantom | 0.05 mm at 0.5 mm voxels | The sampling of its distance map |
| Grid of the scan | Exact from an anatomy made on the scan, the default: the grids share their axes. Between grids turned against each other a fraction is within 0.005 of the exact mean over its voxel up to 11 degrees, 0.01 up to 29 degrees, also for voxels three or four times as long as wide, and the volumes agree to five digits | A voxel of the scan is a parallelepiped in a grid turned against it, which is covered by small boxes; see [Voxel size](artefacts.md#voxel-size) |
| Stored bit depth | 1/255 of a fraction | 8 bit |

## What a result does not show

- **The ground truth is defined by intensity.** It is what the front end read from the scan: a structure brighter than cortex, such as the pallidum or the thalamus, is a mixture of GM and WM in it. A method that knows the anatomy and labels the thalamus as grey matter is "wrong" against this ground truth, and a method that reads intensities like the front end is favoured by it.
- **The image is easier than a scan in some ways.** Inside the brain it has three intensities and their mixtures. The local variation of the tissue intensities that a real scan has within a tissue is gone, except for what the fractions carry.
- **The anatomy is that of one subject.** By default it is made in the size of that subject, and a thickness asked for in mm is one of the subject. With `space="template"` the head is scaled to the template, and so is the thickness asked for.
- **One or two scans stand behind the numbers of this documentation.** They say how `t1sim` behaves there.
