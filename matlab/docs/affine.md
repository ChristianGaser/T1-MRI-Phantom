# Affine registration of the output

Part of the [documentation of the MATLAB tools](../README.md).

`simu.affine` writes the simulated image and its ground truth onto one common
grid instead of the grid of the input, which is what a training set needs. The
registration is applied **after** all anatomical modifications and **before**
the synthesis, so the T1w image is generated from the already resampled label:
the image stays an exact function of the label it ships with, which a pipeline
that resamples an image and segments it a second time cannot offer.

```matlab
% the grid the deepmriprep training scripts use: 336x384x336 at 0.5mm
simu = struct('name', 'sub-01_T1w.nii', 'snrWM', 40, 'affine', 1);
mri_simulate(simu, struct('percent', 0));
% -> sub-01_space-MNI152_res-0p5mm_desc-snr40_T1w.nii
%    sub-01_space-MNI152_res-0p5mm_dseg.nii
```

Field | Meaning (Default)
------|------------------
method | `'deformation'` fits an affine to the overall deformation field of the segmentation, `'spm'` takes `res.Affine` as it is (`'deformation'`)
mask | Region the fit is restricted to: `'brain'`, `'nonbrain'`, `'head'` or `'all'` (`'brain'`)
grid | `'deepmriprep'` for 336x384x336 voxels of 0.5mm on the grid its training scripts use, or `'custom'` through `dim`/`mat` or `bb`/`vx` (`'deepmriprep'`)
dim, mat | Dimensions and voxel-to-mm matrix of a custom grid, the matrix one based as every SPM matrix is
bb, vx | Bounding box in mm and voxel size, as an alternative to `dim`/`mat` (`vx` defaults to 1mm)
interp | Interpolation degree passed to `spm_slice_vol`: `-5` is sinc, a positive value a b-spline (`-5`)
space | Label of the BIDS `space` entity (`'MNI152'`)

## The deepmriprep grid

`2_prep_segment.py` samples the 339x411x339 grid of deepmriprep's
`Template_05mm_bet.nii.gz`, whose origin is (-84,-120,-72)mm, and crops it with
`[1:-2, 15:-12, :336]` to 336x384x336. The crop starts at the voxels (1,15,0),
so the origin of the training grid is that of the template shifted by
(0.5,7.5,0)mm, i.e. (-83.5,-112.5,-72)mm.

The nibabel matrix that `2_prep_segment.py` attaches to the cropped volume still
carries the uncropped origin and has `[0 0 0 0]` as its last row, so it cannot
be used to define the space; the template plus the crop can, and that is what
the voxel content follows. Getting this wrong puts the simulation 7.5mm off
along the anterior-posterior axis, which is 15 voxels.

## Why the affine is fitted to the deformation field

`res.Affine` only *initialises* the unified segmentation. When that initial
registration fails, which happens for an unusual field of view or a strong
tilt, the subsequent nonlinear stage compensates for a part of the error, so
the composed deformation ends up aligned while its affine part does not. The
default therefore fits an affine to the composed field by least squares, which
recovers the global alignment the segmentation actually converged to. The fit
is a closed-form linear problem on at most 200000 sampled voxels and costs
nothing next to the segmentation itself.

The `mask` decides what the affine is optimal *for*. `'brain'` aligns the
brains, which is what a training grid for brain segmentation needs.
`'nonbrain'` fits outside the brain, where the deformation is closer to affine
already, and reproduces the convention of the deepmriprep preprocessing.

## Interpolation and the label

The label is resampled as the scalar PVE image, not as three separate tissue
fractions, and the fractions are rebuilt from it afterwards. Interpolating the
three volumes on their own lets a voxel become a CSF/WM mixture without any
GM, which is exactly the configuration a scalar label cannot represent and
that shows up as spurious GM at a ventricle border.

Sinc interpolation rings at the hard edges of a skull-stripped label: measured
on a 0/3 step it overshoots to about **-0.33 and +3.30**. The resampled label
is therefore clamped back into `[0,3]` (`[0,4]` when WMHs are simulated) before
the fractions are derived from it.

Below a label value of 1, the triangular decomposition describes the fade from
the background into the CSF and its fractions sum to less than one. That sum is
kept rather than normalized to one, so the synthesis blends into the bias
corrected image there, as it does at the brain boundary anyway. Normalizing
would snap the whole ramp to full CSF and dilate the brain by half a voxel.

## Notes

- The phase-encoding axis of `motion` and `ringing` is resolved on the
  orientation of the *input*, while the artefact is applied on the output
  grid. Name `pe` as a voxel axis (`1`, `2`, `3`) if a rotation between the two
  matters for what is being tested.
- The bias corrected image is carried to the new grid and handed to the
  synthesis instead of being rebuilt there. Its DCT basis is defined over the
  field of view of the input and must not be re-evaluated on a grid that is not
  a resampling of that same box.
- The JSON sidecar records the fitted affine, the grid and the interpolation
  under `SimulationParameters.AffineRegistration`.
