# Cleaning the ground truth

Part of the [documentation of the MATLAB tools](../README.md).

`simu.clean` removes from the **label** what is not the tissue it looks like,
and deliberately leaves it in the **image**.

```matlab
simu = struct('name','sub-01_T1w.nii', 'snrWM',40, 'affine',1, 'clean',1);
mri_simulate(simu, struct('percent',0));
% -> sub-01_space-MNI152_res-0p5mm_desc-snr40_T1w.nii        (vessel still bright)
%    sub-01_space-MNI152_res-0p5mm_desc-Clean_dseg.nii       (vessel labelled CSF)
%    sub-01_space-MNI152_res-0p5mm_desc-Clean_label-GM_probseg.nii
```

## Why the image is not cleaned too

The simulated image is `sum_k mn_k * Yp0toC(label,k)`, so removing a vessel
from the label would remove it from the image as well, and the result would be
cleaner than any real scan. Worse, a synthesis whose label describes every
structure of its own image is related to it by a known one-dimensional curve:
the task is invertible, and a network trained on it learns to invert a lookup
table rather than to segment.

Rendering from the uncorrected fractions and writing the corrected ones breaks
that. The T1w still shows the vessel at its true intensity while the ground
truth calls it CSF, which is exactly the "looks like GM but is not GM" signal
that a segmentation has to learn and that a self-consistent synthesis cannot
provide.

## Blood vessels and dura

`cat_vol_partvol` already detects them and writes them into its region label as
`LAB.BV`, using a prior built from the MRA scans of IXI and ICBM
(`cat_bloodvessels.nii`) together with the divergence and the gradient of the
image and a region growing. Nothing is detected again here, the label is only
read out, so the result is CAT12's own and not a second opinion.

A Hessian **sheetness** filter was considered instead and rejected. The
cortical ribbon is itself a sheet of 2-3mm, so such a filter fires on the
cortex as hard as it does on the dura, and at 0.5mm the dura is one to two
voxels and sits at the noise floor. Frangi vesselness is better posed for the
vessels, a tube being separable from a sheet, but the large pial vessels lie
inside the sulcal CSF touching the cortex and merge with it. What separates the
two is not the shape but the position relative to the WM and the anatomical
prior, which is what CAT12 uses and what the divergence terms of
`cat_vol_partvol` already contribute.

## Dura

`clean.dura` is a distance in mm and removes what `clean.bv` does not: the
dural membranes along the inner skull, the falx and the tentorium, which are
sheets of GM intensity that no vessel prior covers.

They cannot be separated from the cortex by their distance to the WM. Measured
on a 0.5mm simulation, 99.5% of all GM labelled voxels lie within 5.3mm of the
WM and the dura is among them, because the dura over a gyral crown is as close
to the WM as the cortex underneath it.

What separates them is that the cortex hangs on the WM while the dura sits on
the far side of the subarachnoid CSF, so a path from the WM to the dura has to
leave the tissue. `cat_vbdist` measures a distance *inside a mask*, and with the
mask set to the tissue the dura ends up far away or unreachable while the
cortex stays within a few mm. It is the same construction `cat_vol_partvol`
uses for its distant blood vessels, where it compares a constrained against a
free distance, and it costs about 1.5s on a 0.5mm volume.

Measured on the same subject, with the threshold in mm and the volume it
removes:

`dura` | removed | comment
-------|---------|--------
3 | 81 cm³ | too much, eats into the cortex
4 | 23 cm³ | default, the fullest coverage of the dural band and the falx
5 | 10 cm³ | matches a true geodesic reference at 96% precision
0 | - | off

Lower is more aggressive. The cerebellum is protected because its WM is a thin
branched tree and its folia are legitimately far from it along the tissue; it
is not dilated, so the tentorium outside it stays detectable. Basal ganglia,
thalamus, hippocampus and brainstem are protected with a 2mm dilation.

Note that this only relabels tissue that the label calls GM. The broad mantle
of subarachnoid space between the cortex and the inner skull is already CSF in
the label, and how much of it is kept at all is decided by the APRG
skull-stripping, not here.

## Periventricular CSF/WM partial volume

A voxel that mixes CSF and WM has an intensity between the two, i.e. the
intensity of GM, so a label built from the intensity calls it GM. This is a
port of the level 3 cleanup of `cat_main_cleanup`: the voxel has to lie next to
the ventricle, the brainstem or the corpus callosum, must not be in the basal
ganglia or the thalamus, must sit between pure WM and pure CSF, and must belong
to a thin structure rather than to a real band of GM.

The correction sets the GM fraction to zero and re-decomposes the voxel as a
pure CSF/WM mixture. With `1*c + 3*(1-c) = label` the CSF fraction is
`c = (3-label)/2`, so **the label comes out exactly as it went in** and only
the fractions change. That is precisely why a scalar label cannot express this
correction, and why the fractions are written as `_probseg` files.

On a phantom with a ventricle rim that reads as GM and a genuine cortical band
of the same label value, the detection puts 100% of its voxels on the rim and
none on the cortex, at 100% rim coverage.

Two deviations from `cat_main_cleanup`: its `Yp0` comes from quantized uint8
posteriors so it can test `Yp0==3` and `Yp0==1`, while the label here is
continuous and thresholds are used instead; and its radii are partly in voxels,
which would silently halve them for a 0.5mm image, so the radii here are in mm
with the voxel size passed to `cat_vol_morph`.

## Cost and notes

- It needs the atlas partitioning of `cat_vol_partvol`, which costs minutes.
  It is computed once and shared with `closeWMHholes` when both are used.
- The label gains a `Clean` tag in its `desc` entity. The simulated image does
  not, because it is bit for bit the same with and without the option.
- `_probseg` files are `uint8` with a scale of 1/255, like CAT12's `p1`/`p2`/`p3`.
- The GM fraction is what a nogm-style correction model needs: the voxels where
  the triangular decomposition of the label overestimates GM are exactly
  `Yp0toC(label,2) - probseg_GM`.
