# mri_simulate
Simulates T1-weighted MR images with optional atrophy, cortical thickness control, WMHs, RF B1 inhomogeneities, and noise (Gaussian or Rician at a target WM SNR). Writes JSON sidecars with simulation metadata.

![Pipeline overview](docs/T1-MRI-Phantom-Scheme.png)

## Overview
`mri_simulate` generates a realistic T1-weighted (T1w) image and its explicit ground truth from a high-quality input (e.g., 0.5 mm Colin27 or a custom T1w). Key steps:

- Start from a segmented T1w volume (GM, WM, CSF, and background) using SPM unified segmentation with the bundled Blaiotta head and neck TPM (`BlaiottaTPM.nii`, seven tissue classes) and its batch job `BlaiottaSegmentJob.m`, both stored next to `mri_simulate.m`. The result is cached as `<name>_seg8.mat` next to the input and is only an initialization.
- Locally normalize tissue intensities with CAT12 Local Adaptive Segmentation (LAS), denoise with SANLM and skull-strip with CAT's adaptive probability region-growing (APRG), then scale CSF/GM/WM to canonical values (1/2/3) to obtain a PVE-like label image.
- Optionally close WM holes to remove native WMHs before adding synthetic lesions.
- Insert user-defined anatomical changes: atlas-based atrophy (e.g., Hammers) and probabilistic WMHs.
- Synthesize a new T1w as the probability-weighted mixture of the tissue means (estimated from the SPM Gaussian mixture) using the modified PVE labels and the optional WMH class as weights; everything that is not brain keeps the intensity of the bias-corrected input, and the two blend continuously at the brain boundary. Optionally modulate with RF bias fields, apply a contrast change and add Rician or Gaussian noise.
- Outputs follow BIDS-like naming with JSON sidecars capturing all simulation parameters.

The tissue means come from the Gaussian mixture of the segmentation, with one exception: GM and WM use the mixing-weighted mean over all their Gaussians, while CSF uses only its darkest Gaussian. The segmentation models CSF with two Gaussians and the brighter one regularly covers GM, so their mean would be far too bright — and that value is not only the CSF intensity of the synthesis but also the low anchor of the LAS correction and the skull stripping. A cached segmentation that describes CSF by a single Gaussian is reported with a warning.

This label-driven synthesis minimizes dependence on the initial segmentation while preserving realistic tissue topology. RF fields can be predefined (MNI A/B/C) or simulated, and contrast-to-noise ratio plus voxel size are user-controlled.

![Example outputs](docs/T1-MRI-Phantom-Examples.png)

## Cortical thickness and PVE simulation
To validate cortical thickness pipelines, the label image can be edited directly:

- Cortical thickness is defined geometrically: the label is smoothed and grey-closed to repair thin WM, then GM is grown outward from the WM with CAT's exact Euclidean distance transform (`cat_bwdist`) up to the target thickness (global or 3-region using the neuromorphometrics atlas). No band is grown around the ventricles, the corpus callosum or any other non-cortical structure, since the CSF facing them is not a sulcus.
- Thickness values can vary across regions (e.g., frontal, occipital, remaining cortex) to produce known ground truth.
- The volume is internally resampled to 0.5 mm for the simulation and written back on the original or requested grid.
- Partial volume is approximated by jittering the tissue boundaries across 15 subvoxel offsets in [-0.25, 0.25] voxels and averaging hard labels (CSF=1, GM=2, WM=3), replacing the original SPM labels in synthesis.
- The original tissue fractions are blended back in where the simulation has nothing to add: inside the excluded structures, and deep inside the WM. The hard labels give a WM fraction of exactly 1 in the interior and therefore a perfectly flat WM, whereas the original fractions carry the local variation that a simulation without thickness manipulation shows. The blend stays a safe distance away from the GM/WM boundary, so the simulated thickness is untouched, and it is capped so that the ground truth cannot fall below the WM range.

![Thickness control](docs/T1-MRI-Phantom-Thickness.png)

## Geometric thickness phantom
`mri_simulate` gives a *real brain* a constant cortical thickness, so its ground truth is only as good as the segmentation it starts from. `thickness_phantom.m` answers the complementary question — what a thickness measure does when the geometry is known exactly — by building the object analytically instead:

- The WM is a sphere with regular folds, the GM is the band of constant thickness around it, and a CSF layer surrounds the GM.
- The phantom is free of artefacts by design and fully deterministic: no noise, no bias field, no random component anywhere. The error a thickness measure shows on it is its own theoretical error and nothing else.
- Outputs are a PVE label image (`_dseg.nii`, CSF=1, GM=2, WM=3) and the ideal T1w image (`_T1w.nii`) that belongs to it, each with a JSON sidecar that records the geometry.

### Why the thickness is exact
The GM is **not** made by offsetting the radius of the folded sphere. A radial offset of `t` gives a normal thickness of `t·cos(alpha)`, where `alpha` is the angle between the radius and the surface normal, so on the flanks of the folds the true thickness would come out too small. Instead the euclidean distance `D` to the WM is computed once with `cat_bwdist`, and **both cortical boundaries are level sets of that one distance map**: the WM/GM boundary is `D = wm_offset` and the GM/CSF boundary is `D = wm_offset + thickness`. Since `D` is 1-Lipschitz with a unit gradient, two of its level sets are exactly `thickness` apart everywhere, whatever the folds look like. This is the same construction `mri_simulate` uses for its cortical band, only starting from an analytic surface rather than a segmentation. The distance map is built on a supersampled grid (factor 3 by default), so the discretization of the surface stays well below the voxel size.

### Partial volume
The PVE follows the idea of `mri_simulate`: the tissue boundary is shifted across a set of sub-voxel offsets, each offset yields a hard label image, and the results are averaged. Here the offsets are applied to the distance map directly, which shifts both boundaries along their normal by a known amount. They are the midpoints of 15 equal intervals over one voxel, so for a locally flat boundary the average reproduces the exact linear PVE ramp — `mri_simulate` uses half that range because it applies the offsets to a label map it smoothed beforehand.

### Evaluation
`thickness_phantom_eval.m` checks two independent things, which is what makes the result interpretable: without the first, an error of the thickness measure cannot be told apart from an error of the phantom.

1. **Fidelity of the phantom.** The distance between the two boundaries is recovered from the label image alone — it is supersampled, `cat_bwdist` gives the distance from the WM surface, and that distance is read out on the GM/CSF isosurface at the sub-voxel crossings along the three axes (a layer of voxels would bias the readout by a good part of a sample). Tissue volumes are compared against the supersampled grid the phantom was built on. It also reports how much of the WM surface faces a *buried sulcus*, i.e. a fold narrower than twice the thickness where the CSF is squeezed out and the two GM banks touch — the situation that makes any thickness measure overestimate.
2. **Accuracy of a thickness measure**, optionally CAT's projection-based thickness `cat_vol_pbtsimple`, which takes exactly this kind of PVE label as input, over all GM voxels and on the central surface.

![Thickness phantom](docs/T1-MRI-Phantom-ThicknessPhantom.png)

### Usage
```matlab
% one phantom with 2.5mm thickness at 0.5mm, then evaluate it
[~, ~, info] = thickness_phantom;
res = thickness_phantom_eval(info);

% a series of thickness values; they share the distance map, so this is
% barely slower than a single one
[~, ~, info] = thickness_phantom(struct('thickness', 1.5:0.5:3.5));
res = thickness_phantom_eval(info, struct('fig','phantom_qc.png'));

% the folding pattern of T1Prep/internal/thickness_phantom.py
thickness_phantom(struct('fold','spherical', 'folds',[6 6], 'radius',25, ...
                         'amplitude',2.5, 'thickness',3, 'csf',1));
```

Main parameters (see `help thickness_phantom` for all of them): `dim`, `vx`, `radius`, `amplitude`, `fold` (`'cartesian'` with a `wavelength`, or `'spherical'` with `folds` — note that the spherical pattern makes its folds arbitrarily fine towards the poles, where no voxel size resolves them), `thickness`, `csf`, `supersample`, `pve_steps`, `pve_range`.

### Results
Default phantom (radius 22 mm, fold amplitude 2.5 mm, wavelength 12 mm) at 0.5 mm, thickness 1.5–3.5 mm. `geometry` is the thickness recovered from the label image, i.e. the fidelity of the phantom; `pbt` is `cat_vol_pbtsimple` with its own defaults:

known | geometry bias | sd | pbt bias | sd | buried sulci
------|---------------|-----|----------|-----|-------------
1.5 mm | −0.014 | 0.036 | −0.068 | 0.027 | 0.0%
2.0 mm | −0.017 | 0.037 | −0.262 | 0.028 | 0.3%
2.5 mm | −0.018 | 0.036 | −0.085 | 0.013 | 4.2%
3.0 mm | −0.016 | 0.039 | −0.258 | 0.029 | 6.9%
3.5 mm | −0.018 | 0.035 | −0.243 | 0.039 | 10.3%

The geometric bias of about −0.017 mm is independent of the thickness and shrinks to −0.006 mm when the evaluation is supersampled by 5 instead of 3, i.e. most of what is left is the readout and not the phantom: the label image carries the requested thickness to well within a twentieth of a voxel. Tissue volumes agree with the supersampled reference to within 0.5% for CSF and WM and 1.5% for GM, the class that consists almost entirely of boundary voxels.

Against that, `cat_vol_pbtsimple` underestimates, and how much depends on the thickness in a way that is not monotonic. Its core algorithm is not the reason — with `pbtopt.supersimple = 1`, which switches the brain-specific refinements off, the bias is a near-constant −0.27 to −0.31 mm across the whole range. The jumps come from those refinements (myelin correction, sulcus/gyrus enhancement, sharpening, blood-vessel and topology correction), which sometimes recover about 0.2 mm of it and sometimes do not.

## Movement artefacts and ringing

`simu.motion` simulates head movement the way it actually corrupts an acquisition. Every phase-encoding line of k-space is sampled at a different time, so a line acquired after the head has moved belongs to a displaced object while the reconstruction treats all of them as one. The simulation therefore transforms the image, Fourier transforms it, and assembles k-space from contiguous blocks of lines taken from differently posed copies. The mismatch between the blocks produces the ringing and ghosting along the phase-encoding direction that is typical for motion.

Note that shifting the samples *inside* k-space, which descriptions of such tools often suggest, is not the same thing: a circular shift of k-space by Δk multiplies the image by `exp(2πi·Δk·r)` and leaves the magnitude image unchanged. What motion does is the dual of that — a translation by Δr multiplies k-space by `exp(-2πi·k·Δr)`. Translations are therefore applied as an exact linear phase ramp with no interpolation, while rotations rotate k-space itself and cannot be written as a phase, so the volume is resampled before its transform.

A scalar gives the severity, and the number of events and both amplitudes grow with `severity^1.5`:

Severity | Events | Max translation | Max rotation
---------|--------|-----------------|-------------
1 (mild) | 1 | 1.0 mm | 1.0°
2 (moderate) | 3 | 2.8 mm | 2.8°
3 (severe) | 5 | 5.2 mm | 5.2°

Each event is one instructed nod: an excursion lasting a single block plus a residual offset of 30% of its amplitude that persists to the end of the scan, since the head rarely returns exactly to its former position. Pitch dominates, as it does for real nodding.

A struct overrides single values:

Field | Meaning (Default)
------|------------------
severity | Base of all the defaults below (required, or `1` if a struct is given without it)
events | Number of motion events (`round(severity^1.5)`)
translation | Maximum translation in mm (`severity^1.5`)
rotation | Maximum rotation in degrees (`severity^1.5`)
blocks | Number of k-space blocks. This is the temporal resolution of the motion and therefore also the duration of one excursion (`32`, about 5–10 s of a typical 5 min MPRAGE, i.e. the time scale of a nod)
continuous | Amplitude of a continuous drift and tremor between the events, as a fraction of `translation` (`0.2`; `0` leaves the pose piecewise constant)
pe | Phase-encoding direction, a voxel axis (`1`,`2`,`3`) or a world axis (`'x'` left-right, `'y'` anterior-posterior, `'z'` inferior-superior) mapped to the closest voxel axis (`'y'`, the usual in-plane direction of a sagittal MPRAGE)
ordering | `'linear'` fills k-space from −kmax to +kmax, so its centre is sampled in the middle of the scan; `'centric'` starts at the centre (`'linear'`)
centre | Place the first event in the block that samples the centre of k-space (`1`)

Two properties are worth knowing when using this for validation:

- **The severity depends far more on whether the centre of k-space is hit than on the amplitude**, because the centre carries most of the energy. `centre=1` therefore forces the first event there, which keeps the severity levels comparable between images instead of leaving them to the random position of the events. Set it to `0` for a purely random time course.
- **The ground truth is unaffected**, since motion does not change the anatomy. `motion` only enters the `desc` tag of the simulated image, so runs that differ only in the motion share one label file. The pose that is subtracted from all blocks is the mean weighted by their k-space energy, which is where the object appears, so the simulated image stays where its label image is. A sub-voxel displacement along the phase-encoding axis is left over and grows with the severity, from about a tenth of a voxel for mild to about half a voxel for severe motion.

The realized motion is written to the JSON sidecar under `SimulationParameters.Motion`, including `Pose`, the full pose time course with one row per k-space block in acquisition order (translations in mm, rotations in degrees, both in world coordinates), which is the ground truth of the movement itself.

### The continuous component

A pose that is piecewise constant with a handful of steps gives a handful of discontinuities in k-space, and therefore ringing along the edges alone. A real motion-corrupted image carries the dense ripple texture of many small ones, because real motion always has drift and tremor between the deliberate movements. `continuous` adds a slow drift and two slow oscillations on top of the events, which is what makes the result look like a rejected scan rather than a blurred one.

It enters as a translation only. A translation is free in k-space (a phase ramp), while a rotation needs the volume resampled once more, and the assembly is grouped so that one resampling serves all the blocks sharing a rotation. The runtime therefore follows the number of *events*, not the number of blocks: measured on a 1 mm volume, a rotation resampling costs 2.9 s, so ~11 distinct rotations are ~30 s whether or not every block has its own translation. Giving the continuous component its own rotations would have made it 32 resamplings, i.e. ~90 s at 1 mm and ~12 min at 0.5 mm.

## Ringing

`simu.ringing` is a separate artefact and can be combined with `motion`. There are **two types**, which are different artefacts and not two settings of one.

Field | Meaning (Default)
------|------------------
strength | `0`=off, `1`/`2`/`3` = mild/moderate/severe, any positive value allowed
type | `'notch'` (default) or `'gibbs'`
pe | Axis the ringing runs along: `1`,`2`,`3`, `'x'`/`'y'`/`'z'`, or `'all'` (`'y'`)
k0 | Band centre for `notch`, in units of the Nyquist frequency (`0.6`, where mriaug puts it)

### `notch` — the pronounced regular ripples

Damps or inverts a narrow band of k-space at `|k| = k0`; the gain inside the band is `1 - strength`, so the band is removed at strength 1 and inverted above it. Narrow in k-space means far-reaching in the image, so a single spatial frequency is laid over the whole image and shows as the strong regular ripple pattern that reads as ringing at first sight — **without** the blurring that truncation brings.

This is what [mriaug](https://github.com/codingfisch/mriaug)'s `ringing3d` does, and it is what to use when the *appearance* of ringing is the goal. It is not what a scanner does: no acquisition removes an isolated band of frequencies.

With `pe` set to a single axis the band is a flat slab, which gives directional stripes. With `'all'` it is a spherical shell, which gives concentric rings around edges.

Implementation note: mriaug builds its band on uncentred array-index coordinates applied to unshifted k-space, so the band lands on one side of k-space only (at about `0.6` Nyquist), and its `torch.fft.irfftn` call then re-symmetrises the result. The version here uses a properly centred, symmetric band at the same `|k|`, which is real by construction. Checked side by side on a 1 mm T1, the two are visually indistinguishable and differ by 0.063 vs 0.058 relative difference at equal depth.

### `gibbs` — the physical one

The ripple that a finite acquisition matrix really produces: keep only the central part of k-space and reconstruct on the same grid. It overshoots a step edge by about 9%, as a rect window should.

Strength | k-space kept | ripple period
---------|--------------|--------------
1 | 0.55 | ~3.6 voxels
2 | 0.40 | ~4.9 voxels
3 | 0.25 | ~8 voxels

Two properties are worth knowing, and both were got wrong in the first version of this option:

- **The overshoot is ~9% whatever the fraction is** — that is the Gibbs constant. Only the *period* of the ripples changes, roughly as `2/fraction` voxels. Above a fraction of about 0.6 the period is two to three voxels, which is the scale of the image texture, so it reads as noise and the only visible effect is the blur. Blurring and ringing cannot be separated here: they are two sides of the same truncation, and buying visible ripples always costs resolution. That is exactly why `notch` exists as a separate type.
- **Only the phase-encoding axis is truncated by default.** That is the axis whose matrix an acquisition actually shortens while the readout is oversampled, and it is why clinical Gibbs ringing appears as bands along one direction. Truncating all three axes equally is an isotropic softening and reads as smoothing, not ringing.

### Both types

The `desc` tag keeps them apart (`Ringing2` for the notch, `Ringing2Gibbs`), so runs of the two types cannot overwrite each other.

**Gibbs ringing and motion ringing are two different things.** Motion ringing comes from the mismatch between k-space lines acquired at different times; Gibbs ringing from k-space ending at a finite frequency. Only the latter is in every image whether the head moved or not. When both are requested they share one k-space: the truncation is applied to the k-space the motion assembled, before the magnitude image is reconstructed, which is the order a scanner produces them in — and that is measurably not the same as applying one after the other.

### Comparison with mriaug

[mriaug](https://github.com/codingfisch/mriaug) is a natural reference, since it also works in k-space. Its `motion3d` is

```python
offset = intensity * fft.fftn(translate3d(x, translate=translate))
return modify_k_space(x, gain=1 - intensity, offset=offset)   # k*gain + offset
```

i.e. `K = (1-α)·FFT(x) + α·FFT(shift(x))`. That gain and offset are uniform over all of k-space, so by linearity it equals `(1-α)·x + α·shift(x)` — a plain alpha blend of the image with a translated copy. Checked numerically on a 1 mm T1, the two agree to `6.7e-16`, i.e. the FFT round-trip does nothing. There is no k-space segmentation, hence no discontinuity and no ringing: of the residual energy only 11% lies in the outer half of k-space, against 48% for the block model here. Visually it is a double exposure rather than a motion artefact.

mriaug's `ringing3d` is not Gibbs ringing either — it damps a narrow band of k-space by `1 - 10*intensity`, i.e. a factor of **−4** at its documented default. That is not what a finite acquisition matrix does, but it is a very effective way to *look* like ringing, so it is available here as the `notch` type (see [Ringing](#ringing)).

So mriaug is not a better starting point for the motion model, but its ringing operator is worth having, which is why `simu.ringing` offers both types.

## Cleaning the ground truth

`simu.clean` removes from the **label** what is not the tissue it looks like,
and deliberately leaves it in the **image**.

```matlab
simu = struct('name','sub-01_T1w.nii', 'snrWM',40, 'affine',1, 'clean',1);
mri_simulate(simu, struct('percent',0));
% -> sub-01_space-MNI152_res-0p5mm_desc-snr40_T1w.nii        (vessel still bright)
%    sub-01_space-MNI152_res-0p5mm_desc-Clean_dseg.nii       (vessel labelled CSF)
%    sub-01_space-MNI152_res-0p5mm_desc-Clean_label-GM_probseg.nii
```

### Why the image is not cleaned too

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

### Blood vessels and dura

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

### Dura

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

### Periventricular CSF/WM partial volume

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

### Cost and notes

- It needs the atlas partitioning of `cat_vol_partvol`, which costs minutes.
  It is computed once and shared with `closeWMHholes` when both are used.
- The label gains a `Clean` tag in its `desc` entity. The simulated image does
  not, because it is bit for bit the same with and without the option.
- `_probseg` files are `uint8` with a scale of 1/255, like CAT12's `p1`/`p2`/`p3`.
- The GM fraction is what a nogm-style correction model needs: the voxels where
  the triangular decomposition of the label overestimates GM are exactly
  `Yp0toC(label,2) - probseg_GM`.

## Affine registration of the output

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

### The deepmriprep grid

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

### Why the affine is fitted to the deformation field

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

### Interpolation and the label

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

### Notes

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

## Limitations of the artefact models

The motion model is retrospective: it splices k-space blocks of a rigidly moved object. Real motion is continuous rather than piecewise constant (which `continuous` mitigates but does not remove), and the spin-history and coil-sensitivity effects of an inversion-recovery sequence like MPRAGE are not reproduced. The appearance and severity of the artefact are realistic; its fine structure is not a substitute for a real motion-corrupted acquisition.

## Requirements
- MATLAB with SPM12 (or SPM25) and CAT >= 26 in the path (`cat_main_LASsimple` is required and is checked for at startup)
- A T1-weighted NIfTI image, `.nii` or `.nii.gz` (default examples use `colin27_t1_tal_hires.nii`)
- `BlaiottaTPM.nii` and `BlaiottaSegmentJob.m` next to `mri_simulate.m`; the run stops with an error if either is missing

No MATLAB toolboxes beyond base MATLAB are needed. Distances use CAT's `cat_bwdist` (and `cat_vbdist` where the index of the nearest voxel is needed) rather than the Image Processing Toolbox, so results do not depend on which toolboxes are installed. The Parallel Computing Toolbox is used only when several input images are given.

## Inputs
### simu: Simulation parameters (struct)

Parameter | Description (Default)
----------|------------------------
name | Input image(s). A single T1w filename, either `.nii` or `.nii.gz`. A compressed image is uncompressed next to the original before the segmentation runs, and the image outputs are written compressed as well (the JSON sidecars are always plain text). (Default: `''`)
snrWM | add Rician magnitude noise with target SNR for white matter. Uses WM mean to derive complex noise sigma; when set, `pn` is ignored. (Default: `40`)
pn | If `>0`,add Gaussian noise as percent of the WM peak. (Default: `0`)
rng | RNG seed for reproducible noise. A fixed number gives every image the same noise pattern; `NaN` or `[]` derive the seed from the filename instead, so each image gets its own reproducible noise. (Default: `0`)
contrast | Power-law exponent for contrast change. Image is normalized to [0,1], transformed as Y.^contrast, then rescaled to original min/max. Meaningful values to simulate contrast are 0.5 (low contrast) and 1.5 (high contrast). (Default: `1`)
motion | Movement artefacts. Scalar severity (`0`=off, `1`/`2`/`3` = mild/moderate/severe, intermediate and larger values allowed), or a struct with `severity`, `events`, `translation` (mm), `rotation` (deg), `blocks`, `continuous`, `pe`, `ordering` and `centre` to override single values. See [Movement artefacts and ringing](#movement-artefacts-and-ringing). (Default: `0`)
ringing | Ringing. `0`=off, `1`/`2`/`3` = mild/moderate/severe, or a struct with `strength`, `type` (`'notch'` default, or `'gibbs'`), `pe` and `k0`. Independent of `motion` and combinable with it. See [Ringing](#ringing). (Default: `0`)
derivative | If `1`, save outputs into BIDS derivatives at the dataset root: `derivatives/mri_simulate-<version>/sub-*/ses-*/...`, mirroring the subject/session path. Thickness simulations use `mri_simulate_thickness-<version>`. (Default: `1`)
resolution | Output voxel size: scalar (applied to x,y,z) or `[x y z]`. `NaN` keeps the original resolution. Ignored when `affine` is active, since the target grid already fixes the voxel size. (Default: `NaN`)
clean | Clean the ground truth label of blood vessels, dura and the periventricular CSF/WM partial volume, while the simulated image keeps them. `0`=off, `1`=on with defaults, or a struct with `bv`, `dura`, `pve` and `probseg`. See [Cleaning the ground truth](#cleaning-the-ground-truth). (Default: `0`)
affine | Write the outputs affinely registered onto a common grid instead of the grid of the input. `0`=off, `1`=on with defaults, `'deformation'`/`'spm'` to pick the method, or a struct with `method`, `mask`, `grid`, `dim`, `mat`, `bb`, `vx`, `interp` and `space`. See [Affine registration of the output](#affine-registration-of-the-output). (Default: `0`)
WMH | Strength of white matter hyperintensities. `0`=off; `1`=mild; `2`=medium; `3`=strong; values `>=1` allowed. Larger values broaden the WMH prior via exponent `1/(WMH-0.8)` and scale the label contribution by `~1/WMH^0.75`. Constrained to (eroded) WM and modulated by a random field. (Default: `0`)
atrophy | Atrophy specification: `{atlasName, roiIds[], factors[]}`; factors >1 increase CSF (reduce GM) within ROIs. Either thickness or atrophy can be simulated. (Default: `[]`)
thickness | Cortical thickness in mm. Scalar = global; 3-vector = `[occipital rest frontal]` using neuromorphometrics atlas masks. The volume is internally resampled to 0.5 mm and written back at the requested resolution. The non-cortical structures of the atlas (subcortical grey matter, cerebellum, brainstem, hippocampus, vessels, basal forebrain) keep their original labels, and no cortical band is grown around them or around the ventricles. Either thickness or atrophy can be simulated. (Default: `0`)
closeWMHholes | Detect and fill existing WMHs in WM so that the simulated image starts from a clean WM, which allows synthetic WMHs to be added via `WMH`. Costs minutes, since it runs CAT's WMH detection. (Default: `0`)
parpool | Number of workers used when several input images are given and the Parallel Computing Toolbox is available; limited to the number of images. (Default: half the available cores)

### rf: RF bias field parameters (struct)

Parameter | Description (Default)
----------|------------------------
percent | Peak-to-peak amplitude in percent; negative values invert the field, i.e. swap its bright and dark areas at the same amplitude. (Default: `30`)
type | `'A'|'B'|'C'` (predefined MNI fields) or numeric `[strength rngSeed]` for a simulated field. Strength in `1..4` (3–4 ~ stronger 7T-like). (Default: `[2 0]`)
save | Save the simulated bias field only when `type` is numeric; ignored for `'A'|'B'|'C'`. (Default: `0`)

## Defaults
If `simu` and/or `rf` are omitted or partially specified, missing fields are filled with defaults. If `simu.name` is empty, a file selection dialog opens.

```matlab
simu = struct('name', '', 'snrWM', 40, 'pn', 0, 'contrast', 1, ...
              'motion', 0, 'ringing', 0, 'affine', 0, 'clean', 0, ...
              'resolution', NaN, 'WMH', 0, 'atrophy', [], 'thickness', 0, ...
              'rng', 0, 'derivative', 1, 'closeWMHholes', 0, ...
              'parpool', feature('numcores')/2);
rf   = struct('percent', 30, 'type', [2 0], 'save', 0);
```

## Outputs
Output names follow the [BIDS](https://bids.neuroimaging.io) filename grammar: entity-value pairs followed by a suffix, where entity labels contain only letters and digits. All option tags are therefore collected into a single camelCase `desc` label, while the output resolution uses the standard `res` entity.

Per input image:
- Simulated image: `{entities}[_res-{vx}mm]_desc-{opts}_T1w.nii`
- Ground-truth PVE label: `{entities}[_res-{vx}mm][_desc-{anatOpts}]_dseg.nii`
- If requested (`rf.save=1`, simulated fields only): `{entities}[_res-{vx}mm]_desc-{anatOpts}Biasfield_T1w.nii`
- A JSON sidecar next to the simulated image and next to the label image

When `derivative=1`, a `dataset_description.json` is written to the root of the pipeline folder, which BIDS requires for a derivative dataset to be valid.

For a `.nii.gz` input all image outputs are written as `.nii.gz`; the JSON sidecars name the compressed input in their `Sources` field. The uncompressed working copy of the input is removed again, while `<name>_seg8.mat` stays next to the input as the segmentation cache.

If the input name ends in `_T1w`, its entities are preserved and the new `res`/`desc` entities are inserted before the suffix (an existing `desc` entity of the input is replaced). For a non-BIDS input the basename is kept unchanged, so the result cannot be fully BIDS-valid — only the entities and the suffix added here follow the specification.

Notes:
- `{opts}` combines the noise tag (`snr30`, or `pn3` when `pn>0` and `snrWM=0`), the bias field (`Rf20A`, `Rf20T2`, `RfNeg20A` for an inverted field), the contrast (`ConLow`/`ConHigh`/`Con1p3`), the movement artefacts (`Motion2`), the ringing (`Ringing2`) and the anatomical tags.
- `{anatOpts}` covers only the anatomical options (`hammersRoi28F2`, `hammersMulti`, `Wmh2`, `Thickness25mm`, `Thickness15to25mm`). The label file omits the noise tag, since the ground truth does not depend on the noise level, so runs that differ only in SNR share one label file.
- Decimal points become `p` (`Con1.3` → `Con1p3`), since BIDS labels must be alphanumeric.
- `res` gives the voxel size in mm with the same `p` convention (`res-0p5mm`, `res-1mm`). Anisotropic voxels are listed per axis (`res-0p5x0p5x1p5mm`) rather than averaged — the previous mean-based `res-08mm` form was both unreadable and lossy, mapping `[0.75 0.75 0.75]` and `[0.5 0.5 1.5]` onto the same name.
- With no options at all the label would be empty, so `desc-simu` is used to keep the output distinct from the input.
- The `dseg` suffix normally implies integer labels; here the values are continuous PVE labels, which the label sidecar documents.
- When thickness is used, the label is PVE-like from the boundary jittering averaging.
- When WMH is used, a 4th label contribution is added (WMH).

## Usage
```matlab
mri_simulate(simu, rf);
```

### 1b) Rician noise at target WM SNR
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 30, ...
              'resolution', NaN, 'rng', 0);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

## Examples

### 1) Basic simulation with specific noise and 0.5 mm voxels
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 10, ...
              'resolution', 0.5, 'atrophy', [], 'rng', 42);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

### 2) Advanced simulation with atrophy (~10% in left middle frontal gyrus and ~15% in right middle frontal gyrus based on Hammers atlas), custom RF field and thicker slices
```matlab
simu = struct('name', 'custom_t1.nii', 'snrWM', 30, ...
              'resolution', [0.5, 0.5, 1.5], 'rng', []);
simu.atrophy = {'hammers', [28, 29], [2, 3]};
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 3) Thickness simulation (region-wise values, original resolution)
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 10, ...
              'resolution', NaN, 'atrophy', [], 'rng', [], ...
              'thickness', [1.5 2.0 2.5]);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

### 4) WMH simulation (medium strength) with simulated RF field
```matlab
simu = struct('name', 'custom_t1.nii', 'snrWM', 30, 'resolution', NaN, ...
              'WMH', 2, 'rng', []);
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 5) Interactive mode for example 4:
```matlab
simu = struct('snrWM', 30, 'resolution', NaN, ...
              'WMH', 2, 'rng', []);
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 6) Compressed BIDS input
```matlab
simu = struct('name', 'bids/sub-01/anat/sub-01_T1w.nii.gz', 'snrWM', 40, ...
              'resolution', NaN, 'rng', 0);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
% -> bids/derivatives/mri_simulate-<version>/sub-01/anat/sub-01_desc-snr40Rf20A_T1w.nii.gz
```

### 7) Apply contrast change (power-law)
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 30, ...
              'contrast', 1.3, 'resolution', NaN, 'rng', 0);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

### 8) Movement artefacts of moderate severity
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40, 'motion', 2);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

### 9) Movement artefacts with full control
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40);
simu.motion = struct('severity', 2, 'events', 3, 'translation', 4, ...
                     'rotation', 4, 'pe', 'y', 'ordering', 'linear');
mri_simulate(simu);
```

### 10) Ringing, alone and combined with movement
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40, 'ringing', 2);
mri_simulate(simu, struct('percent', 0));   % notch, the pronounced ripples
simu.motion = 2;
mri_simulate(simu, struct('percent', 0));

% the physical variant, and a spherical shell giving concentric rings
simu.motion  = 0;
simu.ringing = struct('strength', 2, 'type', 'gibbs');
mri_simulate(simu, struct('percent', 0));
simu.ringing = struct('strength', 2, 'pe', 'all');
mri_simulate(simu, struct('percent', 0));
```

A series that differs only in the severity shares one ground truth label file, which makes it a graded test set for the robustness of a morphometry pipeline against motion:

```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 40);
for severity = [0 1 2 3]
  simu.motion = severity;
  mri_simulate(simu, struct('percent', 0));
end
```

## File Naming Examples

For a BIDS input named `sub-01_T1w.nii`:

```
sub-01_desc-pn3Rf20A_T1w.nii                        % Gaussian noise at 3%
sub-01_desc-snr30Rf20A_T1w.nii                      % Rician noise at SNR=30 in WM
sub-01_res-1mm_desc-snr30Rf20A_T1w.nii              % Resampled to 1.0mm
sub-01_res-0p5x0p5x1p5mm_desc-snr30Rf20A_T1w.nii    % Anisotropic voxels
sub-01_res-1mm_desc-snr30RfNeg20A_T1w.nii           % Inverted bias field
sub-01_desc-snr30Rf20T2Con1p3_T1w.nii               % Simulated field, strength 2, contrast 1.3
sub-01_desc-snr30Rf20AMotion2_T1w.nii               % Moderate movement artefacts
sub-01_desc-snr30Rf20AMotion2Ringing2_T1w.nii       % Movement and notch ringing
sub-01_desc-snr30Rf20ARinging2Gibbs_T1w.nii         % The gibbs ringing type
sub-01_dseg.nii                                     % Shared by every motion/ringing run
sub-01_desc-snr30Rf20AThickness25mm_T1w.nii         % 2.5mm constant thickness
sub-01_desc-Thickness25mm_dseg.nii                  % Label for that thickness run
sub-01_desc-hammersRoi28F2_dseg.nii                 % Label for a single-ROI atrophy run
sub-01_desc-biasfield_T1w.nii                       % Saved RF field
```

Existing entities of the input are preserved and the new ones inserted before the suffix, e.g. `sub-01_ses-1_acq-mprage_T1w.nii` becomes `sub-01_ses-1_acq-mprage_desc-snr30Rf20A_T1w.nii`. A `sub-01_T1w.nii.gz` input gives the same names with a `.nii.gz` extension.

For a non-BIDS input such as `colin27_t1_tal_hires.nii` the basename is kept, so the result is not BIDS-valid even though the added entities are:

```
colin27_t1_tal_hires_desc-snr30Rf20A_T1w.nii
colin27_t1_tal_hires_dseg.nii
```
