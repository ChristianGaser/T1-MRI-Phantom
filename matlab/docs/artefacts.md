# Movement artefacts and ringing

Part of the [documentation of the MATLAB tools](../README.md).

`simu.motion` and `simu.ringing` add the two artefacts that are made in k-space. Neither changes the ground truth.

## Movement artefacts

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

## Limitations of the artefact models

The motion model is retrospective: it splices k-space blocks of a rigidly moved object. Real motion is continuous rather than piecewise constant (which `continuous` mitigates but does not remove), and the spin-history and coil-sensitivity effects of an inversion-recovery sequence like MPRAGE are not reproduced. The appearance and severity of the artefact are realistic; its fine structure is not a substitute for a real motion-corrupted acquisition.
