# T1-MRI Simulator - Quick Overview

## What is it?

`mri_simulate` is a MATLAB tool that generates realistic simulated T1-weighted MRI brain images with controllable artifacts and anatomical variations. Current behavior writes BIDS-like derivatives by default and emits JSON sidecars with simulation metadata. Default noise is Rician at a target WM SNR.

## Requirements

- MATLAB with SPM12 (or SPM25) and CAT >= 26 on the MATLAB path (`cat_main_LASsimple` is required and is checked for at startup)
- No MATLAB toolboxes beyond base MATLAB. Distances use CAT's `cat_bwdist` (and `cat_vbdist` where the index of the nearest voxel is needed) rather than the Image Processing Toolbox, so results do not depend on the installed toolboxes; `parpool` (Parallel Computing Toolbox) is used only for multiple input images
- A T1-weighted NIfTI image, `.nii` or `.nii.gz` (example: `colin27_t1_tal_hires.nii`)
- `BlaiottaTPM.nii` and `BlaiottaSegmentJob.m` next to `mri_simulate.m`, which the segmentation uses; the run stops with an error if either is missing

---

## Basic Workflow

```
Input T1w Image (.nii or .nii.gz, the latter uncompressed first)
      ↓
SPM segmentation with the bundled Blaiotta head/neck TPM
(cached as <name>_seg8.mat)
      ↓
LAS intensity normalization + SANLM denoising + APRG skull stripping
      → tissue maps (GM/WM/CSF)
      ↓
Apply Modifications (optional):
  • Atrophy (reduce GM in specific regions)
  • Thickness (create uniform cortical thickness)
  • WMH (white matter hyperintensities)
      ↓
Smooth tissue fractions with the acquisition PSF (`psf`)
      ↓
Synthesize new T1w image from modified tissues
      ↓
Add RF bias field (optional; predefined A/B/C or simulated)
      ↓
Resample to target resolution
      ↓
Apply contrast change (power-law Y^x, optional)
      ↓
Add noise (Rician at target WM SNR by default; Gaussian %WM if `snrWM<=0`)
      ↓
Outputs: Simulated image + label map (+ RF field) + JSON sidecar
```

---

## Quick Start

### Minimal Example (defaults)
```matlab
% Default: snrWM=40 (Rician), bias field type [2 0], derivatives on
simu = struct('name', 'colin27_t1_tal_hires.nii');
rf = struct('percent', 20, 'type', [2 0]);
mri_simulate(simu, rf);
```

### Gaussian noise instead of SNR
```matlab
% Disable SNR-based noise, use 3 percent Gaussian noise
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 0, 'pn', 3);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

### Common Use Cases

**1. Clean reference image (no artifacts)**
```matlab
simu = struct('name', 'input.nii', 'snrWM', 0, 'pn', 0, 'resolution', 1.0);
rf = struct('percent', 0, 'type', [2 0]);
mri_simulate(simu, rf);
```

**2. Realistic clinical scan (1.5T typical)**
```matlab
simu = struct('name', 'input.nii', 'snrWM', 25, 'resolution', 1.0);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

**3. Aging study with atrophy**
```matlab
simu = struct('name', 'input.nii', 'snrWM', 25, 'resolution', 1.0);
simu.atrophy = {'hammers', [28, 29], [2, 3]};  % ~10% and ~15% GM loss in ROIs
rf = struct('percent', 20, 'type', 'B');
mri_simulate(simu, rf);
```

**4. White matter disease**
```matlab
simu = struct('name', 'input.nii', 'snrWM', 25, 'WMH', 2);  % Medium WMH
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

**5. Cortical thickness study**
```matlab
simu = struct('name', 'input.nii', 'snrWM', 25, 'thickness', 2.5);  % 2.5mm uniform
rf = struct('percent', 15, 'type', 'A');
mri_simulate(simu, rf);
```

---

## Key Parameters

### Simulation Options (simu)

| Parameter | What it does | Typical values |
|-----------|--------------|----------------|
| `name` | Input T1w image file | `.nii` or `.nii.gz` |
| `snrWM` | Target WM SNR (Rician). If >0, overrides `pn` | 20-40 (default 40) |
| `pn` | Gaussian noise % of WM mean (used when `snrWM<=0`) | 0-9 |
| `resolution` | Output voxel size in mm | NaN (keep), 0.5, 1.0, [1 1 3] |
| `contrast` | Contrast-change exponent | 0.5, 1.5, custom |
| `WMH` | White matter lesion strength | 0 (off), >=1 |
| `atrophy` | Regional GM reduction | `{'atlas', [ROIs], [factors]}` |
| `thickness` | Cortical thickness in mm | 1.5-2.5 or `[occ, mid, front]` |
| `rng` | Random seed; a fixed number gives every image the same noise | 0 (default), or `NaN`/`[]` to seed from the filename |
| `derivative` | Save into BIDS `derivatives/mri_simulate-*` | 0/1 (default 1) |
| `closeWMHholes` | Close WMH holes in deep WM (costs minutes) | 0/1 (default 0) |
| `psf` | FWHM of the acquisition PSF in units of the output voxel, i.e. the partial volume / anti-alias filter | 0 (off) to ~1.5 (default 1) |
| `parpool` | Workers used when several input images are given | default half the cores |

### RF Bias Field Options (rf)

| Parameter | What it does | Typical values |
|-----------|--------------|----------------|
| `percent` | Field strength (%) | ±20 to ±100 |
| `type` | Field pattern | `'A'/'B'/'C'` or `[strength, seed]` |
| `save` | Save bias field file (numeric types only) | 0 (no), 1 (yes) |

**RF Field Types:**
- `'A'`, `'B'`, `'C'`: Real MNI-space bias patterns
- `[strength, seed]`: Simulated smooth field (strength 1..4); e.g., `[2, 0]`

---

## What You Get

Each simulation creates **2-3 output files** plus a JSON sidecar, written to the input folder or to `derivatives/mri_simulate-<version>/...` when `derivative=1` (default). Thickness simulations use `derivatives/mri_simulate_thickness-<version>/...`:

1. **Simulated image**: `<entities>[_res-<vx>mm]_desc-<tags>_T1w.nii[.gz]`
   - Full brain with all requested effects

2. **Ground truth labels**: `<entities>[_res-<vx>mm][_desc-<effect>]_dseg.nii[.gz]`
   - CSF=1, GM=2, WM=3 (±WMH=4), with continuous partial volume values in between
   - Useful for training/validation

3. **RF bias field** (if requested and `type` numeric): `<entities>[_res-<vx>mm]_desc-<effect>Biasfield_T1w.nii[.gz]`
   - Saved only for simulated fields when `rf.save=1`

4. **JSON sidecars**: one next to the simulated image and one next to the label image
   - Includes tool info and SimulationParameters (voxel size, pn or snrWM, RF settings, thickness, WMH, atrophy)
   - The label sidecar documents the label values and their partial volume encoding

Plus a `dataset_description.json` at the root of the pipeline folder, which BIDS requires for a valid derivative dataset.

All image outputs are written compressed when the input was a `.nii.gz`, and the JSON sidecars then name the compressed input in their `Sources` field. The uncompressed working copy of the input is removed again; `<name>_seg8.mat` stays next to the input as the segmentation cache.

---

## Main Features

### Tissue Modifications

**Atrophy**
- Reduces gray matter in specific brain regions
- Based on anatomical atlases (Hammers, etc.)
- Increases CSF to simulate volume loss
- Factor 1.5 ≈ 5% reduction, 2.0 ≈ 10%, 3.0 ≈ 15%

**Cortical Thickness**
- Creates uniform cortical ribbon from white matter
- Can vary by region (occipital, middle, frontal)
- Grows the ribbon with CAT's exact distance transform `cat_bwdist` on a 0.5 mm grid and simulates the partial volume with 15 sub-voxel boundary jitters
- Excludes the non-cortical structures of the atlas automatically (subcortical grey matter, cerebellum, brainstem, hippocampus, vessels, basal forebrain), which keep their original labels; no cortical band is grown around them, around the ventricles or around the corpus callosum
- Restores the original tissue fractions in the WM interior, so the WM keeps its natural intensity variation instead of becoming perfectly flat

**White Matter Hyperintensities (WMH)**
- Simulates age/disease-related white matter changes
- Patchy distribution using random fields
- Based on real WMH probability maps
- Strength 1=mild, 2=moderate, 3=severe (>=1 allowed)

### Image Artifacts

**RF Bias Field**
- Smooth intensity inhomogeneity (B1 field)
- Predefined realistic patterns or custom simulated
- Affects image uniformity across space
- Common in clinical scanners

**Noise & Contrast**
- Rician: target SNR in WM (`snrWM`, default 40)
- Gaussian: percentage of WM mean (`pn`) when `snrWM<=0`
- Contrast change: power-law mapping Y^x after normalizing to [0,1], then rescaled; it is applied before the noise is added
- Reproducible with fixed RNG seed (default `rng=0`)

**Resolution Control**
- Isotropic or anisotropic voxels
- Simulates different scanner protocols
- Sinc interpolation for the image, linear interpolation for the label, both preceded by the acquisition PSF as anti-alias prefilter

---

## How It Works (Simplified)

1. **Segmentation**: Uses SPM unified segmentation with the bundled Blaiotta head/neck TPM (seven tissue classes) to identify GM/WM/CSF in the input image, cached as `<name>_seg8.mat`.
2. **Modification**: Alters tissue distributions based on your parameters (atrophy, thickness, WMH).
3. **Synthesis**: Recreates the T1w image as the probability-weighted mixture of the tissue means taken from the segmentation, using the modified tissue maps as weights. GM and WM use the mixing-weighted mean over their Gaussians, CSF only its darkest Gaussian, because the second CSF Gaussian usually covers GM. Non-brain keeps the intensity of the bias-corrected input.
4. **Artifacts**: Applies bias field, optional contrast change, and noise.
5. **Outputs**: Saves the simulated image, the label map, optionally the RF field, and a JSON sidecar (in derivatives by default).

---

## Requirements

- **MATLAB** (R2017b or later recommended)
- **SPM12** or SPM25 (in MATLAB path)
- **CAT >= 26** (in MATLAB path; `cat_main_LASsimple`, `cat_main_APRG`, `cat_bwdist` and `cat_vbdist` are used)
- **Input**: High-quality T1w image, `.nii` or `.nii.gz` (e.g., 0.5mm Colin27 template provided)

---

## Tips & Tricks

### Getting Started
- Use provided `colin27_t1_tal_hires.nii` as template
- Start with default parameters, then customize
- Empty `simu.name` opens file browser (interactive mode)

### Reproducibility
- Default `rng=0` is deterministic and gives every image the same noise pattern; `rng=NaN` or `rng=[]` seed from the filename instead, so each image gets its own reproducible noise
- Set `rf.type = [2, 0]` for a reproducible simulated bias field
- JSON sidecars capture key parameters automatically

### Performance
- First run: slow (needs SPM segmentation ~5-10 min)
- Later runs: fast (reuses `*_seg8.mat` exactly as it is, whatever TPM created it — delete it to segment again with the bundled Blaiotta TPM)
- Atlas warping adds time (atrophy/thickness options); `closeWMHholes=1` costs several more minutes

### Best Practices
- **Test segmentation first**: Run on input without modifications
- **Validate outputs**: Check if effects match expectations
- **Ground truth**: Use label files for quantitative validation
- **Parameter sweep**: Try multiple noise/bias levels

---

## Common Parameter Combinations

| Study | Noise | Resolution | WMH | Atrophy | Thickness | RF Field |
|-------|-------|------------|-----|---------|-----------|----------|
| **Algorithm test** | snrWM=0, pn=0 | Original | 0 | No | No | 0% |
| **Clinical 1.5T** | snrWM=25 | 1.0mm | 0-1 | No | No | 20%, A/B/C |
| **Clinical 3T** | snrWM=30 | 0.8mm | 0-1 | No | No | 20%, A/B/C |
| **7T research** | snrWM=30 | 0.7mm | 0 | No | No | 40%, [4,0] |
| **Aging study** | snrWM=25 | 1.0mm | 1-3 | Yes | No | 20%, A/B/C |

---

## Limitations

- **Single modality**: Only T1w (no T2, FLAIR, etc.)
- **No motion**: Motion artifacts not yet implemented
- **Simplified WMH**: Single intensity class
- **Requires segmentation**: Input must be segmentable by SPM

---

## Output Interpretation

### File Naming
```
<entities>[_res-<vx>mm]_desc-<tags>_T1w.nii[.gz]
<entities>[_res-<vx>mm][_desc-<effect>]_dseg.nii[.gz]
<entities>[_res-<vx>mm]_desc-<effect>Biasfield_T1w.nii[.gz]

BIDS entity labels must be alphanumeric, so the tags for noise/SNR, RF,
contrast, WMH and atrophy/thickness are concatenated into one camelCase
desc label, while the resolution uses the standard res entity. Decimal
points become 'p' and anisotropic voxels are listed per axis, e.g.:
sub-01_res-1mm_desc-snr30Rf20AWmh2_T1w.nii
sub-01_desc-Rf15BConHighThickness15mm_T1w.nii
sub-01_res-0p5x0p5x1p5mm_desc-hammersRoi28F2_dseg.nii
```

### Label Values
- **1.0** = Pure CSF
- **2.0** = Pure GM
- **3.0** = Pure WM
- **4.0** = WMH (if enabled)
- **Intermediate** = Partial volume (e.g., 2.5 = 50% GM + 50% WM)

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Run starts with a long segmentation | Normal on first run - no `*_seg8.mat` exists yet, SPM creates it |
| "the segmentation describes CSF by a single Gaussian" | The cached `*_seg8.mat` comes from another TPM; delete it so that `BlaiottaSegmentJob.m` segments again with two CSF Gaussians |
| "Segmentation job ... was not found" / "TPM ... was not found" | `BlaiottaSegmentJob.m` and `BlaiottaTPM.nii` have to sit next to `mri_simulate.m` |
| "Atlas not found" | Check atlas name spelling |
| "Cannot combine atrophy and thickness" | Choose only one |
| Simulation too slow | Atlas warping takes time; use smaller ROIs or skip |
| Unrealistic results | Check input image quality and parameters |

---

## Interactive Mode

Don't know which file to use? Just run:
```matlab
simu = struct('pn', 3, 'snrWM', 0);  % Gaussian noise only
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```
A file browser will open, and you can select any T1w image, compressed or not.

---

## Further Reading

- **Detailed documentation**: See `README.md` for the full parameter and output description
- **Code comments**: `mri_simulate.m` has extensive inline documentation
- **SPM12 manual**: https://www.fil.ion.ucl.ac.uk/spm/doc/
- **CAT12 manual**: http://www.neuro.uni-jena.de/cat12-help/

---

## Quick Reference Card

### Minimal Working Examples

```matlab
% Default everything (snrWM=40, derivatives on)
mri_simulate();

% Just Gaussian noise (disable snr). An omitted rf is not an absent one:
% the defaults still add a 30% simulated bias field, thus pass rf.percent = 0
simu = struct('name', 'input.nii', 'snrWM', 0, 'pn', 5);
mri_simulate(simu, struct('percent', 0));

% Just bias field
simu = struct('name', 'input.nii', 'snrWM', 0, 'pn', 0);
rf = struct('percent', 30, 'type', 'B');
mri_simulate(simu, rf);

% Everything combined
simu = struct('name', 'input.nii', 'snrWM', 25, 'resolution', 1.0, 'WMH', 2);
simu.atrophy = {'hammers', [28], [2]};
rf = struct('percent', 20, 'type', [3, 42], 'save', 1);
mri_simulate(simu, rf);
```

---

**Version**: 0.10.1  
**Author**: Christian Gaser  
**Repository**: github.com/ChristianGaser/T1-MRI-Phantom  
**License**: See LICENSE file
