# T1-MRI Simulator - Quick Overview

## What is it?

`mri_simulate` is a MATLAB tool that generates realistic simulated T1-weighted MRI brain images with controllable artifacts and anatomical variations. Current behavior writes BIDS-like derivatives by default and emits JSON sidecars with simulation metadata. Default noise is Rician at a target WM SNR.

## Requirements

- MATLAB with SPM12 (or SPM25) and CAT12 >= 12.10 on the MATLAB path
- Image Processing Toolbox (only required when using cortical thickness simulation; uses `bwdist`)
- A T1-weighted NIfTI image (example: `colin27_t1_tal_hires.nii`)

---

## Basic Workflow

```
Input T1w Image
      ↓
SPM12 Segmentation → Extract tissue maps (GM/WM/CSF)
      ↓
Apply Modifications (optional):
  • Atrophy (reduce GM in specific regions)
  • Thickness (create uniform cortical thickness)
  • WMH (white matter hyperintensities)
      ↓
Synthesize new T1w image from modified tissues
      ↓
Add RF bias field (optional; predefined A/B/C or simulated)
      ↓
Resample to target resolution
      ↓
Add noise (Rician at target WM SNR by default; Gaussian %WM if `snrWM<=0`)
      ↓
Apply contrast change (power-law Y^x, optional)
      ↓
Outputs: Simulated image(s) + masked image + label map + JSON sidecars
```

---

## Quick Start

### Minimal Example (defaults)
```matlab
% Default: snrWM=30 (Rician), bias field type [2 0], derivatives on
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
| `name` | Input T1w image file | Any NIfTI file |
| `snrWM` | Target WM SNR (Rician). If >0, overrides `pn` | 20-40 (default 30) |
| `pn` | Gaussian noise % of WM mean (used when `snrWM<=0`) | 0-9 |
| `resolution` | Output voxel size in mm | NaN (keep), 0.5, 1.0, [1 1 3] |
| `contrast` | Contrast-change exponent | 0.5, 1.5, custom |
| `WMH` | White matter lesion strength | 0 (off), >=1 |
| `atrophy` | Regional GM reduction | `{'atlas', [ROIs], [factors]}` |
| `thickness` | Cortical thickness in mm | 1.5-2.5 or `[occ, mid, front]` |
| `rng` | Random seed; 0 for deterministic | 0 (default) or `[]` (MATLAB rng) |
| `derivative` | Save into BIDS `derivatives/mri_simulate-*` | 0/1 (default 1) |
| `closeWMHholes` | Close WMH holes in deep WM | 0/1 (default 1) |

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

Each simulation creates **3-4 output files**, written to the input folder or to `derivatives/mri_simulate-0.9.4/...` when `derivative=1` (default):

1. **Simulated image**: `<name>_desc-<tags>T1w.nii[.gz]`
   - Full brain with all requested effects

2. **Masked image**: `<name>_desc-<tags>masked_T1w.nii[.gz]`
   - Brain-only (skull stripped)

3. **Ground truth labels**: `<name>_desc-<effect>_label-seg.nii[.gz]`
   - CSF=1, GM=2, WM=3 (±WMH=4)
   - Useful for training/validation

4. **RF bias field** (if requested and `type` numeric): `<name>_desc-<effect>_RFfield.nii[.gz]`
   - Saved only for simulated fields when `rf.save=1`
 
5. **JSON sidecars**: one per simulated and masked image
   - Includes tool info and SimulationParameters (voxel size, pn or snrWM, RF settings, thickness, WMH, atrophy)

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
- Uses distance-based growth with PVE smoothing
- Excludes subcortical structures automatically

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
- Rician: target SNR in WM (`snrWM`, default 30)
- Gaussian: percentage of WM mean (`pn`) when `snrWM<=0`
- Contrast change: power-law mapping Y^x after normalizing to [0,1], then rescaled
- Reproducible with fixed RNG seed (default `rng=0`)

**Resolution Control**
- Isotropic or anisotropic voxels
- Simulates different scanner protocols
- High-quality sinc interpolation

---

## How It Works (Simplified)

1. **Segmentation**: Uses SPM12 to identify GM/WM/CSF in input image.
2. **Modification**: Alters tissue distributions based on your parameters (atrophy, thickness, WMH).
3. **Synthesis**: Recreates T1w image from modified tissue maps via the SPM mixture model.
4. **Artifacts**: Applies bias field, optional contrast change, and noise.
5. **Outputs**: Saves simulated image(s), masked version, label map, JSON sidecars (in derivatives by default).

---

## Requirements

- **MATLAB** (R2017b or later recommended)
- **SPM12** (in MATLAB path)
- **CAT12** (in MATLAB path)
- **Input**: High-quality T1w image (e.g., 0.5mm Colin27 template provided)

---

## Tips & Tricks

### Getting Started
- Use provided `colin27_t1_tal_hires.nii` as template
- Start with default parameters, then customize
- Empty `simu.name` opens file browser (interactive mode)

### Reproducibility
- Default `rng=0` is deterministic; set `rng=[]` to use MATLAB rng; `rng=NaN` seeds from filename
- Set `rf.type = [2, 0]` for a reproducible simulated bias field
- JSON sidecars capture key parameters automatically

### Performance
- First run: slow (needs SPM segmentation ~5-10 min)
- Later runs: fast (reuses `*_seg8.mat`)
- Atlas warping adds time (atrophy/thickness options)

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
<name>_desc-<tags>T1w.nii[.gz]
<name>_desc-<tags>masked_T1w.nii[.gz]
<name>_desc-<effect>_label-seg.nii[.gz]
<name>_desc-<effect>_RFfield.nii[.gz]

Tags combine noise/SNR, RF, WMH, atrophy/thickness, and resolution, e.g.:
colin27_desc-snr30_rf20_A_res-10mm_WMH2_T1w.nii
colin27_desc-conHigh_rf15_B_thickness15mm_T1w.nii
colin27_desc-rf20_2_res-10mm_hammers_28_2_label-seg.nii
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
| "seg8.mat not found" | Normal on first run - SPM will create it |
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
A file browser will open, and you can select any T1w image.

---

## Further Reading

- **Detailed documentation**: See `T1_SIMULATION_SCHEME.md` for technical details
- **Code comments**: `mri_simulate.m` has extensive inline documentation
- **SPM12 manual**: https://www.fil.ion.ucl.ac.uk/spm/doc/
- **CAT12 manual**: http://www.neuro.uni-jena.de/cat12-help/

---

## Quick Reference Card

### Minimal Working Examples

```matlab
% Default everything (snrWM=30, derivatives on)
mri_simulate();

% Just Gaussian noise (disable snr)
simu = struct('name', 'input.nii', 'snrWM', 0, 'pn', 5);
mri_simulate(simu);

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

**Version**: 0.9.4  
**Author**: Christian Gaser  
**Repository**: github.com/ChristianGaser/t1-mri-simulator  
**License**: See LICENSE file
