# Development

Part of the [documentation of t1sim](../README.md#documentation).

## Layout of the repository

| Path | Contents |
| --- | --- |
| `src/t1sim/frontend.py` | The front end: `segment`, `fractions_from_label`, `label_from_fractions`, `tissue_means` |
| `src/t1sim/models.py` | Download and check of the model files |
| `tests/` | Tests of both modules |
| `docs/` | Documentation of `t1sim`, its figures and the script that makes them |
| `models/` | The model files, once downloaded; ignored by git |
| `matlab/` | The MATLAB tools with their own documentation in `matlab/docs/` |
| `pyproject.toml` | Packaging; the version is `t1sim.__version__` |

## What t1sim takes from other packages

`t1sim` does not copy the segmentation. The front end imports three things and fails with a clear message if one is missing:

| Function | Package | Used for |
| --- | --- | --- |
| `cat_surf.vol_sanlm` | cat-surf | Denoising |
| `t1prep.unified_seg` (`segment`, `brain_mask`, `working_grid_maps`) | T1Prep | Unified segmentation on the working grid, brain mask |
| `t1prep.apply_LAS` | T1Prep | Local intensity normalisation |

Only `src/t1sim/frontend.py` imports them, inside `_t1prep()`, so the rest of the package works without them.

## Tests

```bash
python -m unittest discover tests
```

The tests need numpy and nibabel only. Neither the model files nor T1Prep are used: the download is tested with local files through a `file://` URL, and the front end with stand-ins for the denoising, the segmentation and LAS. What the real chain does to a real scan is covered by the [comparison with the MATLAB tool](validation.md), which is not automated.

## Figures

`docs/make_figures.py` makes both figures of the documentation from a real scan:

```bash
python docs/make_figures.py --scan sub-01_T1w.nii.gz --matlab DIR --cache slices.npz
```

| Option | Meaning |
| --- | --- |
| `--scan` | The T1w scan the front end runs on |
| `--matlab` | Folder with the results of `matlab/mri_simulate.m` for the same scan: `<name>_desc-Clean_dseg`, `<name>_desc-simu_T1w` and `<name>_desc-snr25Rf45T4_T1w`. They fill the second row of the scheme; without the option that row shows its steps only |
| `--cache` | Keeps the slices in this file, so that a change of the layout does not run the segmentation again |
| `--out` | Where the figures go (`docs/`) |

It needs matplotlib and scipy in addition to the front end.

| Figure | Shows |
| --- | --- |
| `docs/t1sim-scheme.png` | The path from the scan to the simulated image, and which part of it `t1sim` covers |
| `docs/t1sim-frontend.png` | What `frontend.segment` returns: an axial slice and an enlarged detail of each volume |
