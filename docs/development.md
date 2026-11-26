# Development

Part of the [documentation of t1sim](../README.md#documentation).

## Layout of the repository

| Path | Contents |
| --- | --- |
| `src/t1sim/frontend.py` | The front end: a scan becomes a label map and tissue fractions |
| `src/t1sim/atlas.py` | Atlases and other maps of the template space on the working grid |
| `src/t1sim/cleanup.py` | Vessels, dura and periventricular partial volume in the ground truth |
| `src/t1sim/anatomy.py` | Cortical thickness and atrophy |
| `src/t1sim/synthesis.py` | The image of tissue fractions |
| `src/t1sim/artefacts.py` | Bias field, voxel size, contrast, noise |
| `src/t1sim/kspace.py` | Movement artefacts and ringing |
| `src/t1sim/simulate.py` | The pipeline and the command line |
| `src/t1sim/bids.py` | Names, folders and sidecars of the outputs |
| `src/t1sim/phantom.py` | The geometric thickness phantom and its evaluation |
| `src/t1sim/models.py` | Download and check of the model files |
| `src/t1sim/_volume.py` | Smoothing, distances and morphology in mm, shared by the modules |
| `tests/` | Tests of all modules |
| `docs/` | Documentation of `t1sim`, its figures and the script that makes them |
| `models/` | The model files, once downloaded; ignored by git |
| `matlab/` | The MATLAB tools with their own documentation in `matlab/docs/` |
| `pyproject.toml` | Packaging; the version is `t1sim.__version__` |

## What t1sim takes from other packages

`t1sim` does not copy the segmentation or the vessel detection; it imports them where they are needed:

| Function | Package | Used in | For |
| --- | --- | --- | --- |
| `cat_surf.vol_sanlm` | cat-surf | `frontend` | Denoising |
| `t1prep.unified_seg` (`segment`, `brain_mask`, `working_grid_maps`) | T1Prep | `frontend` | Unified segmentation on the working grid, brain mask, deformation |
| `t1prep.apply_LAS` | T1Prep | `frontend` | Local intensity normalisation |
| `t1prep.vessels` (`vessel_weight`, `blood_vessel_prior`, `protected_regions`) | T1Prep | `cleanup` | Vessel detection |
| `t1prep.utils.TEMPLATE_PATH_T1PREP` | T1Prep | `atlas` | The folder of the atlases |
| `cat_surf.vol_dilate_geodesic` | cat-surf | `cleanup` | Growth through the tissue |
| `cat_surf.vol_thickness_pbt` | cat-surf | `phantom` | Optional: a thickness measure to evaluate |
| torch | torch | `kspace`, `anatomy` | FFT and resampling of the movement, deformations |

The imports are inside the functions that need them, so the rest of the package works without these packages.

## Tests

```bash
python -m unittest discover tests
```

Neither the model files nor a real scan are used. The download is tested with local files through a `file://` URL, the front end with stand-ins for the denoising, the segmentation and LAS, and everything else on small made-up objects: a ball of white matter in a shell of cortex for the anatomy, the cleanup and the pipeline. The tests of the dura need cat-surf and those of the deformations and of the movement need torch. What the real chain does to a real scan is in the measurements of the other pages, which are not automated.

## Figures

`docs/make_figures.py` makes the figures of the documentation from a real scan, all of them with `t1sim` itself:

```bash
python docs/make_figures.py --scan sub-01_T1w.nii.gz --cache slices.npz
```

| Option | Meaning |
| --- | --- |
| `--scan` | The T1w scan the front end runs on |
| `--cache` | Keeps the slices in this file, so that a change of the layout does not run the simulation again (about four minutes) |
| `--out` | Where the figures go (`docs/`) |

It needs matplotlib and scipy in addition to the front end.

| Figure | Shows |
| --- | --- |
| `docs/t1sim-scheme.png` | The path from the scan to the simulated image, and which part of it `t1sim` covers |
| `docs/t1sim-frontend.png` | What `frontend.segment` returns: an axial slice and an enlarged detail of each volume |
| `docs/t1sim-examples.png` | One detail of the anatomy under every artefact, and with a thinner cortex |
