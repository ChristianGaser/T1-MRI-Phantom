# Installation

Part of the [documentation of t1sim](../README.md#documentation).

## Requirements

- Python 3.9 or newer.
- For the package itself: numpy, scipy, nibabel and requests, which pip installs with it. This is enough for the geometric phantom, the image-domain artefacts and the model download.
- For everything that starts from a real scan (the extra `frontend`):
  - [cat-surf](https://pypi.org/project/cat-surf/) 1.2.0 or newer, for the unified segmentation, the SANLM filter and the growth through the tissue.
  - [T1Prep](https://github.com/ChristianGaser/T1Prep) 0.9.3 or newer. It provides the segmentation on the working grid with its deformation (`t1prep.unified_seg.segment(..., deformation=True)`), LAS with given tissue intensities (`t1prep.apply_LAS(..., peaks=...)`), the atlases and the vessel detection.
  - torch, which T1Prep needs as well; it runs the movement artefacts and the deformations.

Without them `t1sim` can be imported, but `frontend.segment` stops with an error that names what is missing. `t1sim` checks T1Prep by what it can do (the two functions above), not by its version number.

## Installing

`t1sim` is not on PyPI yet. Install it from a clone of the repository:

```bash
git clone https://github.com/ChristianGaser/T1-MRI-Phantom.git
cd T1-MRI-Phantom
pip install -e ".[frontend]"     # with the front end
pip install -e .                 # or without it
```

If pip stops because it finds no T1Prep 0.9.3, that version is not released yet. Install the three packages of the front end by hand, T1Prep from its repository, and `t1sim` without the extra:

```bash
pip install "cat-surf>=1.2.0" torch "git+https://github.com/ChristianGaser/T1Prep.git"
pip install -e .
```

This works with any T1Prep that has the two functions, whatever version it reports.

## Model files

The unified segmentation takes its tissue priors from five shooting templates with six classes (GM, WM, CSF, bone, soft tissue, air), 44 MB in all. They are not part of T1Prep's PyPI package, so `t1sim` downloads them itself, from the commit of the T1Prep repository that added them, and checks every file against its SHA-256. A file is kept only if its hash is right.

`frontend.segment` downloads what is missing on its first run. By hand:

```bash
python -m t1sim.models            # download what is missing
python -m t1sim.models --verify   # hash the files that are there, download those that differ
python -m t1sim.models --force    # download all of them again
python -m t1sim.models --where    # print the folder and exit
python -m t1sim.models --dir DIR  # use another folder
```

After the installation the same command is available as `t1sim-download-models`.

### Where they are kept

| Case | Folder |
| --- | --- |
| `T1SIM_MODEL_DIR` is set | that folder |
| `t1sim` runs from a clone of the repository | `models/` in the clone, which git ignores |
| `t1sim` is installed as a package | the user cache: `~/Library/Caches/t1sim/models` on macOS, `$XDG_CACHE_HOME/t1sim/models` or `~/.cache/t1sim/models` on Linux, `%LOCALAPPDATA%\t1sim\Cache\models` on Windows |

`T1SIM_MODEL_URL` names another folder to download from; a `file://` URL works, for machines without internet access.
