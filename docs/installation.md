# Installation

Part of the [documentation of t1sim](../README.md#documentation).

## Requirements

- Python 3.9 or newer.
- For the package itself: numpy, nibabel and requests, which pip installs with it.
- For the front end, in addition:
  - [cat-surf](https://pypi.org/project/cat-surf/) 1.2.0 or newer, for the unified segmentation and the SANLM filter.
  - [T1Prep](https://github.com/ChristianGaser/T1Prep) 0.9.2 or newer. This is the first version in which the segmentation on the working grid (`t1prep.unified_seg`) and LAS with given tissue intensities (`t1prep.apply_LAS(..., peaks=...)`) are public.

Without the two packages of the front end, `t1sim` can be imported and the model files can be downloaded, but `frontend.segment` stops with an error that names what is missing.

## Installing

`t1sim` is not on PyPI yet. Install it from a clone of the repository:

```bash
git clone https://github.com/ChristianGaser/T1-MRI-Phantom.git
cd T1-MRI-Phantom
pip install -e ".[frontend]"     # with the front end
pip install -e .                 # or without it
```

If pip finds no T1Prep 0.9.2, PyPI does not have that version yet. Install T1Prep from its repository first and repeat the command:

```bash
pip install "git+https://github.com/ChristianGaser/T1Prep.git"
```

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
