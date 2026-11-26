# t1sim

Simulated T1-weighted MR images with an exact ground truth, for validating segmentation and cortical thickness methods and for training data.

A real scan is turned into tissue fractions of CSF, GM and WM. These fractions are the ground truth, and the simulated image is rendered from them, with a bias field, noise, movement and other artefacts added. Because the ground truth generated the image, it is exact for that image, which the segmentation of a real scan can never be.

![From the scan to the simulated image](docs/t1sim-scheme.png)

## Two implementations

- **`t1sim`**, the Python package in `src/t1sim/`. It is in development and covers the first row of the figure: the ground truth from a real scan.
- **The MATLAB tools** in [`matlab/`](matlab/README.md). They cover both rows, and a geometric phantom with an exact cortical thickness. The images of the second row were made by `matlab/mri_simulate.m` from the same scan.

| Step | `t1sim` | MATLAB tools |
| --- | --- | --- |
| Ground truth from a real scan: denoising, segmentation, intensity normalisation | yes | yes |
| Cleanup of the label: blood vessels, dura, periventricular partial volume | not yet | yes |
| Changed anatomy: cortical thickness, atrophy, white matter hyperintensities | not yet | yes |
| Synthesis of the T1w image from the fractions | not yet | yes |
| Artefacts: bias field, resolution, contrast, movement, ringing, noise | not yet | yes |
| BIDS output with JSON sidecars | not yet | yes |
| Geometric thickness phantom and its evaluation | not yet | yes |

## Quick start

```bash
git clone https://github.com/ChristianGaser/T1-MRI-Phantom.git
cd T1-MRI-Phantom
pip install -e ".[frontend]"
```

```python
from t1sim import frontend

seg = frontend.segment("sub-01_T1w.nii.gz")

seg["label"]       # 0 outside the brain, CSF = 1, GM = 2, WM = 3, partial volumes in between
seg["fractions"]   # (X, Y, Z, 3): fractions of CSF, GM and WM
seg["corrected"]   # the denoised and bias-corrected scan, whole head
seg["means"]       # {"CSF": ..., "GM": ..., "WM": ...}, on the scale of seg["corrected"]
```

The first run downloads the model files (44 MB). A run takes about 100 s on an M2 MacBook Air. If the installation stops at T1Prep, see [Installation](docs/installation.md).

![What the front end returns](docs/t1sim-frontend.png)

## Documentation

| Document | Contents |
| --- | --- |
| [Installation](docs/installation.md) | Requirements, installation, the model files and where they are kept |
| [The front end](docs/frontend.md) | How the ground truth is made from a scan, and everything `frontend.segment` returns |
| [Agreement with the MATLAB tool](docs/validation.md) | The label of `t1sim` against that of `mri_simulate.m` on the same scan |
| [Development](docs/development.md) | Tests, how the figures are made, layout of the repository |
| [MATLAB tools](matlab/README.md) | `mri_simulate.m` and the thickness phantom, with their own documentation |

## License

GNU General Public License, version 2; see [LICENSE](LICENSE).
