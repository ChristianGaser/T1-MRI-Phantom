# t1sim

Simulated T1-weighted MR images with a known ground truth, for validating segmentation and cortical thickness methods and for training data.

A real scan is turned into tissue fractions of CSF, GM and WM. The simulated image is rendered from these fractions, with a bias field, noise, movement and other artefacts added, and the same fractions are its ground truth. Because the ground truth generated the image, it is exact for that image in a way the segmentation of a real scan can never be. [Where that holds and where it does not](#what-the-ground-truth-is-and-is-not) is listed below.

![From the scan to the simulated image and its ground truth](docs/t1sim-scheme.png)

## Two implementations

- **`t1sim`**, the Python package in `src/t1sim/`, which this documentation describes.
- **The MATLAB tools** in [`matlab/`](matlab/README.md), from which it was ported.

| Step | `t1sim` | MATLAB tools |
| --- | --- | --- |
| Tissue fractions from a real scan: denoising, segmentation, intensity normalisation | yes | yes |
| Cleanup of the ground truth: blood vessels, dura, periventricular partial volume | yes | yes |
| Changed anatomy: cortical thickness | yes | yes |
| Changed anatomy: atrophy | yes, by a deformation | by scaling the CSF fraction |
| Changed anatomy: white matter hyperintensities | not yet | yes |
| Synthesis of the T1w image from the fractions | yes | yes |
| Artefacts: contrast, bias field, voxel size, movement, ringing, noise | yes | yes |
| BIDS output with JSON sidecars | yes | yes |
| Output on the grid of the scan | yes, the default; template space with `grid="working"` | yes, the default |
| Geometric thickness phantom and its evaluation | yes | yes |

[Agreement with the MATLAB tools](docs/matlab-agreement.md) lists where the two give the same result and where `t1sim` differs on purpose.

## Quick start

```bash
git clone https://github.com/ChristianGaser/T1-MRI-Phantom.git
cd T1-MRI-Phantom
pip install -e ".[frontend]"
```

If pip stops at T1Prep, see [Installation](docs/installation.md): the version it asks for may not be released yet.

```bash
t1sim sub-01_T1w.nii.gz --snr 30 --bias 20                   # on the grid of the scan
t1sim sub-01_T1w.nii.gz --snr 30 --bias 20 --grid working    # in template space, 0.5 mm
```

```python
from t1sim import simulate

truth = simulate.ground_truth("sub-01_T1w.nii.gz")     # once per anatomy, about two minutes
for snr in (20, 40):
    acquired = simulate.acquire(truth, snr=snr)        # seconds each
    simulate.write(truth, acquired)
```

Every run writes the simulated image, the label map of its ground truth, the three tissue fractions and JSON sidecars, by default on the grid of the scan: the image lies on the scan voxel by voxel. They go to `derivatives/t1sim-<version>/` of a BIDS dataset, or to a `derivatives` folder next to the folder of the scan; the command prints the name of every image it writes. The first run downloads the model files (44 MB).

![What the acquisition and a change of the anatomy do](docs/t1sim-examples.png)

## What the ground truth is, and is not

- **It is the fractions the image was mixed from.** Inside the brain every voxel of the image is the sum of the three tissue intensities weighted by its fractions. Contrast, bias field and noise do not change that; a coarser voxel size averages the fractions over the new voxels.
- **It is defined by intensity, not by anatomy.** The fractions come from the intensities of the scan. The deep grey nuclei, which are brighter than cortex, are mixtures of GM and WM in it, as they are for any method that reads intensities.
- **The cleanup makes it differ from the image on purpose.** Blood vessels, dura and the rims of the ventricles look like grey matter. They are corrected in the ground truth and stay in the image, 29 ml on one scan and 23 ml on another. This is the part a method has to get right without reading it off the intensity.
- **The anatomy is made in template space.** It is made on a 0.5 mm grid that is affinely registered to MNI152NLin2009cAsym, which scales the head: on one scan volumes are 1.26 times those of the subject, and a thickness or a region is defined there. By default the image and its ground truth are then written on the grid of the scan, in the size of the subject, where the ground truth is the mean of the fractions over every voxel of the scan. `grid="working"` writes the anatomy as it was made.
- **Movement leaves a small offset.** After a simulated movement the image lies 0.1 to 0.5 mm from its ground truth, and a rotation of the apparent object is not removed.
- **Most of this was measured on one or two scans.** The numbers in this documentation say how the tool behaves there, not how it varies between scanners.

[Validating a method](docs/validation.md) says what follows from this for a comparison.

## Documentation

| Document | Contents |
| --- | --- |
| [Installation](docs/installation.md) | Requirements, installation, the model files and where they are kept |
| [Simulating images](docs/simulation.md) | The three steps, all options, the files that are written, the command line, timings |
| [Validating a method](docs/validation.md) | What to compare with what, which settings to use, and how exact the ground truth is |
| [Training data for deepmriprep](docs/training-data.md) | Grid, file pairing and the settings that matter for training |
| [The front end](docs/frontend.md) | How the tissue fractions are made from a scan, and the grid they lie on |
| [Cleaning the ground truth](docs/cleanup.md) | Vessels, dura and periventricular partial volume |
| [Changing the anatomy](docs/anatomy.md) | Cortical thickness and atrophy, and a comparison of three kinds of atrophy |
| [Artefacts of the acquisition](docs/artefacts.md) | Contrast, bias field, voxel size, movement, ringing and noise, and what was checked |
| [Geometric thickness phantom](docs/phantom.md) | A folded phantom whose cortical thickness is known by construction |
| [Agreement with the MATLAB tools](docs/matlab-agreement.md) | Where `t1sim` gives the result of `mri_simulate.m`, and where it differs on purpose |
| [Development](docs/development.md) | Layout of the repository, tests, how the figures are made |
| [MATLAB tools](matlab/README.md) | `mri_simulate.m` and the thickness phantom, with their own documentation |

## License

GNU General Public License, version 2; see [LICENSE](LICENSE).
