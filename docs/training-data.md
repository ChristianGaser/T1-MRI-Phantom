# Training data for deepmriprep

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim` writes image and ground truth in the form the training scripts of [deepmriprep-train](https://github.com/wwu-mmll/deepmriprep-train) read: affinely registered to the template at 0.5 mm, on their 336 x 384 x 336 grid, with the file names their pairing expects.

## A set of training pairs

```bash
t1sim sub-*_T1w.nii.gz --grid deepmriprep \
      --snr 15 30 60 --bias 30 --bias-type 2 scan \
      --resolution 1 --upsample
```

or, with full control:

```python
from t1sim import simulate

for scan in scans:
    truth = simulate.ground_truth(scan)                      # once per subject
    for snr in (15, 30, 60):
        acquired = simulate.acquire(truth, snr=snr, bias_type=(2, None),
                                    resolution=1.0, upsample=True, grid="deepmriprep")
        simulate.write(truth, acquired, out_dir="data/simu")
```

Then point `simu_path` of `2b_prep_simulate.py` at the folder.

Keep the training data in a folder of its own. A validation set of the same scans on the full working grid has the same file names, and `simulate.write` refuses to put a ground truth of one grid over that of the other.

## What matters for training

| Point | Setting | Why |
| --- | --- | --- |
| Grid | `grid="deepmriprep"` | The 336 x 384 x 336 part of the working grid that `2_prep_segment.py` cuts out; `2b_prep_simulate.py` refuses another shape |
| Resolution | `resolution=1, upsample=True` | Most real data on that grid is a 1 mm scan that was resampled. The image is acquired at 1 mm, with noise that is white there, and brought back to 0.5 mm, where the noise is correlated (0.69 between neighbours) as in such a scan. The ground truth keeps its 0.5 mm |
| Noise | `seed=None` (the default) | Every output file gets noise of its own, derived from the folder and the name of the scan and from the settings. A fixed seed gives every subject the same noise pattern, which a network can learn |
| Bias field | `bias_type=(strength, None)` | Every subject gets a field of its own. The default `(2, 0)` is one field for all |
| Contrast | `contrast=0.7` to `1.4`, for example | Another ratio between the tissue intensities, with a ground truth that stays exact: the image is mixed again from the changed intensities |
| Ground truth | `clean=True` (the default) | Vessels, dura and periventricular partial volume are no grey matter in the ground truth while the image still shows them |
| Variety of anatomy | `atrophy=...` | A thinner cortex or a smaller region gives a new anatomy with an exact ground truth; see [Changing the anatomy](anatomy.md) |

## How the files are paired

`pair_simulations` of deepmriprep-train splits a name into what precedes `desc`, the `desc` label and the suffix. An image belongs to the ground truth with the same first part whose `desc`, without `Clean`, ends the `desc` of the image. `t1sim` writes exactly that.

`2b_prep_simulate.py` also needs the grey matter fraction, `..._label-GM_probseg.nii.gz` next to the label map. From it and the label map it derives its `nogm` target: the voxels in which the label map, read as a mixture of the two tissues next to its value, holds more grey matter than the fraction says. Those are the mixtures of CSF and WM around the ventricles, whose label lies between 1 and 3 while they hold no grey matter (see [Cleaning the ground truth](cleanup.md)).

Checked on simulations of one scan with the functions of deepmriprep-train itself: every image was paired with its ground truth and its grey matter fraction, and the shape was 336 x 384 x 336.

## Keep out of the training set

A scan that is used to validate must not be trained on: a network that has seen an anatomy with its ground truth is no longer tested by it. T1Prep's phantom test is built on the scan HR075, so do not simulate training data from that scan.
