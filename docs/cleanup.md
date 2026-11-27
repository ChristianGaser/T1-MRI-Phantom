# Cleaning the ground truth

Part of the [documentation of t1sim](../README.md#documentation).

The fractions are made from the intensity of the scan, so everything as bright as grey matter is grey matter in them. `t1sim.cleanup` finds three kinds of voxels for which that is wrong and corrects the ground truth. The image is not corrected: it is rendered from the fractions as they were, so a vessel stays bright in the simulated scan while the ground truth calls it CSF. That is what a segmentation has to get right, and what a validation should test.

- [Options](#options)
- [What is found](#what-is-found)
- [Measured](#measured)
- [Limits](#limits)

## Options

```python
truth = simulate.ground_truth(scan)                                   # cleaned by default
truth = simulate.ground_truth(scan, clean={"dura": 5.0})              # other settings
truth = simulate.ground_truth(scan, clean=False)                      # as segmented
```

| Option | Meaning (default) |
| --- | --- |
| `vessels` | Strength of the vessel detection in 0..1; 0 switches it off (`0.5`) |
| `dura` | Length in mm of the path through the tissue beyond which grey matter is dura; 0 switches it off, lower removes more (`4.0`) |
| `pve` | Correct the partial volume of CSF and WM around the ventricles (`True`) |

`truth["fractions"]` and `truth["label"]` are the cleaned ground truth, and `truth["image_fractions"]` is what the image is rendered from. The sidecar lists the volume of each of the three corrections (`LabelCleanup`, in ml of the working grid).

After a simulated [cortical thickness](anatomy.md#cortical-thickness) the cleanup leaves the rebuilt cortex alone: it holds no vessel and no dura, and its label no longer follows the intensities the detection reads. Only the structures that keep their labels are cleaned.

## What is found

**Blood vessels** become CSF. Two sources are combined:

- T1Prep's detection (`t1prep.vessels`), its port of CAT12's. It admits a voxel that is bright, has CSF around it, lies away from the white matter body and where a prior from MRA scans allows it. It is built not to touch thin gyral white matter and is careful accordingly.
- Voxels that are clearly brighter than white matter. LAS brings white matter to 3 everywhere, so a brain voxel above 3.5 is no brain tissue; such voxels are taken together with the voxels above 3.2 connected to them, outside the white matter body. The brightest of them, from 4 on, have no fractions at all before the cleanup (see [the front end](frontend.md#how-the-fractions-are-made)) and would be holes in the ground truth.

**Dura** becomes CSF. It cannot be told from cortex by its shape or by its distance to the white matter. What differs is that the cortex hangs on the white matter while the dura lies beyond the subarachnoid CSF. The white matter is therefore grown through the tissue only, and grey matter it does not reach within `dura` mm is dura. The cerebellum and the deep grey matter are spared.

**Periventricular partial volume** loses its grey matter. A voxel that mixes CSF and WM has the intensity of GM. Next to a ventricle, the brainstem or the corpus callosum, between pure WM and pure CSF, and thin, it is decomposed again as the mixture of CSF and WM that has its intensity in the rendered image. GM is not half way between CSF and WM, so that is not the mixture with its label: with tissue intensities of 148, 380 and 501, a voxel as bright as GM is 34% CSF and 66% WM, and its label becomes 2.3. The detection follows the level 3 cleanup of CAT12's `cat_main_cleanup`.

## Measured

Volumes that were corrected, converted to ml of the scan (the sidecar gives ml of the working grid, 1.26 and 1.42 times as much for these two scans):

| Scan | Vessels | Dura | Periventricular |
| --- | --- | --- | --- |
| 0.75 mm MPRAGE with few visible vessels | 0.4 | 15.4 | 13.3 |
| 0.7 mm scan with many bright vessels | 8.9 | 3.9 | 10.4 |

- **Vessels.** On the second scan 8.4 ml of the brain are brighter than 3.3. T1Prep's detection alone took 19% of them (2.0 ml), nearly all of it real: 99% within 1 mm of such a voxel. With the bright voxels added, 92% are taken, among them all 3.3 ml that had no fractions. On the first scan the voxels above 3.5 are 0.05 ml, against 6.3 ml on the second.
- **Against `mri_simulate.m`** on the first scan, on the working grid. MATLAB relabelled 16.0 ml as CSF, of which 4.8 ml are CSF or outside the brain for `t1sim` anyway. (The 19.7 ml in [Agreement with the MATLAB tools](matlab-agreement.md#the-ground-truth-of-one-scan) are another count: every voxel whose label its cleanup changed by more than 0.1, on the grid of the scan.) Inside MATLAB's brain mask the voxels the two relabel as CSF agree with a Dice of 0.47. Much of what only `t1sim` relabels lies at the rim that MATLAB's skull-stripping had already removed; much of what only MATLAB relabels is as dark as grey matter (median label 1.8).

The cleanup takes 35 to 40 s.

## Limits

- **Vessels as dark as grey matter are not found.** They give none of the cues the detection uses, and CAT12's atlas partitioning, which finds part of them, is not available in Python.
- **The dura rule is a rule of distance.** Cortex that lies further than `dura` mm from white matter along a path through the tissue becomes CSF as well. On the two scans this did not happen visibly, but it was not measured on scans whose thin white matter blades are lost to a coarse resolution.
- **Two scans.** The numbers above say what the rules do there. How much of what they take is really vessel or dura was checked by eye and by the statistics given, not against a manual labelling.
