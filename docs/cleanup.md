# Cleaning the ground truth

Part of the [documentation of t1sim](../README.md#documentation).

The label map is made from the intensity of the scan, so everything as bright as grey matter is grey matter in it. `t1sim.cleanup` finds three kinds of voxels for which that is wrong and corrects the ground truth. The image is not corrected: it is rendered from the fractions as they were, so a vessel stays bright in the simulated scan while the ground truth calls it CSF.

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

## What is found

**Blood vessels** become CSF. Two sources are combined:

- T1Prep's detection (`t1prep.vessels`), its port of CAT12's. It admits a voxel that is bright, has CSF around it, lies away from the white matter body and where a prior from MRA scans allows it. It is built not to touch thin gyral white matter and is careful accordingly.
- Voxels that are clearly brighter than white matter. LAS brings white matter to 3 everywhere, so a brain voxel above 3.5 is no brain tissue; such voxels are taken together with the voxels above 3.2 connected to them, outside the white matter body.

**Dura** becomes CSF. It cannot be told from cortex by its shape or by its distance to the white matter. What differs is that the cortex hangs on the white matter while the dura lies beyond the subarachnoid CSF. The white matter is therefore grown through the tissue only, and grey matter it does not reach within `dura` mm is dura. The cerebellum and the deep grey matter are spared.

**Periventricular partial volume** keeps its label and loses its grey matter. A voxel that mixes CSF and WM has the intensity of GM. Next to a ventricle, the brainstem or the corpus callosum, between pure WM and pure CSF, and thin, it is decomposed again as CSF and WM only. The same label value describes both mixtures, so only the fractions carry this correction. This follows the level 3 cleanup of CAT12's `cat_main_cleanup`.

## Measured

Volumes that were relabelled, in ml of the scan:

| Scan | Vessels | Dura | Periventricular |
| --- | --- | --- | --- |
| 0.75 mm MPRAGE with few visible vessels | 0.4 | 15.4 | 13.3 |
| 0.7 mm scan with many bright vessels | 5.6 | 3.9 | 10.4 |

- **Vessels.** On the second scan 8.4 ml of the brain are brighter than 3.3. T1Prep's detection alone took 19% of them (2.0 ml), nearly all of it real: 99% within 1 mm of such a voxel. With the bright voxels added, 52% are taken. On the first scan the voxels above 3.5 are 0.05 ml, against 6.3 ml on the second.
- **Not found:** vessels as dark as grey matter. They give none of the cues the detection uses, and CAT12's atlas partitioning, which finds part of them, is not available in Python.
- **Against `mri_simulate.m`** on the first scan: MATLAB relabelled 16.0 ml as CSF on the working grid, of which 4.8 ml are CSF or outside the brain for `t1sim` anyway. Inside MATLAB's brain mask the two agree with a Dice of 0.47. Much of what only `t1sim` relabels lies at the rim that MATLAB's skull-stripping had already removed; much of what only MATLAB relabels is as dark as grey matter (median label 1.8).

The cleanup takes 35 to 40 s.
