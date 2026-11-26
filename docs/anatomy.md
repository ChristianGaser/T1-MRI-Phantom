# Changing the anatomy

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim.anatomy` changes the ground truth before the image is rendered from it, so the simulated image shows an anatomy that is known: a cortex of a chosen thickness, or a region with atrophy.

- [Cortical thickness](#cortical-thickness)
- [Atrophy](#atrophy)
- [Which kind of atrophy](#which-kind-of-atrophy)

## Cortical thickness

```python
truth = simulate.ground_truth(scan, thickness=2.5)               # 2.5 mm everywhere
truth = simulate.ground_truth(scan, thickness=(1.5, 2.0, 2.5))   # occipital, rest, frontal
```

The cortex is rebuilt as a band of the chosen thickness around the white matter. This is a port of `simulate_thickness` of `matlab/mri_simulate.m`:

1. The label map is smoothed with 1.25 mm and thin white matter is repaired by a grey-scale closing.
2. The outer boundary of the cortex is a level set of the Euclidean distance to the white matter surface. The distance is seeded with the depth of every white matter voxel below that surface, so its zero lies on the surface and not on the voxel centres.
3. Partial volumes come from repeating this with the boundaries shifted by 15 sub-voxel offsets and averaging the hard labels.
4. Non-cortical structures of the Neuromorphometrics atlas keep their labels, and no band is grown around the ventricles. Deep inside the white matter the original fractions are restored, so rendered white matter keeps its texture.

Measured:

| Object | Asked for | Result |
| --- | --- | --- |
| Sphere at 0.5 mm, thickness along rays through the centre | 1.5 / 2.5 / 3.5 mm | 1.50 / 2.49 / 3.49 mm |
| The same, from the volumes of GM and WM | 1.5 / 2.5 / 3.5 mm | 1.46 / 2.47 / 3.47 mm |
| Real brain (0.75 mm MPRAGE), projection-based thickness of cat-surf, median and 10 to 90% | 2.5 mm | 2.66 mm, 2.55 to 2.71 |
| The same | 1.5 mm | 1.66 mm, 1.62 to 1.73 |
| The same, before the simulation | - | 2.89 mm, 2.27 to 3.48 |

The projection-based thickness reads 0.15 to 0.18 mm more than was asked for on the sphere as well (1.68 / 2.65 / 3.65 mm), so the offset on the real brain is that of the measure, not of the simulation: by default it adds 0.25 mm (`correct_thickness`), a correction for the label maps of cat-surf's own pipeline. A run takes 7 to 8 minutes on the 0.5 mm grid, most of it in 75 distance transforms.

For a thickness that is exact by construction, without a real brain, see the [geometric phantom](phantom.md).

## Atrophy

```python
truth = simulate.ground_truth(scan, atrophy={"atlas": "Hammers", "regions": [28], "change": -0.10})
truth["jacobian"]                               # the ground truth of the volume change
truth["parameters"]["Atrophy"][0]["Achieved"]   # what was reached
```

The tissue of a region changes its volume by `change`, by a smooth deformation that is applied to the fractions and to the image. Boundaries move, and the Jacobian determinant of the deformation is the ground truth of the volume change. Two kinds:

- **`mode="thinning"`** (default). Every point of the cortex of the region moves towards the white matter surface in proportion to its distance from it. The cortex keeps its area and gets thinner; the CSF beyond it widens. The white matter surface stays.
- **`mode="volume"`**. The tissue of the region is the sink and the CSF within `reach` the source of a curl-free flow. The tissue shrinks the same in all directions.

In both the strength is adjusted until the tissue of the region has changed by `change`. Several regions are deformed one after the other. A positive `change` thickens or enlarges.

## Which kind of atrophy

`matlab/mri_simulate.m` multiplies the CSF fraction inside the region by a factor and renormalises. The three were compared on the left middle frontal gyrus (Hammers 28) of a 0.75 mm MPRAGE, on the 0.5 mm grid. The thickness is the median projection-based thickness of cat-surf in the region, 2.62 mm before the change (2.87 mm as it reports it, less the 0.25 mm it adds by default).

| Method | GM in the region | Thickness | Pial boundary | WM boundary | Intensity of the cortex |
| --- | --- | --- | --- | --- | --- |
| Before | 26.7 ml | 2.62 mm | - | - | - |
| Thinning, -10% | -9.7% | -0.26 mm (-10.0%) | -0.34 mm | -0.05 mm | -1.5% |
| Volume, -10% | -6.1% | -0.09 mm (-3.4%) | -0.18 mm | -0.01 mm | -0.5% |
| Scaling of CSF, factor 2 | -9.5% | -0.32 mm (-12.1%) | -0.29 mm | 0 | -2.1% |
| Scaling of CSF, factor 3 | -15.8% | -0.53 mm (-20.4%) | -0.48 mm | 0 | -3.3% |

The same anatomy averaged to 1 mm voxels before the change:

| Method | GM in the region at 0.5 mm | at 1 mm |
| --- | --- | --- |
| Thinning, -10% | -9.7% | -9.6% |
| Volume, -10% | -6.1% | -6.0% |
| Scaling of CSF, factor 2 | -9.5% | -10.6% |
| Scaling of CSF, factor 3 | -15.8% | -17.5% |

What this shows:

- **Scaling does thin the cortex**, because nearly every cortical voxel of a label map made from intensities holds some CSF. But its effect is set by a factor without a unit, it depends on the voxel size, and it darkens the cortex, since the voxels lose grey matter where they are instead of the boundary moving. It gives no deformation to validate against. It is kept as `anatomy.atrophy_by_scaling` for comparison only.
- **A change of volume** is the right model when the volume is what is tested, and its Jacobian is the ground truth for deformation-based morphometry. It thins the cortex by only a third of the volume change, and in a fixed region less grey matter is lost than the tissue of the region lost, because neighbouring tissue moves in.
- **Thinning** is the model of cortical atrophy: the loss of volume is a loss of thickness, the same at any voxel size, with the intensity of the cortex nearly kept.

Limits of the thinning: where two banks of a sulcus touch without CSF between them, they cannot move apart, and the cortex thins less there; the thin CSF of a narrow sulcus is stretched instead, so the Jacobian determinant is large in single voxels (up to 5 for -10%, up to 8 for -15%).
