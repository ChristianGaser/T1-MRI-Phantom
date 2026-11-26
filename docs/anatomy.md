# Changing the anatomy

Part of the [documentation of t1sim](../README.md#documentation).

`t1sim.anatomy` changes the fractions before the image is rendered from them, so the simulated image shows an anatomy that is known: a cortex of a chosen thickness, or a region with atrophy. Both are options of `simulate.ground_truth`.

- [Cortical thickness](#cortical-thickness)
- [Atrophy](#atrophy)
- [What `change` means](#what-change-means)
- [Which kind of atrophy](#which-kind-of-atrophy)

All lengths and volumes on this page are those of the output grid, on which the head has the size of the template ([The working grid](frontend.md#the-working-grid)).

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

The [cleanup](cleanup.md) leaves the rebuilt cortex alone and acts only on the structures that kept their labels.

Measured:

| Object | Asked for | Result |
| --- | --- | --- |
| Sphere at 0.5 mm, thickness along rays through the centre | 1.5 / 2.5 / 3.5 mm | 1.50 / 2.49 / 3.49 mm |
| The same, from the volumes of GM and WM | 1.5 / 2.5 / 3.5 mm | 1.46 / 2.47 / 3.47 mm |
| A larger sphere (white matter radius 22 mm), projection-based thickness of cat-surf without its correction, median | 1.5 / 2.5 / 3.5 mm | 1.41 / 2.40 / 3.40 mm |
| Spheres of that size written down analytically, the same measure | 1.5 / 2.5 / 3.5 mm | 1.47 / 2.43 / 3.42 mm |
| Real brain (0.75 mm MPRAGE), the same measure with its default correction, median and 10 to 90% | 2.5 mm | 2.66 mm, 2.55 to 2.71 |
| The same | 1.5 mm | 1.66 mm, 1.62 to 1.73 |
| The same, before the simulation | - | 2.89 mm, 2.27 to 3.48 |

What this shows:

- **The band has the thickness asked for**, to 0.01 mm along rays on a sphere.
- **Its boundaries are soft.** The partial volume at both boundaries is 1.15 mm wide on the 0.5 mm grid, because the label map is smoothed first. The ground truth of a simulated thickness is smoother than its grid.
- **The projection-based thickness reads less than the truth on a ball, simulated or not.** By default it adds 0.25 mm (`correct_thickness`), a correction for the label maps of cat-surf's own pipeline; with it the real brain reads 2.66 mm for 2.5. Without it, it reads 0.03 to 0.09 mm too little on analytic spheres, whatever the width of their partial volume (0.5 to 1.5 mm were tried), and 0.02 to 0.06 mm less again on the simulated ones. On the folded [geometric phantom](phantom.md) the same measure is within 0.02 mm. So most of the offset is the measure on a convex ball; the last 0.02 to 0.06 mm are not explained.

A run takes 7 to 8 minutes on the 0.5 mm grid, most of it in 75 distance transforms. For a thickness that is known by construction, without a real brain, see the [geometric phantom](phantom.md).

## Atrophy

```python
truth = simulate.ground_truth(scan, atrophy={"atlas": "Hammers", "regions": [28], "change": -0.10})
truth["jacobian"]                               # Jacobian determinant of the deformation
truth["parameters"]["Atrophy"][0]               # what was asked for and what was reached
```

The tissue of a region changes its volume by `change`, by a smooth deformation that is applied to the ground truth, to the fractions of the image and to the scan outside the brain. Boundaries move, and the Jacobian determinant of the deformation says where volume was lost. Two kinds:

- **`mode="thinning"`** (default). Every point of the cortex of the region moves towards the white matter surface in proportion to its distance from it. The cortex keeps its area and gets thinner; the CSF beyond it widens. The white matter surface stays.
- **`mode="volume"`**. The tissue of the region is the sink and the CSF within `reach` the source of a curl-free flow. The tissue shrinks the same in all directions.

The options are listed in [Simulating images](simulation.md#options-of-the-anatomy). A positive `change` thickens or enlarges.

- **The deformation is made from the cleaned ground truth**, which knows what is no brain tissue; made from the fractions of the image, it drew the dura in front of the gyrus in by 0.4 mm instead of 0.2 mm, further than the pial surface itself.
- **Several regions are deformed one after the other.** From the second on, the atlas has moved with the anatomy, and the Jacobian determinants are chained: the second is read where the first took every voxel. Their plain product is wrong where two regions meet. For two neighbouring gyri it was off by 0.04 (RMS) in the 9.6 ml where both act, and by more than 0.05 in 11% of those voxels.
- **The Jacobian determinant lies on the anatomy before the atrophy.** Voxel by voxel it is the factor by which the tissue that was there changed its volume; below 1 is a loss. It belongs to a simulation of the same scan without atrophy, the baseline of a longitudinal comparison, and not to the label map of the atrophied anatomy, in which the same tissue lies up to half a mm away. Its sidecar says so.

## What `change` means

`change` is the relative change of the volume of the tissue that was in the region, followed with that tissue: the integral of the Jacobian determinant over it. The strength of the deformation is adjusted until that number is reached, and the sidecar holds it as `Achieved`. Two other ways of counting are just as legitimate and give other numbers, because tissue crosses the border of the region and because neighbouring cortex is deformed a little as well. The sidecar holds them too:

| Entry of `Atrophy` in the sidecar | Counts | Thinning, -10% | Volume, -10% |
| --- | --- | --- | --- |
| `Achieved` | The tissue that was in the region, wherever it is afterwards (Jacobian) | -10.0% | -10.0% |
| `ChangeInRegion` | The tissue inside the fixed region, before and after | -9.7% | -6.1% |
| `TissueChangeMl`, GM, relative to the GM of the region | All grey matter of the ground truth, before and after | -3.0 ml = -11.1% | -2.0 ml = -7.6% |

(Left middle frontal gyrus, Hammers 28, with 26.7 ml of grey matter, on a 0.75 mm MPRAGE.) A method that measures the volume of a region in the same atlas should find `ChangeInRegion`; one that measures the whole brain should find `TissueChangeMl`. With the volume flow the three differ most: neighbouring tissue moves into the region, and the cortex next to it gains a little, because the CSF that takes up the space is partly the partial volume of that cortex. In both kinds the white matter next to the region loses a little as well (0.5 and 0.8 ml), and the three tissues together change by 0.2 and 1.0 ml of 2037 ml, because the outer border of the brain is moved slightly.

The exact ground truth of any measure is the difference between the fractions of the two simulations, with and without the atrophy: both images are rendered from their own fractions.

## Which kind of atrophy

`matlab/mri_simulate.m` multiplies the CSF fraction inside the region by a factor and renormalises. The three were compared on the left middle frontal gyrus (Hammers 28) of a 0.75 mm MPRAGE, on the 0.5 mm grid, all starting from the cleaned ground truth. The thickness is the median projection-based thickness of cat-surf in the region, with the 0.25 mm it adds by default.

| Method | GM in the fixed region | All GM, relative to the region | Thickness | Pial boundary | WM boundary | Intensity of the cortex | GM in the region at 1 mm |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Before | 26.7 ml | - | 2.88 mm | - | - | - | - |
| Thinning, -10% | -9.7% | -11.1% | -0.26 mm (-9.9%) | -0.35 mm | -0.05 mm | -1.5% | -9.7% |
| Volume, -10% | -6.1% | -7.6% | -0.09 mm (-3.4%) | -0.18 mm | -0.01 mm | -0.5% | -6.1% |
| Scaling of CSF, factor 2 | -9.5% | -9.5% | -0.32 mm (-12.2%) | -0.29 mm | 0 | -2.1% | -10.6% |
| Scaling of CSF, factor 3 | -15.8% | -15.8% | -0.53 mm (-20.2%) | -0.48 mm | 0 | -3.3% | -17.5% |

The percentages of the thickness refer to the 2.63 mm that are left without the 0.25 mm. The boundary shifts are the change of the volume beyond a boundary divided by its area. The last column is the same anatomy at 1 mm voxels: the deformed fractions averaged to 1 mm, and for the scaling the factor applied to the fractions at 1 mm, which is how that method would be used there.

What this shows:

- **Scaling does thin the cortex**, because nearly every cortical voxel of a label map made from intensities holds some CSF. But its effect is set by a factor without a unit, it depends on the voxel size, and it darkens the cortex, since the voxels lose grey matter where they are instead of the boundary moving. It gives no deformation to validate against. It is kept as `anatomy.atrophy_by_scaling` for comparison only.
- **A change of volume** shrinks the tissue the same in all directions, so it thins the cortex by only about a third of the volume change, and in the fixed region less grey matter is lost than was asked for, because neighbouring tissue moves in.
- **Thinning** is the model of cortical atrophy: the loss of volume is a loss of thickness, the same at any voxel size, with the intensity of the cortex nearly kept.

Limits of the thinning: where two banks of a sulcus touch without CSF between them, they cannot move apart, and the cortex thins less there; the thin CSF of a narrow sulcus is stretched instead, so the Jacobian determinant is large in single voxels (up to 5 for -10%).
