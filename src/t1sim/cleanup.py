"""Cleanup of the ground truth: what looks like grey matter and is not.

The label map is made from the intensity of the scan, so everything as
bright as grey matter is grey matter in it: blood vessels, the dura, and the
voxels around the ventricles that mix CSF and white matter.  This module
finds the three and corrects the ground truth.

The image is deliberately not corrected.  It is rendered from the fractions
as they were (:mod:`t1sim.synthesis` runs first), so a vessel keeps its
intensity in the simulated scan while the ground truth calls it CSF.  That is
what a segmentation has to learn and what a validation has to test: a
synthesis whose label describes every structure of its own image poses a
problem that is solved by reading the intensity.

**Blood vessels** come from two sources.  T1Prep's detection
(``t1prep.vessels``), its port of CAT12's, admits a voxel that is bright, has
CSF around it, lies away from the white matter body and where a prior built
from MRA scans allows it.  It is careful not to touch thin gyral white
matter, and on a scan full of bright vessels it found a fifth of the voxels that
are clearly brighter than white matter.  Those are added directly: LAS brings
white matter to 3 everywhere, so a voxel well above that is no brain tissue.
Voxels above :data:`HYPER_SEED` are taken, together with the voxels above
:data:`HYPER_GROW` that are connected to them, outside the white matter
body.  Measured on that scan and on a scan without visible vessels, the
seeds are 6.3 ml against 0.05 ml.

The brightest of these voxels have no label at all.  From an intensity of 4
on, the decomposition of the front end gives a brain voxel neither tissue nor
CSF (:func:`t1sim.frontend.fractions_from_label`), so that the image keeps
the scan there: 3.7 ml on the scan with the vessels, as holes of zeros inside
the ground truth.  They are vessels by the same rule and become CSF with the
others.

**Dura** cannot be told from cortex by its shape or by its distance to the
white matter; both are thin sheets of grey matter intensity.  What differs is
that the cortex hangs on the white matter while the dura lies beyond the
subarachnoid CSF, so a path from the white matter to the dura has to leave
the tissue.  The white matter is therefore grown through the tissue only, and
grey matter it does not reach within a few mm is dura.  Structures whose
shape makes a long path normal are spared: the cerebellum, whose white matter
is a thin tree, and the deep grey matter.

**Periventricular partial volume** is found as the level 3 cleanup of CAT12's
``cat_main_cleanup`` finds it: next to a ventricle, the brainstem or the
corpus callosum, not in the basal ganglia or the thalamus, between pure white
matter and pure CSF, and thin.  Such a voxel is decomposed again as the
mixture of CSF and WM that has its intensity in the rendered image.  That is
not the mixture with its label: GM is not half way between CSF and WM (148,
380 and 501 on the scan this was measured on), so a voxel as bright as GM is
34% CSF and 66% WM, and its label becomes 2.3.
"""

from __future__ import annotations

import numpy as np

from . import atlas
from ._volume import (bounding_box, box_mean, dilate, erode, grid_voxel_size, opening)

#: Order of the fractions: CSF, GM, WM.
_CSF, _GM, _WM = 0, 1, 2

#: Intensity on the scale of the label (WM = 3) above which a voxel of the
#: brain is a vessel, and above which a voxel connected to such a one is.
HYPER_SEED = 3.5
HYPER_GROW = 3.2
#: Intensity from which on the front end gives a voxel no label and no fractions.
NO_LABEL = 4.0


def _in(regions_on_grid: np.ndarray, numbers) -> np.ndarray:
    return np.isin(regions_on_grid, numbers)


def detect(label: np.ndarray, normalised: np.ndarray, grid_affine: np.ndarray,
           regions_on_grid: np.ndarray, vessels: float = 0.5, dura: float = 4.0,
           pve: bool = True) -> dict:
    """Find vessels, dura and periventricular partial volume in a label map.

    Args:
        label: (X, Y, Z) label map, CSF = 1, GM = 2, WM = 3.
        normalised: (X, Y, Z) LAS-normalised image on the scale of the label
            and 0 outside the brain (``frontend.segment(...)["normalised"]``);
            the vessel detection reads the intensities from it.
        grid_affine: Voxel-to-world matrix of the working grid.
        regions_on_grid: (X, Y, Z) Neuromorphometrics atlas on the grid
            (``atlas.to_grid("Neuromorphometrics", iy, template_mat)``).
        vessels: Strength of the vessel detection in 0..1; 0 switches it off.
            A voxel is a vessel where T1Prep's vessel weight exceeds
            ``1 - vessels``, or where it is clearly brighter than white
            matter (see the module notes).
        dura: Length in mm of the path through the tissue beyond which grey
            matter is dura; 0 switches it off.  Lower removes more.
        pve: Detect the periventricular partial volume.

    Returns:
        ``{"vessels": ..., "dura": ..., "pve": ...}``, boolean masks.
    """
    if not 0 <= vessels <= 1:
        raise ValueError("vessels must lie in 0..1")
    if dura < 0:
        raise ValueError("dura is a length in mm and cannot be negative")
    vx = grid_voxel_size(grid_affine)
    shape = label.shape
    found = {name: np.zeros(shape, dtype=bool) for name in ("vessels", "dura", "pve")}
    box = bounding_box(label > 0.5, int(np.ceil(4 / vx.min())))
    lab = label[box]
    regions = regions_on_grid[box]

    if vessels > 0:
        from t1prep import vessels as t1prep_vessels

        weight = t1prep_vessels.vessel_weight(
            np.ascontiguousarray(normalised, dtype=np.float32) / 3.0,
            np.ascontiguousarray(label, dtype=np.float32), tuple(float(v) for v in vx),
            bv_prior=t1prep_vessels.blood_vessel_prior(grid_affine, shape),
            protect=t1prep_vessels.protected_regions(grid_affine, shape))
        import cat_surf

        norm = normalised[box]
        body = erode(lab >= 2.5, 1.5, vx)                # hyperintense white matter is no vessel
        # the image is 0 outside the brain, so every bright voxel is one of the brain
        seeds = (norm > HYPER_SEED) & ~body
        bright = cat_surf.vol_dilate_geodesic(seeds, region=(norm > HYPER_GROW) & ~body, iterations=4)
        found["vessels"] = weight > 1 - vessels
        found["vessels"][box] |= np.asarray(bright, dtype=bool)
        # A vessel the label already calls CSF needs no correction.  One that
        # is too bright to have a label does: it is a hole in the ground truth.
        found["vessels"] &= (label > 1.5) | (normalised >= NO_LABEL)

    if dura > 0:
        import cat_surf

        tissue = lab > 1.5
        steps = int(round(dura / vx.mean()))
        reached = cat_surf.vol_dilate_geodesic(lab > 2.5, region=tissue, iterations=steps)
        spared = dilate(_in(regions, atlas.BASAL_GANGLIA + atlas.THALAMUS + atlas.HIPPOCAMPUS
                            + atlas.BRAINSTEM), 2.0, vx) | _in(regions, atlas.CEREBELLUM)
        found["dura"][box] = tissue & (lab < 2.5) & ~np.asarray(reached, dtype=bool) & ~spared

    if pve:
        ventricle = _in(regions, atlas.VENTRICLES)
        near_ventricle = dilate(ventricle, 2.0, vx)
        brainstem = _in(regions, atlas.BRAINSTEM) & (lab > 2)
        cortex = regions >= atlas.CORTEX_FROM
        right = (cortex & (regions % 2 == 0)) | (regions == atlas.CEREBRAL_WM[0])
        left = (cortex & (regions % 2 == 1)) | (regions == atlas.CEREBRAL_WM[1])
        # the corpus callosum: where the two hemispheres meet, next to the ventricle
        callosum = dilate(left, 3.0, vx) & dilate(right, 3.0, vx) & near_ventricle
        deep_grey = box_mean(_in(regions, atlas.BASAL_GANGLIA + atlas.THALAMUS)) > 0.3
        between = (lab > 1) & (lab < 3)
        thin = box_mean(between & ~opening(between, 1.5, vx)) > 0.1
        found["pve"][box] = ((dilate(ventricle | brainstem, 2.0, vx) | callosum) & ~deep_grey
                             & dilate(lab > 2.9, 2.0, vx) & dilate(lab < 1.1, 2.0, vx)
                             & between & thin)
    return found


def apply(fractions: np.ndarray, label: np.ndarray, to_csf: np.ndarray = None,
          pve: np.ndarray = None, means=None) -> tuple:
    """Correct the ground truth in the detected voxels.

    Args:
        fractions: (X, Y, Z, 3) fractions of CSF, GM and WM.
        label: (X, Y, Z) label map.
        to_csf: Voxels that become CSF: vessels and dura.
        pve: Voxels that become the mixture of CSF and WM with the intensity
            they are rendered at: ``csf = (WM - intensity) / (WM - CSF)``
            with ``intensity = fractions @ means``, ``wm = 1 - csf``, no GM.
        means: The tissue intensities the image is rendered with: a dict
            with ``CSF``, ``GM`` and ``WM``, or three values in that order.
            Needed with *pve*.

    Returns:
        ``(fractions, label)``, new arrays; the label is that of the new
        fractions.
    """
    fractions = np.array(fractions, dtype=np.float32)
    label = np.array(label, dtype=np.float32)
    if pve is not None and pve.any():
        if means is None:
            raise ValueError("the periventricular correction needs the tissue intensities (means)")
        if isinstance(means, dict):
            means = [means[t] for t in ("CSF", "GM", "WM")]
        means = np.asarray(means, dtype=np.float32)
        if not means[_CSF] < means[_WM]:
            raise ValueError("white matter has to be brighter than CSF")
        intensity = fractions[pve] @ means
        csf = np.clip((means[_WM] - intensity) / (means[_WM] - means[_CSF]), 0, 1)
        fractions[pve] = np.stack([csf, np.zeros_like(csf), 1 - csf], axis=-1)
        label[pve] = csf + 3 * (1 - csf)
    if to_csf is not None and to_csf.any():
        fractions[to_csf] = (1.0, 0.0, 0.0)
        label[to_csf] = 1.0
    return fractions, label


def clean(seg: dict, vessels: float = 0.5, dura: float = 4.0, pve: bool = True,
          fractions: np.ndarray = None, label: np.ndarray = None) -> dict:
    """Detect and correct in one step.

    Args:
        seg: Result of ``frontend.segment(..., deformation=True)``.
        vessels, dura, pve: See :func:`detect`.
        fractions, label: The ground truth to clean if it is no longer that
            of *seg*, e.g. after a change of the anatomy.

    Returns:
        A dict with the cleaned ``fractions`` and ``label``, the ``masks`` of
        :func:`detect` and ``volumes_ml``, the volume of each mask in ml of
        the working grid.
    """
    for key in ("iy", "template_mat"):
        if key not in seg:
            raise ValueError("the cleanup needs the deformation: frontend.segment(..., deformation=True)")
    fractions = seg["fractions"] if fractions is None else fractions
    label = seg["label"] if label is None else label
    affine = seg["grid"].affine
    regions_on_grid = atlas.to_grid("Neuromorphometrics", seg["iy"], seg["template_mat"])
    masks = detect(label, seg["normalised"], affine, regions_on_grid, vessels, dura, pve)
    new_fractions, new_label = apply(fractions, label, masks["vessels"] | masks["dura"], masks["pve"],
                                     means=seg["means"])
    voxel_ml = float(np.prod(grid_voxel_size(affine))) / 1000
    return {"fractions": new_fractions, "label": new_label, "masks": masks,
            "volumes_ml": {name: float(mask.sum()) * voxel_ml for name, mask in masks.items()}}
