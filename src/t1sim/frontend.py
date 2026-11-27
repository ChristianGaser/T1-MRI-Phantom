"""The front end of t1sim: a real scan becomes a label map and tissue fractions.

The same three steps ``matlab/mri_simulate.m`` takes, with CAT-Surface's
unified segmentation in place of SPM's:

1. SANLM denoising, at the resolution of the scan (``cat_surf.vol_sanlm``);
2. unified segmentation with the six-class shooting templates
   (``t1prep.unified_seg.segment``, i.e. ``cat_surf.vol_unified_seg``), which
   gives the affine to the template, the bias-corrected image, the brain mask
   and a first label map on T1Prep's 0.5 mm working grid;
3. LAS, the local intensity normalisation (``t1prep.apply_LAS``), which maps
   the corrected brain onto the scale CSF = 1, GM = 2, WM = 3.

The label map is that normalised image, and the tissue fractions are its
triangular decomposition, as in the MATLAB version: the ground truth is a
function of the image it was made from and of nothing else.  No second
segmentation (AMAP, a network) takes part.

Everything is returned on the RAS-ordered working grid.  The intensities of
the three tissues, which the synthesis needs, are those of the mixture of
Gaussians the segmentation estimated; they are not fitted a second time.

The same holds for the intensity LAS takes for pure GM: it is the GM mean of
that mixture, as ``mri_simulate.m`` hands SPM's to ``cat_main_LASsimple``,
and not the peak LAS would estimate itself.  A voxel as bright as the GM of
the mixture then gets the label 2 and is rendered at that very intensity.
Left to itself LAS put pure GM 6% lower on a 0.75 mm MPRAGE (358 against
380), which moved the GM/WM boundary 0.18 mm outward and made WM 9% larger
than in the label of ``mri_simulate.m``; with the GM mean of the mixture the
boundaries of the two agree to 0.04 mm.  CSF and WM keep LAS's own estimates,
medians over the voxels the first label map calls pure; ``mri_simulate.m``
prefers a median for WM as well, so that white matter hyperintensities do not
pull that anchor down.
"""

from __future__ import annotations

import inspect
import os

import nibabel as nib
import numpy as np

from . import models

#: Order of the tissue fractions, and their values in the label map.
TISSUES = ("CSF", "GM", "WM")

#: Lowest value of the normalised image inside the brain.  0 is what lies
#: outside of it, and has no fractions.
DARKEST = 0.01


def _t1prep():
    """T1Prep's public unified segmentation and LAS, or an error that names what is missing."""
    try:
        import cat_surf
        from t1prep import apply_LAS, unified_seg
    except ImportError as exc:
        raise RuntimeError(
            "The front end of t1sim needs cat-surf >= 1.2.0 and a T1Prep with the "
            "public modules t1prep.unified_seg and t1prep.apply_LAS (newer than "
            f"0.9.1, or a source checkout): {exc}"
        ) from exc
    if "peaks" not in inspect.signature(apply_LAS).parameters:
        raise RuntimeError(
            "The front end of t1sim needs a T1Prep whose apply_LAS takes the tissue "
            "peaks (apply_LAS(..., peaks=...)): newer than 0.9.1, or a source checkout."
        )
    return cat_surf, apply_LAS, unified_seg


def _in_storage_order(data: np.ndarray, affine: np.ndarray, stored) -> nib.Nifti1Image:
    """A volume in RAS order as an image on the grid the scan is stored on.

    Args:
        data: (X, Y, Z) volume in the order of ``nib.as_closest_canonical``.
        affine: Its voxel-to-world matrix.
        stored: The image as it was loaded; its storage order and affine are
            taken.

    Returns:
        The volume with the axes and directions of *stored*.
    """
    there = nib.orientations.io_orientation(stored.affine)
    here = nib.orientations.io_orientation(affine)
    data = nib.orientations.apply_orientation(data, nib.orientations.ornt_transform(here, there))
    return nib.Nifti1Image(np.ascontiguousarray(data), stored.affine)


def fractions_from_label(label: np.ndarray) -> np.ndarray:
    """Triangular decomposition of a label map into CSF, GM and WM.

    ``1 - min(1, |label - c|)`` for c = 1, 2, 3, divided by their sum, as
    ``mri_simulate.m`` does it.  Between two tissues the fractions are their
    linear mixture.  Below 1 and above 3 only one fraction is left, and the
    division makes it 1: a brain voxel darker than CSF is CSF, one brighter
    than WM is WM.  From 4 on, and at 0, all fractions are zero: such a voxel
    is not brain, and the synthesis keeps the image there.  Inside the brain
    this happens to the brightest vessels, which :mod:`t1sim.cleanup` makes
    CSF; :func:`segment` keeps the darkest voxels of the brain above 0.

    Args:
        label: (X, Y, Z) array on the scale CSF = 1, GM = 2, WM = 3.

    Returns:
        (X, Y, Z, 3) float32, in the order of :data:`TISSUES`.
    """
    label = np.asarray(label, dtype=np.float32)
    out = np.empty(label.shape + (3,), dtype=np.float32)
    for k in range(3):
        np.subtract(label, k + 1, out=out[..., k])
        np.abs(out[..., k], out=out[..., k])
        np.clip(1 - out[..., k], 0, 1, out=out[..., k])
    total = out.sum(axis=-1, keepdims=True)
    np.divide(out, total, out=out, where=total > 0)
    return out


def label_from_fractions(fractions: np.ndarray) -> np.ndarray:
    """The label map of tissue fractions: CSF + 2 GM + 3 WM."""
    return fractions[..., 0] + 2 * fractions[..., 1] + 3 * fractions[..., 2]


def tissue_means(params: dict) -> dict:
    """Intensity of CSF, GM and WM in the bias-corrected image.

    Read from the mixture of Gaussians of the unified segmentation
    (``params`` of ``cat_surf.vol_unified_seg``), as ``mri_simulate.m`` reads
    them from SPM's: for GM and WM the mean of their Gaussians weighted by
    the mixing proportions, for CSF its darkest Gaussian -- a second CSF
    Gaussian regularly covers GM, and their mean would be far too bright.

    Args:
        params: ``lkp`` (class of every Gaussian, 0-based in the order of
            ``classes``), ``mg`` (mixing proportions) and ``mn`` (means).

    Returns:
        ``{"CSF": ..., "GM": ..., "WM": ...}``.
    """
    lkp = np.asarray(params["lkp"])
    mg = np.asarray(params["mg"], dtype=np.float64)
    mn = np.asarray(params["mn"], dtype=np.float64)
    classes = list(params["classes"])
    means = {}
    for tissue in TISSUES:
        of_class = lkp == classes.index(tissue)
        if not of_class.any():
            raise ValueError(f"the segmentation has no Gaussian for {tissue}")
        if tissue == "CSF":
            means[tissue] = float(mn[of_class].min())
        else:
            means[tissue] = float(np.sum(mg[of_class] * mn[of_class]) / np.sum(mg[of_class]))
    return means


def segment(t1, template_dir=None, verbose: bool = False, deformation: bool = False) -> dict:
    """Label map, tissue fractions and corrected image of a real T1w scan.

    Args:
        t1: The scan, a file name or a nibabel image, in any storage order.
        template_dir: Folder with the six-class shooting templates; default:
            :func:`t1sim.models.ensure_models`, which downloads them once.
        verbose: Print the stages of the segmentation.
        deformation: Also return the deformation to the templates, which
            :mod:`t1sim.atlas` needs to bring an atlas onto the working grid.

    Returns:
        A dict, all volumes on the RAS-ordered working grid (``grid``):

        - ``label``: (X, Y, Z) float32, 0 outside the brain, CSF = 1, GM = 2,
          WM = 3 with partial volumes in between;
        - ``fractions``: (X, Y, Z, 3) float32, CSF, GM and WM
          (:func:`fractions_from_label` of the normalised image).  They sum
          to 1 in the brain, except in voxels of 4 or more in ``normalised``,
          the brightest vessels, where they are 0;
        - ``normalised``: (X, Y, Z) float32, the LAS-normalised brain on the
          scale of the label, up to 5 for what is brighter than WM; its 2 is
          the GM mean of the mixture (``means["GM"]``);
        - ``corrected``: the denoised and bias-corrected scan, whole head;
        - ``mask``: the brain mask, boolean;
        - ``means``: :func:`tissue_means`, on the scale of ``corrected``;
        - ``grid``: the working grid as an image; ``affine``, ``world``,
          ``native_mat`` and ``native_shape``: how it relates to the scan
          (see ``t1prep.unified_seg.segment``);
        - ``scan``: the denoised scan as an image on its own grid, in the
          storage order of the input.  ``simulate.acquire(grid="scan")``
          writes onto that grid and takes from it what the scan shows beyond
          the working grid;
        - ``params``: the parameters the segmentation estimated;
        - with *deformation*: ``iy``, (X, Y, Z, 3) float32, for every voxel
          the template voxel it maps to (0-based, NaN where undefined), and
          ``template_mat``, the voxel-to-world matrix of the templates.
    """
    cat_surf, apply_LAS, unified_seg = _t1prep()
    if template_dir is None:
        template_dir = models.ensure_models(verbose=verbose)

    stored = nib.load(os.fspath(t1)) if isinstance(t1, (str, os.PathLike)) else t1
    # RAS storage order: SANLM and the estimation follow the stored lattice
    image = nib.as_closest_canonical(stored)
    data = np.nan_to_num(np.squeeze(image.get_fdata(dtype=np.float32)))
    if data.ndim != 3:
        raise ValueError(f"the scan must be a single 3-D volume, not {data.shape}")
    raw = nib.Nifti1Image(data, image.affine)
    denoised = nib.Nifti1Image(
        np.asarray(cat_surf.vol_sanlm(data), dtype=np.float32), image.affine)
    scan = _in_storage_order(np.asarray(denoised.dataobj), image.affine, stored)

    extra = {}
    if deformation:
        if "deformation" not in inspect.signature(unified_seg.segment).parameters:
            raise RuntimeError(
                "The deformation needs a T1Prep whose unified_seg.segment takes "
                "deformation=True (newer than 0.9.2, or a source checkout).")
        extra["deformation"] = True
    seg = unified_seg.segment(denoised, raw, template_dir=template_dir, verbose=verbose, **extra)
    mask = unified_seg.brain_mask(seg.pop("prob"))
    _, brain, first_label = unified_seg.working_grid_maps(seg, mask)
    means = tissue_means(seg["params"])
    # Pure GM is where the mixture has it; CSF and WM are left to LAS
    las = apply_LAS(brain, first_label, verbose=verbose, peaks=(None, means["GM"], None))
    normalised = 3.0 * np.asarray(las.dataobj, dtype=np.float32)
    # A brain voxel as dark as the background would get no fractions: it is CSF
    np.maximum(normalised, DARKEST, out=normalised)
    normalised[~mask] = 0.0

    fractions = fractions_from_label(normalised)
    maps = {"iy": seg["iy"], "template_mat": seg["template_mat"]} if deformation else {}
    return {
        **maps,
        "label": label_from_fractions(fractions),
        "fractions": fractions,
        "normalised": normalised,
        "corrected": np.nan_to_num(seg["corrected"]).astype(np.float32, copy=False),
        "mask": mask,
        "means": means,
        "grid": seg["grid"],
        "affine": seg["affine"],
        "world": seg["world"],
        "native_mat": seg["native_mat"],
        "native_shape": seg["native_shape"],
        "scan": scan,
        "params": seg["params"],
    }
