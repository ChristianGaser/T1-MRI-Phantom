"""Operations on volumes that several modules share.

Smoothing as SPM does it, distances and the morphology built on them, with
all sizes in mm.  ``matlab/mri_simulate.m`` takes these from SPM
(``spm_smooth``) and CAT12 (``cat_bwdist``, ``cat_vol_morph``); the functions
here follow their definitions, so that a port of an algorithm keeps its
parameters.

Volumes are 3-D numpy arrays indexed ``[x, y, z]``; ``vx`` is the voxel size
in mm, a scalar or one value per axis.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.special import erf


def voxel_size(vx) -> np.ndarray:
    """*vx* as three floats."""
    vx = np.atleast_1d(np.asarray(vx, dtype=np.float64))
    if vx.size == 1:
        vx = np.repeat(vx, 3)
    if vx.size != 3 or np.any(vx <= 0):
        raise ValueError(f"the voxel size must be one or three positive values, not {vx}")
    return vx


def grid_voxel_size(affine) -> np.ndarray:
    """Voxel size in mm of a grid, from its voxel-to-world matrix."""
    return np.sqrt((np.asarray(affine, dtype=np.float64)[:3, :3] ** 2).sum(0))


def bounding_box(mask: np.ndarray, margin: int = 0) -> tuple:
    """Slices of the smallest box that holds *mask*, widened by *margin* voxels."""
    if not mask.any():
        raise ValueError("the mask is empty")
    box = []
    for axis in range(3):
        other = tuple(a for a in range(3) if a != axis)
        index = np.where(mask.any(axis=other))[0]
        box.append(slice(max(int(index[0]) - margin, 0),
                         min(int(index[-1]) + margin + 1, mask.shape[axis])))
    return tuple(box)


def _smoothing_kernel(fwhm: float) -> np.ndarray:
    """SPM's smoothing kernel for a FWHM in voxels (``spm_smoothkern``, linear).

    A Gaussian convolved with the triangle of a linear interpolation, sampled
    out to six standard deviations and scaled to a sum of one.
    """
    sd = fwhm / np.sqrt(8 * np.log(2))
    x = np.arange(-round(6 * sd), round(6 * sd) + 1, dtype=np.float64)
    var = sd ** 2 + np.finfo(np.float64).eps
    w1, w2, w3 = 0.5 * np.sqrt(2 / var), -0.5 / var, np.sqrt(var / 2 / np.pi)
    kernel = (0.5 * (erf(w1 * (x + 1)) * (x + 1) + erf(w1 * (x - 1)) * (x - 1) - 2 * erf(w1 * x) * x)
              + w3 * (np.exp(w2 * (x + 1) ** 2) + np.exp(w2 * (x - 1) ** 2) - 2 * np.exp(w2 * x ** 2)))
    kernel[kernel < 0] = 0
    return kernel / kernel.sum()


def smooth(volume: np.ndarray, fwhm, vx=1.0) -> np.ndarray:
    """Gaussian smoothing with SPM's kernel (``spm_smooth``).

    Args:
        volume: The volume.
        fwhm: Full width at half maximum in mm, a scalar or one value per axis.
        vx: Voxel size in mm.

    Returns:
        The smoothed volume, float32.  Outside the volume it is taken as 0.
    """
    fwhm = voxel_size(fwhm) / voxel_size(vx)
    out = np.asarray(volume, dtype=np.float32)
    for axis in range(3):
        if fwhm[axis] > 0:
            out = ndimage.correlate1d(out, _smoothing_kernel(fwhm[axis]).astype(np.float32),
                                      axis=axis, mode="constant", cval=0.0)
    return out


def distance(mask: np.ndarray, vx=1.0) -> np.ndarray:
    """Euclidean distance in mm of every voxel to the nearest voxel of *mask*.

    Exact, between voxel centres, and 0 inside the mask (CAT12's
    ``cat_bwdist``).  Infinite everywhere if the mask is empty.
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.full(mask.shape, np.inf, dtype=np.float32)
    return ndimage.distance_transform_edt(~mask, sampling=voxel_size(vx)).astype(np.float32)


def nearest_index(mask: np.ndarray, vx=1.0) -> tuple:
    """For every voxel the index of the nearest voxel of *mask* (three index arrays)."""
    index = ndimage.distance_transform_edt(~np.asarray(mask, dtype=bool), sampling=voxel_size(vx),
                                           return_distances=False, return_indices=True)
    return tuple(index)


def dilate(mask: np.ndarray, mm: float, vx=1.0) -> np.ndarray:
    """All voxels within *mm* of the mask (``cat_vol_morph`` 'dd')."""
    mask = np.asarray(mask, dtype=bool)
    if mm <= 0 or not mask.any():
        return mask.copy()
    # only a box around the mask can change
    box = bounding_box(mask, int(np.ceil(mm / voxel_size(vx).min())) + 1)
    out = np.zeros(mask.shape, dtype=bool)
    out[box] = distance(mask[box], vx) <= mm
    return out


def erode(mask: np.ndarray, mm: float, vx=1.0) -> np.ndarray:
    """All voxels of the mask further than *mm* from its outside ('de')."""
    mask = np.asarray(mask, dtype=bool)
    if mm <= 0 or mask.all():
        return mask.copy()
    box = bounding_box(mask, 1) if mask.any() else None
    out = np.zeros(mask.shape, dtype=bool)
    if box is not None:
        out[box] = distance(~mask[box], vx) > mm
    return out


def closing(mask: np.ndarray, mm: float, vx=1.0) -> np.ndarray:
    """Dilation by *mm*, then erosion by *mm* ('dc')."""
    mask = np.asarray(mask, dtype=bool)
    if mm <= 0 or not mask.any():
        return mask.copy()
    pad = int(np.ceil(mm / voxel_size(vx).min())) + 1
    padded = np.pad(mask, pad)          # the dilated mask must not touch the border
    return erode(dilate(padded, mm, vx), mm, vx)[pad:-pad, pad:-pad, pad:-pad]


def opening(mask: np.ndarray, mm: float, vx=1.0) -> np.ndarray:
    """Erosion by *mm*, then dilation by *mm* ('do')."""
    return dilate(erode(mask, mm, vx), mm, vx)


def largest_component(mask: np.ndarray) -> np.ndarray:
    """The largest 6-connected component of the mask ('l')."""
    labels, n = ndimage.label(np.asarray(mask, dtype=bool))
    if n < 2:
        return labels > 0
    return labels == (np.bincount(labels.ravel())[1:].argmax() + 1)


def _ball(radius: int) -> np.ndarray:
    axis = np.arange(-radius, radius + 1)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    return (x ** 2 + y ** 2 + z ** 2) <= radius ** 2


def grey_close(volume: np.ndarray, radius: int) -> np.ndarray:
    """Grey-scale closing with a ball of *radius* voxels ('gc').

    The maximum over the ball, then the minimum over it: bright structures
    thinner than the ball's diameter that nearly touch are joined, and dark
    gaps of that size inside a bright region are filled.
    """
    ball = _ball(int(radius))
    out = ndimage.maximum_filter(np.asarray(volume, dtype=np.float32), footprint=ball, mode="nearest")
    return ndimage.minimum_filter(out, footprint=ball, mode="nearest")


def box_mean(volume: np.ndarray) -> np.ndarray:
    """Mean over the 3 x 3 x 3 neighbourhood (MATLAB's ``smooth3`` with its defaults)."""
    return ndimage.uniform_filter(np.asarray(volume, dtype=np.float32), size=3, mode="nearest")


def percentile(volume: np.ndarray, p: float) -> float:
    """Percentile of the finite values, as MATLAB's ``prctile`` defines it."""
    values = np.asarray(volume)
    values = values[np.isfinite(values)]
    return float(np.percentile(values, p, method="hazen"))
