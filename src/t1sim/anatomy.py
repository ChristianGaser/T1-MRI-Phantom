"""Changes of the anatomy: a cortex of a chosen thickness, and atrophy.

Both change the ground truth before the image is rendered from it, so the
simulated image shows an anatomy that is known.

**Cortical thickness** (:func:`thickness`) is a port of ``simulate_thickness``
of ``matlab/mri_simulate.m``.  The cortex is rebuilt as the band of a chosen
thickness around the white matter: its outer boundary is a level set of the
Euclidean distance to the white matter surface.  Partial volumes come from
repeating that with the boundaries shifted by 15 sub-voxel offsets and
averaging the hard labels.  Non-cortical structures keep their labels, and
no band is grown around the ventricles or the corpus callosum.

**Atrophy** (:func:`atrophy`) is not a port.  ``mri_simulate.m`` multiplies
the CSF fraction inside a region and renormalises.  No point of the anatomy
moves, but nearly every cortical voxel of a label map made from intensities
holds some CSF, so the cortex does read thinner afterwards (12% for a factor
of 2 on the scan this was measured on).  The effect is set by a factor
without a unit, depends on the voxel size (a loss of 9.5% of the grey
matter at 0.5 mm was one of 10.6% at 1 mm), darkens the cortex, and leaves
no deformation to validate against.  That method is kept as
:func:`atrophy_by_scaling` for comparison.

:func:`atrophy` deforms the anatomy instead.  A smooth deformation is
integrated from a velocity field and applied to the fractions and to the
image: boundaries move, the cortex thins and the sulci widen, and the
Jacobian determinant of the deformation is the ground truth of the volume
change.  The velocity either scales the distance of the cortex from the
white matter surface (``"thinning"``) or is the curl-free field whose sinks
are the tissue of the region and whose sources are the CSF around it
(``"volume"``).

A deformation is described by three maps, all in voxels of the grid:

- ``displacement`` (backward): where every voxel of the new anatomy comes
  from, relative to itself.  It is what resamples a volume (:func:`warp`).
- ``forward``: where every voxel of the old anatomy goes, relative to itself.
- ``jacobian``: the Jacobian determinant of the forward map, so it lies on
  the old anatomy: the factor by which the tissue that was at a voxel changes
  its volume.  :func:`compose` chains the last two for deformations that
  follow each other.
"""

from __future__ import annotations

import numpy as np
from scipy import fft as sp_fft
from scipy import ndimage
from scipy.optimize import brentq

from . import atlas
from ._volume import (bounding_box, closing, dilate, distance, erode, grey_close, grid_voxel_size,
                      largest_component, smooth)

_CSF, _GM, _WM = 0, 1, 2

#: Sub-voxel offsets of the boundaries, in units of the label, whose hard
#: labels are averaged into partial volumes.
PVE_OFFSETS = np.linspace(-0.25, 0.25, 15)


def _round(values: np.ndarray) -> np.ndarray:
    """Round half up, as MATLAB rounds the non-negative labels."""
    return np.floor(values + 0.5)


# ------------------------------------------------------------------- thickness
def thickness(fractions: np.ndarray, label: np.ndarray, regions_on_grid: np.ndarray,
              grid_affine: np.ndarray, mm, verbose: bool = False, return_rebuilt: bool = False) -> tuple:
    """Give the cortex a chosen thickness.

    Args:
        fractions: (X, Y, Z, 3) fractions of CSF, GM and WM.
        label: (X, Y, Z) label map, CSF = 1, GM = 2, WM = 3.
        regions_on_grid: (X, Y, Z) Neuromorphometrics atlas on the grid.
        grid_affine: Voxel-to-world matrix of the grid.
        mm: The thickness in mm: one value, or three for the mainly occipital
            parcels, the rest and the mainly frontal parcels.
        verbose: Print the progress.
        return_rebuilt: Also return where the result is the rebuilt anatomy.

    Returns:
        ``(fractions, label)`` with the new cortex; with *return_rebuilt*
        also an (X, Y, Z) boolean map of the brain voxels whose fractions
        come from the rebuilt bands and no longer from the scan.  The other
        brain voxels, in the structures that are kept and deep in the white
        matter, have the fractions they had.
    """
    mm = np.atleast_1d(np.asarray(mm, dtype=np.float64))
    if mm.size not in (1, 3) or np.any(mm <= 0):
        raise ValueError("the thickness is one positive value in mm, or three (occipital, rest, frontal)")
    vx = grid_voxel_size(grid_affine)
    box = bounding_box(label > 0.5, int(np.ceil((mm.max() + 4) / vx.min())))
    frac0 = np.asarray(fractions[box], dtype=np.float32)
    lab = np.asarray(label[box], dtype=np.float32)

    # Regions: white matter loses its own, so that deep white matter takes
    # the region of the structure next to it
    regions = atlas.fill_from_nearest(regions_on_grid[box], vx, hole=_round(lab) == 3)
    keep = (regions > 0) & (regions < atlas.CORTEX_FROM) \
        & ~np.isin(regions, atlas.VENTRICLES + atlas.CEREBRAL_WM + (46,))
    keep = dilate(keep, 1.0, 1.0)                       # by one voxel
    # no band around the ventricles either: the CSF facing them is no sulcus
    no_band = keep | np.isin(regions, atlas.VENTRICLES)
    border_fwhm = 2.0                                   # width in mm of the blend into the kept labels
    keep_soft = np.clip(smooth(keep, border_fwhm, vx), 0, 1)
    if mm.size == 3:
        occipital, frontal = np.isin(regions, atlas.OCCIPITAL), np.isin(regions, atlas.FRONTAL)
        zones = [occipital, ~occipital & ~frontal, frontal]
    inside = _round(lab) > 0

    # A wider partial volume ramp, and thin white matter repaired
    lab = grey_close(smooth(lab, 1.25, vx), 2)
    # Depth of a voxel below the white matter surface in mm, to first order:
    # (label - threshold) / |gradient|.  It seeds the distance below, so the
    # zero of the distance lies on the surface and not on the voxel centres.
    gradient = np.sqrt(sum(g ** 2 for g in np.gradient(lab, *vx)))
    gradient = np.maximum(gradient, np.finfo(np.float32).eps)
    levels = np.linspace(0, 0.5 * vx.mean(), 5)

    out = np.zeros(lab.shape + (3,), dtype=np.float32)
    for step, offset in enumerate(PVE_OFFSETS):
        if verbose:
            print(f"  thickness: boundary offset {step + 1} of {PVE_OFFSETS.size}", flush=True)
        hard = _round(lab + offset)
        wm = largest_component(hard == 3)
        wm[keep] = False                                # kept structures must not seed the band
        depth = np.minimum((lab - (2.5 - offset)) / gradient, levels[-1])
        to_wm = np.full(lab.shape, np.inf, dtype=np.float32)
        for level in levels:
            seeds = wm & (depth >= level)
            if seeds.any():
                np.minimum(to_wm, distance(seeds, vx) - level, out=to_wm)
        bands = []
        for value in mm:
            band = hard.copy()
            band[~wm] = 1
            band[~inside] = 0
            band[(band == 1) & (to_wm <= value) & ~no_band] = 2
            band[keep] = hard[keep]
            bands.append(band)
        if mm.size == 1:
            result = bands[0]
        else:
            result = np.zeros_like(hard)
            for zone, band in zip(zones, bands):
                result[zone] = band[zone]
        for tissue in range(3):
            out[..., tissue] += (result == tissue + 1) / np.float32(PVE_OFFSETS.size)

    # the kept structures get their own fractions back, through a soft border
    out = out * (1 - keep_soft[..., None]) + frac0 * keep_soft[..., None]
    # So does the interior of the white matter: the hard labels made it
    # exactly 1 there, and its small deviations are what makes rendered white
    # matter look like tissue.  The blend stays 2 mm away from the boundary
    # that defines the thickness, and is capped so that the ground truth
    # cannot drop below white matter.
    core = np.clip(smooth(erode(out[..., _WM] > 0.99, 2.0, vx), border_fwhm, vx), 0, 1)
    core *= np.minimum(1, np.maximum(out[..., _WM] - 0.5, 0)
                       / np.maximum(out[..., _WM] - frac0[..., _WM], np.finfo(np.float32).eps))
    out = out * (1 - core[..., None]) + frac0 * core[..., None]

    new_label = out[..., 0] + 2 * out[..., 1] + 3 * out[..., 2]
    # holes that are left inside the brain are CSF
    brain = new_label > 0.5
    holes = brain != closing(brain, 4.0, 1.0)
    out[holes] = (1.0, 0.0, 0.0)
    new_label[holes] = 1.0

    all_fractions = np.array(fractions, dtype=np.float32)
    all_label = np.array(label, dtype=np.float32)
    all_fractions[box] = out
    all_label[box] = new_label
    if not return_rebuilt:
        return all_fractions, all_label
    # share of the original fractions in the two blends above
    original = 1 - (1 - keep_soft) * (1 - core)
    rebuilt = np.zeros(all_label.shape, dtype=bool)
    rebuilt[box] = (inside | holes) & ((original < 0.5) | holes)
    return all_fractions, all_label, rebuilt


# --------------------------------------------------------------------- atrophy
def atrophy_by_scaling(fractions: np.ndarray, roi: np.ndarray, factor: float) -> np.ndarray:
    """The atrophy of ``mri_simulate.m``: more CSF where there is CSF.

    The CSF fraction inside the region is multiplied by *factor* and the
    fractions are scaled back to their sum.  Kept for comparison; see the
    module notes for what it does and does not do.

    Args:
        fractions: (X, Y, Z, 3) fractions of CSF, GM and WM.
        roi: (X, Y, Z) boolean region.
        factor: Greater than 1.

    Returns:
        The new fractions.
    """
    out = np.array(fractions, dtype=np.float32)
    total = out.sum(-1)
    part = out[roi]
    part[:, _CSF] *= factor
    new_total = part.sum(-1, keepdims=True)
    part *= np.where(new_total > 0, total[roi][:, None] / np.maximum(new_total, 1e-12), 1)
    out[roi] = part
    return out


def _sample(volume: np.ndarray, coords: np.ndarray, nearest: bool = False) -> np.ndarray:
    """Samples of (C, X, Y, Z) *volume* at (3, ...) voxel coordinates: trilinear, or the nearest voxel."""
    import torch
    import torch.nn.functional as F

    vol = torch.from_numpy(np.ascontiguousarray(volume, dtype=np.float32))[None]
    size = torch.tensor(volume.shape[1:], dtype=torch.float32)
    grid = torch.from_numpy(np.ascontiguousarray(coords, dtype=np.float32))
    grid = grid.reshape(3, -1).T / (size - 1) * 2 - 1           # to -1..1, corners aligned
    grid = grid.flip(-1).reshape(1, *coords.shape[1:], 3) if coords.ndim == 4 \
        else grid.flip(-1).reshape(1, -1, 1, 1, 3)
    out = F.grid_sample(vol, grid, mode="nearest" if nearest else "bilinear", padding_mode="border",
                        align_corners=True)
    return out[0].reshape(volume.shape[0], *coords.shape[1:]).numpy()


def _identity(shape) -> np.ndarray:
    return np.stack(np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij"))


def _exponential(velocity: np.ndarray, steps: int = 7) -> np.ndarray:
    """Displacement of the flow of a stationary velocity field over unit time.

    Scaling and squaring: the field is divided by 2**steps, which makes it a
    small displacement, and composed with itself *steps* times.  Both are
    (3, X, Y, Z) in voxels.
    """
    grid = _identity(velocity.shape[1:])
    disp = velocity / np.float32(2 ** steps)
    for _ in range(steps):
        disp = disp + _sample(disp, grid + disp)
    return disp


def _jacobian(disp: np.ndarray) -> np.ndarray:
    """Jacobian determinant of x -> x + disp(x), (3, X, Y, Z) in voxels."""
    grad = np.stack([np.stack(np.gradient(disp[c]), 0) for c in range(3)], 0)    # d disp_c / d x_a
    for c in range(3):
        grad[c, c] += 1
    return (grad[0, 0] * (grad[1, 1] * grad[2, 2] - grad[1, 2] * grad[2, 1])
            - grad[0, 1] * (grad[1, 0] * grad[2, 2] - grad[1, 2] * grad[2, 0])
            + grad[0, 2] * (grad[1, 0] * grad[2, 1] - grad[1, 1] * grad[2, 0]))


def _velocity(divergence: np.ndarray, vx: np.ndarray, pad: int = 12) -> np.ndarray:
    """The curl-free velocity field with a given divergence, in voxels.

    Solves the Poisson equation in k-space.  The field is padded so that the
    periodic copies of the solution stay away.
    """
    div = np.pad(divergence.astype(np.float64), pad)
    div -= div.mean()
    spectrum = sp_fft.rfftn(div, workers=-1)
    k = [2 * np.pi * (sp_fft.rfftfreq(n, d=vx[a]) if a == 2 else sp_fft.fftfreq(n, d=vx[a]))
         for a, n in enumerate(div.shape)]
    kx, ky, kz = np.meshgrid(*k, indexing="ij")
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    k2[0, 0, 0] = 1.0
    potential = -spectrum / k2
    potential[0, 0, 0] = 0
    out = []
    for a, ka in enumerate((kx, ky, kz)):
        v = sp_fft.irfftn(1j * ka * potential, s=div.shape, workers=-1)          # mm per unit time
        out.append(v[pad:-pad, pad:-pad, pad:-pad] / vx[a])
    return np.stack(out).astype(np.float32)


def _volume_flow(fractions, roi, fine_vx, change, tissue_index, reach, resolution, verbose):
    """Backward and forward displacement and Jacobian on the fine grid for a change of volume.

    The tissue of the region is a sink and the CSF around it the source of a
    curl-free velocity field, found on a coarser grid; the flow shrinks the
    tissue the same in all directions.
    """
    shape = np.array(fractions.shape[:3])
    coarse = np.maximum(np.round(shape * fine_vx / resolution).astype(int), 8)
    vx = shape * fine_vx / coarse
    zoom = coarse / shape

    def to_coarse(volume):
        return ndimage.zoom(np.asarray(volume, dtype=np.float32), zoom, order=1, mode="nearest",
                            grid_mode=True)

    roi_c = np.clip(to_coarse(roi), 0, 1)
    sink = roi_c * sum(to_coarse(fractions[..., t]) for t in tissue_index)
    source = to_coarse(fractions[..., _CSF]) * np.clip(1 - distance(roi_c > 0.5, vx) / reach, 0, 1)
    if sink.sum() <= 0:
        raise ValueError("the region holds none of the tissues")
    if source.sum() <= 0:
        raise ValueError(f"there is no CSF within {reach} mm of the region to take up the volume")

    def flow(amplitude):
        # the CSF expands by as much as the tissue shrinks
        lost = np.sum(np.expm1(amplitude * sink))
        gain = brentq(lambda b: np.sum(np.expm1(b * source)) + lost, -50, 50)
        velocity = _velocity(smooth(amplitude * sink + gain * source, 1.5 * vx, vx), vx)
        jac = _jacobian(_exponential(velocity))
        return velocity, jac, float(np.sum(sink * (jac - 1)) / np.sum(sink))

    # The divergence stays where it was while the tissue moves, so the change
    # that is reached differs from the amplitude; adjust it
    amplitude = np.log1p(change)
    for iteration in range(4):
        velocity, jac, achieved = flow(amplitude)
        if verbose:
            print(f"  atrophy: amplitude {amplitude:+.4f} gives {100 * achieved:+.2f}% "
                  f"(aim {100 * change:+.2f}%), Jacobian {jac.min():.3f}..{jac.max():.3f}", flush=True)
        if abs(achieved - change) < 0.002 * max(abs(change), 0.05) or iteration == 3:
            break
        amplitude *= np.log1p(change) / np.log1p(achieved)

    # to the fine grid: positions and displacements both scale with the grids
    scale = zoom[:, None, None, None].astype(np.float32)
    at_coarse = (_identity(shape) + 0.5) * scale - 0.5
    backward = _sample(_exponential(-velocity), at_coarse) / scale
    forward = _sample(_exponential(velocity), at_coarse) / scale
    return backward, forward, _sample(jac[None], at_coarse)[0], achieved


def _thinning_flow(fractions, label, roi, fine_vx, change, tissue_index, reach, verbose):
    """Backward and forward displacement and Jacobian on the fine grid for a thinner cortex.

    Every point of the cortex moves along the normal of the white matter
    surface with a speed proportional to its distance from that surface, so
    those distances, and with them the thickness, are scaled while the white
    matter surface stays.  The CSF beyond the pial surface follows it and
    takes up the space within *reach*.  The speed is adjusted until the
    tissue of the region has changed by *change*: where two banks of a sulcus
    face each other closely, their opposite movements are smoothed into each
    other and the cortex thins less than the speed says.
    """
    shape = fractions.shape[:3]
    weight = np.clip(smooth(np.asarray(roi, dtype=np.float32), 3.0, fine_vx), 0, 1)
    box = bounding_box(weight > 0.01, int(np.ceil((reach + 4) / fine_vx.min())))
    lab = label[box]
    wm = lab >= 2.5
    if not wm.any():
        raise ValueError("there is no white matter near the region to measure the cortex from")
    depth = smooth(distance(wm, fine_vx), 1.0, fine_vx)          # mm above the WM surface
    normal = np.stack(np.gradient(depth, *fine_vx))
    normal /= np.maximum(np.sqrt((normal ** 2).sum(0)), 1e-6)
    tissue = lab >= 1.5
    index = ndimage.distance_transform_edt(~tissue, sampling=fine_vx, return_indices=True)
    beyond, nearest = index[0].astype(np.float32), tuple(index[1])
    # in the CSF the speed of the pial surface below, fading out with the distance from it
    speed = np.where(tissue, depth, depth[nearest] * np.exp(-beyond / (reach / 3)))
    speed[wm] = 0
    unit = weight[box] * speed * normal                           # mm per unit of the rate
    unit = np.stack([smooth(unit[a], 0.75, fine_vx) / fine_vx[a] for a in range(3)])
    sink = np.asarray(roi, dtype=np.float32)[box] * sum(fractions[box][..., t] for t in tissue_index)
    if sink.sum() <= 0:
        raise ValueError("the region holds none of the tissues")

    rate = np.log1p(change)
    for iteration in range(5):
        velocity = rate * unit
        forward = _exponential(velocity)
        jac = _jacobian(forward)
        achieved = float(np.sum(sink * (jac - 1), dtype=np.float64) / np.sum(sink, dtype=np.float64))
        if verbose:
            print(f"  atrophy: rate {rate:+.4f} gives {100 * achieved:+.2f}% "
                  f"(aim {100 * change:+.2f}%), Jacobian {jac.min():.3f}..{jac.max():.3f}", flush=True)
        if abs(achieved - change) < 0.01 * abs(change) or iteration == 4:
            break
        rate *= np.log1p(change) / np.log1p(achieved)
    # outside the box nothing moves
    full = (slice(None),) + box
    backward_fine = np.zeros((3,) + shape, dtype=np.float32)
    forward_fine = np.zeros((3,) + shape, dtype=np.float32)
    jac_fine = np.ones(shape, dtype=np.float32)
    backward_fine[full] = _exponential(-velocity)
    forward_fine[full] = forward
    jac_fine[box] = jac
    return backward_fine, forward_fine, jac_fine, achieved


def atrophy(fractions: np.ndarray, corrected: np.ndarray, grid_affine: np.ndarray,
            roi: np.ndarray, change: float, mode: str = "thinning", tissues=("GM",),
            reach: float = 10.0, resolution: float = 1.0, verbose: bool = False) -> dict:
    """Atrophy of a region by a deformation.

    The tissue of the region changes its volume by *change*, in one of two
    ways:

    - ``"thinning"``: every point of the cortex of the region moves towards
      the white matter surface in proportion to its distance from it, so the
      cortex keeps its area and gets thinner, and the CSF beyond it widens.
      This is what cortical atrophy mostly is: the loss of volume is a loss
      of thickness.
    - ``"volume"``: the tissue of the region is the sink and the CSF within
      *reach* the source of a curl-free flow.  The tissue shrinks the same
      in all directions, so a loss of 10% of volume is a loss of about 3% of
      thickness.

    In both the strength is adjusted until the tissue has changed by
    *change*.

    Args:
        fractions: (X, Y, Z, 3) fractions of CSF, GM and WM.
        corrected: (X, Y, Z) bias-corrected scan; it is deformed with the
            fractions, because it carries everything that is not brain.  The
            fractions are interpolated linearly, the image by a cubic
            B-spline.
        grid_affine: Voxel-to-world matrix of the grid.
        roi: (X, Y, Z) region, boolean or a weight in 0..1.
        change: Relative change, e.g. -0.1 for a loss of 10%.  Between -0.5
            and 0.5; a positive value thickens or enlarges.
        mode: ``"thinning"`` or ``"volume"``.
        tissues: The tissues whose volume is prescribed (``"volume"``) and
            whose change is reported: ``("GM",)`` or ``("GM", "WM")``.
        reach: Distance in mm within which the CSF takes up the space.
        resolution: ``"volume"`` only: voxel size in mm at which the
            deformation is computed before it is applied on the grid.
        verbose: Print the steps of the adjustment.

    Returns:
        A dict with the deformed ``fractions`` and ``corrected``, the
        ``label`` of the new fractions, and the deformation (see the module
        notes): ``jacobian`` ((X, Y, Z), on the anatomy that was passed in:
        below 1 where volume is lost), ``displacement`` and ``forward`` ((3,
        X, Y, Z) in voxels).  ``achieved`` is the relative volume change of
        the tissue of the region, followed with the tissue: the integral of
        the Jacobian over the tissue that was in the region.  What is found
        in the fixed region afterwards differs from it, because tissue moves
        across its border.
    """
    if not -0.5 <= change <= 0.5 or change == 0:
        raise ValueError("change is a relative change between -0.5 and 0.5, and not 0")
    if mode not in ("thinning", "volume"):
        raise ValueError("mode is 'thinning' or 'volume'")
    index = {"GM": _GM, "WM": _WM}
    try:
        tissue_index = [index[t] for t in tissues]
    except KeyError as exc:
        raise ValueError(f"tissues are 'GM' and 'WM', not {exc}") from None
    fractions = np.asarray(fractions, dtype=np.float32)
    if not np.any(roi):
        raise ValueError("the region is empty")
    fine_vx = grid_voxel_size(grid_affine)
    if mode == "volume":
        disp, forward, jac, achieved = _volume_flow(fractions, roi, fine_vx, change, tissue_index,
                                                    reach, resolution, verbose)
    else:
        label = fractions[..., 0] + 2 * fractions[..., 1] + 3 * fractions[..., 2]
        disp, forward, jac, achieved = _thinning_flow(fractions, label, roi, fine_vx, change,
                                                      tissue_index, reach, verbose)
    if jac.min() <= 0:
        raise RuntimeError("the deformation folds; choose a smaller change or a larger reach")

    moving = np.abs(disp).max(0) > 0.01
    if not moving.any():
        raise ValueError("the deformation moves nothing by more than 0.01 voxels: the region is too "
                         "small for this change")
    moved = bounding_box(moving, 2)
    stack = np.moveaxis(fractions, -1, 0)
    sub = (slice(None),) + moved
    grid = _identity([s.stop - s.start for s in moved])
    # The fractions by linear interpolation: they stay between 0 and 1 and
    # keep their sum, which a spline does not grant
    warped = stack.copy()
    warped[sub] = _sample(stack[sub], grid + disp[sub])
    new_fractions = np.ascontiguousarray(np.moveaxis(warped, 0, -1))
    # The image by a cubic B-spline, which blurs less; what it swings below 0 is set to 0
    new_corrected = np.array(corrected, dtype=np.float32)
    new_corrected[moved] = np.maximum(ndimage.map_coordinates(
        new_corrected[moved], grid + disp[sub], order=3, mode="nearest"), 0)
    return {
        "fractions": new_fractions,
        "label": new_fractions[..., 0] + 2 * new_fractions[..., 1] + 3 * new_fractions[..., 2],
        "corrected": new_corrected,
        "jacobian": jac,
        "displacement": disp,
        "forward": forward,
        "achieved": achieved,
    }


def warp(volume: np.ndarray, displacement: np.ndarray, nearest: bool = False) -> np.ndarray:
    """Read a volume at the places a displacement points to.

    With the ``displacement`` of :func:`atrophy` this deforms a volume as the
    anatomy was deformed; with its ``forward`` it brings a map of the new
    anatomy back onto the old one.

    Args:
        volume: (X, Y, Z), or (X, Y, Z, C) with a last axis of classes.
        displacement: (3, X, Y, Z) in voxels: the result at a voxel is the
            volume at that voxel plus its displacement.
        nearest: Take the nearest voxel instead of interpolating linearly:
            for an atlas, whose numbers must not be mixed.

    Returns:
        The volume in its shape: float32, or its own type with *nearest*.
    """
    dtype = np.asarray(volume).dtype if nearest else np.float32
    volume = np.asarray(volume, dtype=np.float32)
    classes = volume.ndim == 4
    stack = np.moveaxis(volume, -1, 0) if classes else volume[None]
    out = stack.copy()
    if np.abs(displacement).max() > 0.01:
        moved = bounding_box(np.abs(displacement).max(0) > 0.01, 2)
        sub = (slice(None),) + moved
        out[sub] = _sample(stack[sub], _identity([s.stop - s.start for s in moved]) + displacement[sub],
                           nearest=nearest)
    out = np.ascontiguousarray(np.moveaxis(out, 0, -1)) if classes else out[0]
    return out.astype(dtype, copy=False)


def compose(jacobian: np.ndarray, forward: np.ndarray, step: dict) -> tuple:
    """Jacobian determinant and forward displacement of two deformations, one after the other.

    The Jacobian of the second deformation lies on the anatomy the first one
    made.  It is read where the first one took every voxel, and only then
    multiplied: the plain product of the two maps is wrong wherever both
    change the volume.

    Args:
        jacobian: (X, Y, Z) Jacobian determinant of the first deformation.
        forward: (3, X, Y, Z) its forward displacement.
        step: The result of the :func:`atrophy` that followed.

    Returns:
        ``(jacobian, forward)`` of both together, on the anatomy before the
        first.
    """
    then = np.moveaxis(warp(np.moveaxis(step["forward"], 0, -1), forward), -1, 0)
    return jacobian * warp(step["jacobian"], forward), forward + then


def roi_from_atlas(regions_on_grid: np.ndarray, numbers) -> np.ndarray:
    """The voxels of one or several regions of an atlas on the grid."""
    numbers = np.atleast_1d(numbers)
    roi = np.isin(regions_on_grid, numbers)
    if not roi.any():
        raise ValueError(f"the atlas has none of the regions {numbers.tolist()} on the grid")
    return roi
