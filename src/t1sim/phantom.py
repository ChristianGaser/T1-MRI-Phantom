"""A folded sphere with an exactly known cortical thickness, and its evaluation.

The Python counterpart of ``matlab/thickness_phantom.m`` and
``matlab/thickness_phantom_eval.m``.  The phantom answers what a thickness
measure does when the geometry is known exactly: the white matter is a sphere
with regular folds, the grey matter is the band of constant thickness around
it, and a layer of CSF surrounds the grey matter.  There is no noise, no bias
field and no random number anywhere, so the error a thickness measure shows
on it is its own.

Why the thickness is exact:
    The grey matter is not made by offsetting the radius of the folded
    sphere.  A radial offset of t gives a normal thickness of t cos(alpha),
    with alpha the angle between the radius and the surface normal, which is
    too thin on the flanks of the folds.  Instead the Euclidean distance D to
    the white matter is computed once, and both cortical boundaries are level
    sets of that one distance map: the WM/GM boundary is D = wm_offset and the
    GM/CSF boundary is D = wm_offset + thickness.  D is 1-Lipschitz with a
    gradient of unit length, so two of its level sets are exactly
    ``thickness`` apart everywhere, whatever the folds look like.  The
    distance map is built on a supersampled grid, which keeps the
    discretisation of the surface well below the voxel size of the phantom.

Partial volume:
    The boundaries are shifted across a set of sub-voxel offsets, each offset
    gives a hard label image, and these are averaged.  The offsets are
    applied to the distance map itself, which moves all boundaries along
    their normal by a known amount.  They are the midpoints of ``pve_steps``
    equal intervals over ``pve_range`` voxels, so a locally flat boundary
    gets the exact linear ramp of one voxel width.

Evaluation:
    :func:`evaluate_phantom` checks two independent things, and the second
    cannot be read without the first.  (1) The fidelity of the phantom: the
    distance between the two boundaries is recovered from the label image
    alone, which is all a thickness pipeline ever sees, and the tissue
    volumes are compared with those of the supersampled grid.  (2) The
    accuracy of a thickness measure on that label image, CAT-Surface's
    projection-based thickness or any thickness map that is passed in.

Accuracy:
    Exact are the level sets of the distance map.  Two things are not, and
    for one reason: the distance to a voxelised object takes few values next
    to it (0, 1, sqrt(2), ... samples), and its level sets are the envelope
    of the spheres around the outermost samples, not their average.

    The sub-voxel offsets sweep the levels 0 to 2 wm_offset, the lowest of
    which the distance map does not resolve.  With the default wm_offset of
    1.5 fine samples the label of a sphere without folds holds a cortex that
    is 0.04 mm thicker than requested at 0.5 mm voxels and 0.07 mm at 1 mm,
    going by its tissue volumes; with a wm_offset of 6 fine samples it is
    0.01 and 0.02 mm.

    The readout of :func:`evaluate_phantom` takes half a sample off the
    distance to the nearest sample, which is right for a surface along the
    grid and too much for an oblique one.  On a label written down
    analytically, two concentric spheres with linear ramps, it reads 0.4
    samples of its grid too little: 0.07 mm with the defaults.  In the
    geometry bias of the default phantom (-0.018 mm) the two nearly cancel.

Precision:
    The distance is exact: scipy's Euclidean feature transform names the
    nearest white matter sample, and the squared distance to it is summed in
    integers.  It is kept in single precision, as CAT's ``cat_bwdist``
    returns it.  The analytic surface and the label arithmetic follow
    ``thickness_phantom.m`` operation by operation, each one rounded to
    single precision (:func:`_single`), because a sample within rounding of
    a boundary would otherwise fall on its other side than in MATLAB.  With
    that, label and T1w image of the default phantom and of 25 others were
    those of the MATLAB version bit by bit.

    Two further phantoms differed, and any will whose shifted levels coincide
    exactly with a distance the grid can take: e.g. 1 mm voxels with the
    default offsets, where 14 fine samples plus a third of a voxel are the
    5 mm of the outer CSF boundary.  Whether such a sample counts for that one
    sub-voxel offset is decided by the last bit of its distance, in MATLAB as
    here.  ``cat_bwdist`` sums the squared distance in single precision and
    is one float32 step off in 15 to 20% of the samples, the distance here is
    rounded correctly, so these voxels can differ by ``1 / pve_steps``.
"""

from __future__ import annotations

import json
import math
import os
import re
import warnings

import nibabel as nib
import numpy as np
from scipy import ndimage

from . import __version__

#: Order of the tissues in ``intensities`` and in the volumes, and their label.
LABELS = {"CSF": 1, "GM": 2, "WM": 3}

#: Folding patterns of the white matter surface.
FOLDS = ("cartesian", "spherical")


def _single(values) -> np.ndarray:
    """The result of one MATLAB operation that has a single precision operand.

    MATLAB computes such an operation in double precision and rounds the
    result to single; for two single operands this is the single precision
    operation itself, and with a double scalar (the radius, the wave number,
    a sub-voxel offset) it is not the same as rounding that scalar first.
    The values come back as float64, so that the next operation and the
    comparisons with the double thresholds need no further conversion.

    The elementary functions (cos, sin, acos, atan2) are rounded correctly
    this way.  MATLAB's single precision versions are not, and differ in the
    last bit for 2% (cos) to 16% (atan2) of the arguments.  That step cannot
    be reproduced; it matters only for a sample that lies within about
    1e-6 mm of the analytic surface.
    """
    return np.asarray(values, dtype=np.float64).astype(np.float32).astype(np.float64)


def _round(value: float) -> int:
    """MATLAB's ``round`` for a non-negative value: halves are rounded up."""
    return math.floor(value + 0.5)


def _bids_label(tag: str) -> str:
    """A valid BIDS entity label: 0.5 becomes 0p5, all else but letters and digits is dropped."""
    return re.sub(r"[^a-zA-Z0-9]", "", tag.replace(".", "p"))


def _distance_to(obj: np.ndarray, vx: float) -> np.ndarray:
    """Euclidean distance in mm to the nearest sample of ``obj``, 0 inside it.

    The quantity CAT's ``cat_bwdist(single(obj), vx*[1 1 1])`` returns: the
    distance between sample centres.  Only the feature transform is taken
    from scipy.  The distance scipy returns itself is float64 and is built
    from temporaries of the size of the feature transform, at a peak of 64
    bytes per sample or 3.5 GB for the default phantom.  Summing the squared
    offsets slab by slab needs 26 bytes per sample (1.5 GB), takes the same
    time, and gives the squared distance as an exact integer.

    Args:
        obj: (X, Y, Z) boolean array.
        vx: Isotropic sample distance in mm.

    Returns:
        (X, Y, Z) float32; infinite everywhere if ``obj`` is empty.
    """
    if not obj.any():
        return np.full(obj.shape, np.inf, dtype=np.float32)
    nearest = np.empty((3,) + obj.shape, dtype=np.int32)
    # scipy measures to the nearest zero, hence the inverted object
    ndimage.distance_transform_edt(
        ~obj, sampling=(vx,) * 3, return_distances=False, return_indices=True, indices=nearest)
    dist = np.empty(obj.shape, dtype=np.float32)
    j = np.arange(obj.shape[1], dtype=np.int32)[:, None]
    k = np.arange(obj.shape[2], dtype=np.int32)[None, :]
    for i in range(obj.shape[0]):
        squared = (nearest[0, i] - i) ** 2 + (nearest[1, i] - j) ** 2 + (nearest[2, i] - k) ** 2
        dist[i] = vx * np.sqrt(squared)
    return dist


def _count_within(dist: np.ndarray, level: float) -> int:
    """Number of samples of a single precision map that do not exceed ``level``.

    Compared in double precision, as MATLAB compares a single with a double:
    rounding the level to single first would move it by up to half a step.
    """
    return int(np.count_nonzero(np.less_equal(dist, level, signature="dd->?")))


def _folded_sphere(dimf, vxf, r0, amplitude, fold, wavelength, folds) -> np.ndarray:
    """The samples of the fine grid that lie inside the folded sphere.

    Args:
        dimf: Dimensions of the fine grid.
        vxf: Its sample distance in mm.
        r0: Mean radius of the analytic surface in mm.
        amplitude: Amplitude of the folds in mm.
        fold: ``"cartesian"`` or ``"spherical"``.
        wavelength: Fold wavelength in mm, for ``"cartesian"``.
        folds: Number of folds in theta and phi, for ``"spherical"``.

    Returns:
        Boolean array of shape ``dimf``.
    """
    # Coordinates in mm, centred on the volume.  The grid of the phantom has
    # the same origin, so each of its voxel centres coincides with a fine sample.
    x, y, z = (_single((np.arange(1, n + 1) - (n + 1) / 2) * vxf) for n in dimf)

    # The in-plane terms do not depend on the slice.  The volume is filled
    # slice by slice: coordinate volumes on the fine grid would need several
    # times the memory of the distance map.
    rxy2 = _single(_single(x * x)[:, None] + _single(y * y)[None, :])
    if fold == "cartesian":
        wave = 2 * np.pi / wavelength
        cxcy = _single(_single(np.cos(_single(wave * x)))[:, None]
                       * _single(np.cos(_single(wave * y)))[None, :])
        in_plane = _single(amplitude * cxcy)
    else:
        phi = _single(np.arctan2(y[None, :], x[:, None]))
        in_plane = _single(np.cos(_single(folds[1] * phi)))

    inside = np.empty(dimf, dtype=bool)
    tiny = float(np.finfo(np.float32).eps)
    for k in range(dimf[2]):
        r = _single(np.sqrt(_single(rxy2 + _single(z[k] * z[k]))))
        if fold == "cartesian":
            surface = _single(r0 + _single(in_plane * _single(np.cos(_single(wave * z[k])))))
        else:
            # theta is undefined at the centre, where r is 0.  Any value
            # does, the centre is deep inside the object.
            theta = _single(np.arccos(_single(z[k] / np.maximum(r, tiny))))
            along = _single(amplitude * _single(np.sin(_single(folds[0] * theta))))
            surface = _single(r0 + _single(along * in_plane))
        inside[:, :, k] = r <= surface
    return inside


def thickness_phantom(
    thickness=2.5,
    *,
    dim=(128, 128, 128),
    vx: float = 0.5,
    radius: float = 22.0,
    amplitude: float = 2.5,
    fold: str = "cartesian",
    wavelength: float = 12.0,
    folds=(6, 6),
    csf: float = 2.0,
    wm_offset: float | None = None,
    supersample: int = 3,
    pve_steps: int = 15,
    pve_range: float = 1.0,
    intensities=(30.0, 100.0, 150.0),
    outdir=None,
    name: str = "thicknessPhantom",
    verbose: bool = False,
):
    """Folded sphere phantom with an exactly known grey matter thickness.

    Args:
        thickness: GM thickness in mm.  A sequence gives one phantom per
            value; they share the distance map, by far the most expensive
            step, so a series is barely slower than a single phantom.
        dim: Dimensions of the phantom in voxels, one number or three.
        vx: Isotropic voxel size in mm.
        radius: Mean radius of the WM surface in mm.
        amplitude: Amplitude of the folds in mm, half their peak to peak
            height.
        fold: Folding pattern.  ``"cartesian"`` is
            ``amplitude * cos(kx) cos(ky) cos(kz)`` with
            ``k = 2 pi / wavelength``; the feature size is the same
            everywhere, which is what a phantom for a thickness measure
            needs.  ``"spherical"`` is
            ``amplitude * sin(n1 theta) cos(n2 phi)``, the pattern of
            ``internal/thickness_phantom.py`` of T1Prep.  Its folds get
            arbitrarily fine towards the poles, where no voxel size resolves
            them: use it to reproduce that phantom, not to measure a bias.
        wavelength: Fold wavelength in mm for ``"cartesian"``.
        folds: Number of folds ``(n1, n2)`` in theta and phi for
            ``"spherical"``.
        csf: Thickness of the CSF layer around the GM in mm.
        wm_offset: Level of the distance map that is the WM/GM boundary, in
            mm; default (None): half a voxel.  It has to stay clearly above
            the sampling of the distance map, since the level sets of a
            distance map are scalloped at that scale close to the object.
            The radius of the analytic surface is corrected for it, so
            ``radius`` stays the mean radius of the WM surface, to about a
            third of a sample of the distance map.
        supersample: Factor by which the grid is supersampled for the
            distance map.  An even factor is increased by one: an odd one
            puts a fine sample on every voxel centre of the phantom, and the
            distance is read there without interpolation.  Memory and time
            grow with its third power.
        pve_steps: Number of sub-voxel offsets that are averaged for the
            partial volume effect.
        pve_range: Total width of those offsets in voxels.  1 gives the
            exact linear ramp of a flat boundary.
        intensities: Intensity of CSF, GM and WM in the T1w image; the
            background is 0.
        outdir: Folder to write the images and their sidecars to with
            :func:`write_phantom`; default (None): nothing is written.
        name: Basename of those files.
        verbose: Print the progress.

    Returns:
        ``(label, t1w, info)``; three lists with one entry per value if
        ``thickness`` is a sequence.

        - ``label``: (X, Y, Z) float32, CSF = 1, GM = 2, WM = 3 with partial
          volumes in between, 0 outside;
        - ``t1w``: (X, Y, Z) float32, the tissue intensities mixed by the
          partial volume fractions;
        - ``info``: the geometry (``thickness``, ``vx``, ``dim``, ``radius``,
          ``amplitude``, ``fold``, ``wavelength``, ``folds``, ``csf``,
          ``wm_offset``, ``supersample``, ``pve_steps``, ``pve_range``,
          ``intensities``), the tissue volumes in ml of the supersampled
          grid (``volume_ref``) and of the label (``volume_pve``), which has
          to reproduce the former at a much coarser voxel size, and the
          names of the written ``files``.

    Raises:
        ValueError: A parameter is out of range, or the phantom does not fit
            into the field of view.
    """
    several = np.ndim(thickness) > 0
    values = [float(t) for t in np.atleast_1d(thickness)]
    if not values or not all(t > 0 for t in values):
        raise ValueError("thickness has to be positive")
    dim = [_round(n) for n in np.atleast_1d(dim)]
    if len(dim) == 1:
        dim = dim * 3
    if len(dim) != 3 or min(dim) < 1:
        raise ValueError(f"dim has to be one or three positive numbers, not {dim}")
    fold = str(fold).lower()
    if fold not in FOLDS:
        raise ValueError(f"unknown fold '{fold}', use 'cartesian' or 'spherical'")
    if fold == "spherical" and len(folds) != 2:
        raise ValueError("folds has to hold the number of folds in theta and in phi")
    if not vx > 0:
        raise ValueError("vx has to be positive")
    if fold == "cartesian" and not wavelength > 0:
        raise ValueError("wavelength has to be positive")
    if csf < 0:
        raise ValueError("csf must not be negative")
    if pve_steps != int(pve_steps) or pve_steps < 1:
        raise ValueError("pve_steps has to be a positive integer")
    if len(intensities) != 3:
        raise ValueError("intensities has to hold the values of CSF, GM and WM")
    pve_steps = int(pve_steps)

    ss = max(1, _round(supersample))
    if ss % 2 == 0:
        ss += 1
        warnings.warn(f"supersample has to be odd and was increased to {ss}", stacklevel=2)
    vxf = vx / ss
    dimf = [n * ss for n in dim]

    # Half a voxel is far enough from the object for the scalloping of the
    # level sets to stay below a tenth of a voxel.
    if wm_offset is None or not np.isfinite(wm_offset):
        wm_offset = 0.5 * vx
    wm_offset = float(wm_offset)
    if wm_offset <= vxf:
        warnings.warn(
            f"wm_offset of {wm_offset:g} mm is not larger than the sampling of the distance "
            f"map ({vxf:g} mm), so the WM surface will be visibly scalloped: increase "
            "wm_offset or supersample", stacklevel=2)

    # With the CSF layer cut open at the border the label image would no
    # longer describe the documented geometry.
    r_max = radius + amplitude + max(values) + csf
    if r_max + vx > min(dim) * vx / 2:
        raise ValueError(
            f"the phantom needs a radius of {r_max:g} mm but the field of view only offers "
            f"{min(dim) * vx / 2:g} mm: increase dim or reduce the geometry")

    # The object of the distance map is the set of samples inside the
    # analytic surface.  Its level set at wm_offset lies wm_offset outside of
    # it, and the voxelisation moves it back in by half a sample on average.
    # That average holds along a line.  In the volume the level set is the
    # envelope of the spheres around the outermost samples, and on a sphere
    # without folds it ends up 0.3 samples outside of radius.
    r0 = radius - wm_offset + 0.5 * vxf

    if verbose:
        print("Build folded sphere on a {}x{}x{} grid ({:.3f} mm).".format(*dimf, vxf))
    inside = _folded_sphere(dimf, vxf, r0, amplitude, fold, wavelength, folds)
    if not inside.any():
        raise ValueError("the folded sphere is empty, please check radius and amplitude")

    if verbose:
        print("Compute Euclidean distance map.")
    dist = _distance_to(inside, vxf)
    del inside

    # Tissue volumes of the fine grid in ml, the reference for the label.
    # The three tissues are nested, so each is counted as what a level adds
    # to the one below.
    n_wm = _count_within(dist, wm_offset)
    volume_ref = []
    for t in values:
        n_gm = _count_within(dist, wm_offset + t) - n_wm
        n_csf = _count_within(dist, wm_offset + t + csf) - n_gm - n_wm
        volume_ref.append({key: n * vxf ** 3 / 1000
                           for key, n in (("csf", n_csf), ("gm", n_gm), ("wm", n_wm))})

    # The distance on the grid of the phantom: the fine samples that sit on
    # the voxel centres
    half = (ss - 1) // 2
    dist = dist[half::ss, half::ss, half::ss].astype(np.float64)

    # Sub-voxel offsets of the boundary: the midpoints of pve_steps equal
    # intervals.  The midpoint rule integrates the uniform distribution of
    # the boundary position exactly, while equally spaced offsets that
    # include both ends would weight the two extremes like all others and
    # bias the ramp.
    offsets = pve_range * vx * ((np.arange(1, pve_steps + 1) - 0.5) / pve_steps - 0.5)
    if offsets[0] < -wm_offset:
        # the distance is 0 in the whole object, there is nothing left to shift
        warnings.warn(
            f"the sub-voxel offsets reach {-offsets[0]:g} mm inwards, further than wm_offset "
            f"({wm_offset:g} mm), so no voxel will be pure WM: reduce pve_range or increase "
            "wm_offset", stacklevel=2)

    # Number of offsets at which a voxel lies within the outer boundary of
    # the WM, of the GM and of the CSF.  The shifted distance and the WM do
    # not depend on the thickness, which is why the offsets are the outer loop.
    if verbose:
        print(f"Average {pve_steps} sub-voxel offsets for {len(values)} thickness value(s).")
    in_wm = np.zeros(dim, dtype=np.int32)
    in_gm = [np.zeros(dim, dtype=np.int32) for _ in values]
    in_csf = [np.zeros(dim, dtype=np.int32) for _ in values]
    for offset in offsets:
        # Subtracting the offset moves all three boundaries outwards by that
        # amount along their normal
        d = _single(dist - offset)
        in_wm += d <= wm_offset
        for i, t in enumerate(values):
            in_gm[i] += d <= wm_offset + t
            in_csf[i] += d <= wm_offset + t + csf
    del dist

    labels, images, infos = [], [], []
    for i, (t, ref) in enumerate(zip(values, volume_ref)):
        f_wm = _single(in_wm / pve_steps)
        f_gm = _single((in_gm[i] - in_wm) / pve_steps)
        f_csf = _single((in_csf[i] - in_gm[i]) / pve_steps)

        label = _single(_single(f_csf + _single(2 * f_gm)) + _single(3 * f_wm))
        # The intensity is the mixture of the tissue intensities weighted by
        # the partial volume fractions, the background contributes 0
        t1w = _single(_single(_single(intensities[0] * f_csf) + _single(intensities[1] * f_gm))
                      + _single(intensities[2] * f_wm))

        info = {
            "thickness": t,
            "vx": float(vx),
            "dim": list(dim),
            "radius": float(radius),
            "amplitude": float(amplitude),
            "fold": fold,
            "wavelength": float(wavelength),
            "folds": [int(n) if n == int(n) else float(n) for n in folds],
            "csf": float(csf),
            "wm_offset": wm_offset,
            "supersample": ss,
            "pve_steps": pve_steps,
            "pve_range": float(pve_range),
            "intensities": [float(v) for v in intensities],
            "volume_ref": ref,
            "volume_pve": {key: float(f.sum()) * vx ** 3 / 1000
                           for key, f in (("csf", f_csf), ("gm", f_gm), ("wm", f_wm))},
            "files": {"label": None, "t1w": None},
        }
        label = label.astype(np.float32)
        t1w = t1w.astype(np.float32)
        if outdir is not None:
            info["files"] = write_phantom(label, t1w, info, outdir, name=name, verbose=verbose)
        labels.append(label)
        images.append(t1w)
        infos.append(info)

    if several:
        return labels, images, infos
    return labels[0], images[0], infos[0]


def write_phantom(label, t1w, info: dict, outdir=".", name: str = "thicknessPhantom",
                  verbose: bool = False) -> dict:
    """Write the label image, the T1w image and their JSON sidecars.

    The names follow the BIDS filename grammar in the way ``mri_simulate.m``
    uses it: entity-value pairs and a suffix, the voxel size in the ``res``
    entity and the thickness in an alphanumeric ``desc`` label::

        <name>_res-<vx>mm_desc-Thickness<10 t>mm_dseg.nii.gz   label image
        <name>_res-<vx>mm_desc-Thickness<10 t>mm_T1w.nii.gz    ideal T1w image

    with a ``.json`` next to each.  The sidecars record the geometry: the
    whole point of the phantom is that its thickness is known, and
    :func:`evaluate_phantom` reads it back from there.  The world origin is
    the centre of the volume, which is the centre of the sphere.

    Args:
        label: (X, Y, Z) label image of :func:`thickness_phantom`.
        t1w: (X, Y, Z) T1w image that belongs to it.
        info: The ``info`` of that phantom.
        outdir: Output folder, created if it is missing.
        name: Basename of the files.
        verbose: Print the file names.

    Returns:
        ``{"label": ..., "t1w": ...}``, the names of the two images.
    """
    label = np.asarray(label, dtype=np.float32)
    t1w = np.asarray(t1w, dtype=np.float32)
    if label.ndim != 3 or label.shape != t1w.shape:
        raise ValueError(f"label and t1w have to be 3-D volumes of one shape, "
                         f"not {label.shape} and {t1w.shape}")
    outdir = os.fspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    vx, t = info["vx"], info["thickness"]
    base = "{}_res-{}mm_desc-{}".format(
        name, _bids_label(f"{vx:g}"), _bids_label(f"Thickness{_round(10 * t):02d}mm"))
    files = {"label": os.path.join(outdir, base + "_dseg.nii.gz"),
             "t1w": os.path.join(outdir, base + "_T1w.nii.gz")}

    affine = np.diag([vx, vx, vx, 1.0])
    affine[:3, 3] = -vx * (np.array(label.shape) - 1) / 2
    for key, data in (("label", label), ("t1w", t1w)):
        if verbose:
            print(f"Save {files[key]}")
        image = nib.Nifti1Image(data, affine)
        image.header.set_data_dtype(np.float32)
        image.header.set_xyzt_units("mm", "sec")
        image.header["descrip"] = b"thickness_phantom"
        image.set_sform(affine, code="aligned")
        image.set_qform(affine, code="aligned")
        nib.save(image, files[key])

    generated = [{"Name": "t1sim", "Version": __version__}]
    parameters = {
        "Thickness": t,
        "CSFLayer": info["csf"],
        "Radius": info["radius"],
        "Amplitude": info["amplitude"],
        "Fold": info["fold"],
    }
    if info["fold"] == "cartesian":
        parameters["FoldWavelength"] = info["wavelength"]
    else:
        parameters["FoldNumber"] = list(info["folds"])
    parameters.update({
        "WMOffset": info["wm_offset"],
        "VoxelSize": [vx, vx, vx],
        "Supersample": info["supersample"],
        "PVESteps": info["pve_steps"],
        "PVERange": info["pve_range"],
    })
    sidecars = {
        "label": {
            "GeneratedBy": generated,
            "Manual": False,
            "Description": (
                "Ground truth partial volume label image of a folded sphere with a constant "
                "GM thickness. Values are continuous and interpolate between the tissue "
                "labels, e.g. 2.5 is an equal mixture of GM and WM."),
            "Labels": dict(LABELS),
            "PhantomParameters": parameters,
        },
        "t1w": {
            "GeneratedBy": generated,
            "Sources": ["bids::" + os.path.basename(files["label"])],
            "Description": (
                "Ideal T1w image of the phantom, i.e. the tissue intensities mixed by the "
                "partial volume fractions. It carries no noise, no bias field and no other "
                "artefact, use mri_simulate for those."),
            "PhantomParameters": dict(parameters, Intensities=list(info["intensities"])),
        },
    }
    for key, meta in sidecars.items():
        with open(_sidecar(files[key]), "w", encoding="utf-8") as handle:
            json.dump(meta, handle, indent=2)
            handle.write("\n")
    return files


def _sidecar(fname: str) -> str:
    """Name of the JSON sidecar of an image."""
    return re.sub(r"\.nii(\.gz)?$", "", fname) + ".json"


def _along(axis: int, *item) -> tuple:
    """Index of a volume: ``slice(*item)`` along one axis and all of the other two."""
    index = [slice(None)] * 3
    index[axis] = slice(*item)
    return tuple(index)


def _upsample(label: np.ndarray, ss: int) -> np.ndarray:
    """Trilinear interpolation of a volume on a grid that is ``ss`` times finer.

    The fine grid is centred as in :func:`thickness_phantom`: with an odd
    ``ss`` one fine sample sits on every voxel centre.  The interpolation
    adds no geometry, it only lets a distance map resolve the boundaries
    better than the voxel size would allow.  The fine samples beyond the
    outermost voxel centres are 0, i.e. background; the phantom keeps a
    margin to the border, so nothing is lost there.

    Args:
        label: (X, Y, Z) volume.
        ss: Odd supersampling factor.

    Returns:
        (ss X, ss Y, ss Z) float32.
    """
    fine = np.asarray(label, dtype=np.float32)
    half = (ss - 1) // 2
    for axis in range(3):
        n = fine.shape[axis]
        shape = list(fine.shape)
        shape[axis] = n * ss
        out = np.zeros(shape, dtype=np.float32)
        out[_along(axis, half, None, ss)] = fine
        for step in range(1, ss):
            # the fine samples that lie step/ss of the way from a voxel to the next
            weight = np.float32(step / ss)
            out[_along(axis, half + step, half + step + ss * (n - 1), ss)] = (
                (1 - weight) * fine[_along(axis, 0, n - 1)] + weight * fine[_along(axis, 1, n)])
        fine = out
    return fine


def _sample_at_isolevel(fine: np.ndarray, level: float, values: np.ndarray) -> np.ndarray:
    """Read a map on an isosurface of a second map with sub-voxel accuracy.

    Each pair of neighbouring samples along the three axes whose values
    straddle the isolevel holds one crossing.  Its position between the two
    samples follows from linear interpolation, and ``values`` is interpolated
    there as well.

    A layer of voxels around the isosurface would not do: it has a thickness
    of its own, its centres sit a fraction of a sample inside the surface,
    and that fraction depends on the orientation of the surface, which
    biases what is read on it.

    Args:
        fine: (X, Y, Z) map whose isosurface says where to sample.
        level: The isolevel.
        values: (X, Y, Z) map to read.

    Returns:
        float64 vector with one entry per crossing.
    """
    above = fine >= level
    samples = []
    for axis in range(3):
        first, second = _along(axis, 0, -1), _along(axis, 1, None)
        crossing = above[first] != above[second]
        a = fine[first][crossing].astype(np.float64)
        b = fine[second][crossing].astype(np.float64)
        v1 = values[first][crossing].astype(np.float64)
        v2 = values[second][crossing].astype(np.float64)
        # the two values straddle the level, so they differ and the weight lies in [0, 1)
        samples.append(v1 + (level - a) / (b - a) * (v2 - v1))
    return np.concatenate(samples)


def _stats(values, thickness: float) -> dict:
    """Summary of measured thickness values against the known thickness.

    Args:
        values: Measured values in mm; those that are not finite are dropped.
        thickness: Known thickness in mm.

    Returns:
        ``n``, ``mean``, ``median``, ``sd``, ``bias`` (mean - thickness),
        ``rmse``, ``p5`` and ``p95``; NaN if nothing was measured.
    """
    values = np.asarray(values, dtype=np.float64).ravel()
    values = np.sort(values[np.isfinite(values)])
    n = values.size
    if n == 0:
        nan = float("nan")
        return {"n": 0, "mean": nan, "median": nan, "sd": nan, "bias": nan, "rmse": nan,
                "p5": nan, "p95": nan}

    def quantile(q):
        # the sample of that rank, as thickness_phantom_eval.m takes it
        return float(values[max(1, _round(q * n)) - 1])

    mean = float(values.mean())
    return {
        "n": n,
        "mean": mean,
        "median": quantile(0.5),
        "sd": float(values.std(ddof=1)) if n > 1 else 0.0,
        "bias": mean - thickness,
        "rmse": float(np.sqrt(np.mean((values - thickness) ** 2))),
        "p5": quantile(0.05),
        "p95": quantile(0.95),
    }


def _each(value, n: int, what: str, shared: bool = False) -> list:
    """One entry of an argument per phantom."""
    if value is None:
        return [None] * n
    if isinstance(value, (list, tuple)) or (shared and np.ndim(value) > 0):
        if len(value) != n:
            raise ValueError(f"{what} has {len(value)} entries but {n} phantoms were given")
        return list(value)
    if n > 1 and not shared:
        raise ValueError(f"{what} has to be a sequence with one entry per phantom")
    return [value] * n


def _read_phantom(source, info, thickness, vx, name: str) -> dict:
    """The label image of one phantom with its voxel size and known thickness.

    An explicit ``thickness`` or ``vx`` wins over ``info``, which wins over
    the sidecar and the header of a file.
    """
    info = info or {}
    job = {"name": name, "thickness": thickness, "vx": vx, "volume_ref": info.get("volume_ref")}
    if job["thickness"] is None:
        job["thickness"] = info.get("thickness")
    if job["vx"] is None:
        job["vx"] = info.get("vx")

    if isinstance(source, (str, os.PathLike)):
        fname = os.fspath(source)
        image = nib.load(fname)
        job["name"] = os.path.basename(fname)
        job["label"] = np.squeeze(np.asarray(image.dataobj, dtype=np.float32))
        sizes = np.sqrt((image.affine[:3, :3] ** 2).sum(axis=0))
        if job["vx"] is None:
            job["vx"] = float(sizes.mean())
            if np.any(np.abs(sizes - job["vx"]) > 1e-6):
                warnings.warn("{} is not isotropic ({:g} {:g} {:g} mm), the mean is used".format(
                    job["name"], *sizes), stacklevel=3)
        if job["thickness"] is None and os.path.isfile(_sidecar(fname)):
            with open(_sidecar(fname), encoding="utf-8") as handle:
                job["thickness"] = json.load(handle).get("PhantomParameters", {}).get("Thickness")
    else:
        job["label"] = np.asarray(source, dtype=np.float32)

    if job["label"].ndim != 3:
        raise ValueError(f"{job['name']} must be a single 3-D volume, not {job['label'].shape}")
    if job["vx"] is None:
        raise ValueError("vx is required when the label image is passed directly")
    if job["thickness"] is None or not np.isfinite(job["thickness"]):
        raise ValueError(
            f"the thickness of {job['name']} is unknown: pass thickness or the info of "
            "thickness_phantom, or keep the JSON sidecar next to the label image")
    job["thickness"] = float(job["thickness"])
    job["vx"] = float(job["vx"])
    return job


def _pbt(label: np.ndarray, vx: float, options: dict | None):
    """CAT-Surface's projection-based thickness, or None if it is not installed."""
    try:
        import cat_surf
    except ImportError:
        warnings.warn("cat_surf was not found, the thickness measure is skipped", stacklevel=3)
        return None
    gmt, ppm = cat_surf.vol_thickness_pbt(
        np.ascontiguousarray(label, dtype=np.float32), voxelsize=[vx] * 3, **(options or {}))[:2]
    return np.asarray(gmt), np.asarray(ppm)


def evaluate_phantom(
    label,
    info=None,
    *,
    thickness=None,
    vx: float | None = None,
    supersample: int = 3,
    pbt: bool = True,
    pbt_options: dict | None = None,
    thickness_map=None,
    position_map=None,
    fig=None,
    verbose: bool = False,
):
    """Fidelity of a thickness phantom, and accuracy of a thickness measure on it.

    1. Fidelity of the phantom.  The label image is all a thickness pipeline
       ever sees, so the question is not whether the requested geometry is
       right but whether the label image still carries it at its voxel size.
       The label is supersampled by trilinear interpolation and thresholded
       at 2.5 and at 1.5, the two boundaries of a partial volume label.  The
       Euclidean distance from the WM surface is read on the GM/CSF surface
       at the crossings of the isolevel along the three axes
       (``geom["outer"]``).  This models nothing, so a deviation is an error
       of the phantom and not of a thickness method.

       The opposite direction (``geom["inner"]``) is not redundant.  The
       distance from the WM outwards is the thickness by construction, but
       the distance from the GM/CSF surface back to the WM is only the same
       while there is CSF on the other side of the sulcus.  Where a fold is
       narrower than twice the thickness the two banks of GM touch, the CSF
       is squeezed out, and this direction measures how far away the nearest
       remaining CSF is.  That is the situation in which a thickness measure
       overestimates, so the share of the WM surface that faces such a
       buried sulcus is reported (``geom["buried"]``).

       The tissue volumes of the label are compared with those of the
       supersampled grid the phantom was built on.  They test the partial
       volume label as a whole, including the voxels no boundary passes
       through.

    2. Accuracy of a thickness measure: CAT-Surface's projection-based
       thickness (``cat_surf.vol_thickness_pbt``), which takes exactly this
       kind of label image, or a thickness map that is passed in.  It is
       compared with the known thickness over all GM voxels and, where the
       percentage position is known, on the central surface.

    Args:
        label: The label image: an (X, Y, Z) array, the name of a ``_dseg``
            file of :func:`write_phantom`, or a list of either for a series
            of phantoms.
        info: The ``info`` of :func:`thickness_phantom` (a list for a
            series).  It supplies the thickness, the voxel size and the
            reference volumes.
        thickness: Known thickness in mm (a sequence for a series).  It wins
            over that of ``info``, which wins over the sidecar of a file.
        vx: Isotropic voxel size in mm.  It wins over that of ``info``, which
            wins over the header of a file.
        supersample: Factor of the geometric measurement.  The boundaries
            are recovered from a thresholded image, so that measurement is
            accurate to about half a sample of this grid; on an analytic
            sphere it reads 0.4 samples too little (see the module
            docstring).  An even factor is increased by one.
        pbt: Run ``cat_surf.vol_thickness_pbt`` if no ``thickness_map`` is
            given.  Skipped with a warning if cat_surf is not installed.
        pbt_options: Keyword arguments passed on to it; default: its own
            defaults.  One of them adds 0.25 mm to the thickness for the
            border shift of a real segmentation, which the phantom does not
            have: ``{"correct_thickness": 0.0}`` shows the measure itself.
        thickness_map: (X, Y, Z) thickness map in mm of any other measure on
            the grid of the label, 0 where it is undefined (a list for a
            series).
        position_map: Percentage position map that belongs to
            ``thickness_map``, 0.5 on the central surface (a list for a
            series).
        fig: Name of a QC figure (PNG) with a central slice of each phantom
            and the histograms of the thickness; needs matplotlib.
        verbose: Print the report.

    Returns:
        A dict, or a list of them if ``label`` is a list:

        - ``name``, ``thickness`` and ``vx`` of the phantom;
        - ``geom``: ``outer`` and ``inner``, the statistics of the two
          distances, and ``buried``, the share of the WM surface whose
          distance to the CSF exceeds the thickness by more than half a
          voxel;
        - ``volume``: tissue volumes in ml of the ``label`` and of the
          reference grid (``ref``, None without ``info``);
        - ``measure``: None, or the ``name`` of the measure and its
          statistics over all ``gm`` voxels and on the ``central`` surface
          (None without a position map).

        All statistics hold ``n``, ``mean``, ``median``, ``sd``, ``bias``
        (mean - thickness), ``rmse``, ``p5`` and ``p95``.
    """
    several = isinstance(label, (list, tuple))
    sources = list(label) if several else [label]
    n = len(sources)
    infos = _each(info, n, "info")
    known = _each(thickness, n, "thickness", shared=True)
    maps = _each(thickness_map, n, "thickness_map")
    positions = _each(position_map, n, "position_map")

    ss = max(1, _round(supersample))
    if ss % 2 == 0:
        ss += 1
    if fig is not None:
        # a missing matplotlib should stop the evaluation before it starts, not after
        import matplotlib  # noqa: F401

    jobs, results = [], []
    for i, source in enumerate(sources):
        job = _read_phantom(source, infos[i], known[i], vx,
                            f"label image {i + 1}" if several else "label image")
        y, size, t = job["label"], job["vx"], job["thickness"]
        gmt, ppm, measured_by = maps[i], positions[i], "thickness map"
        if ppm is not None and gmt is None:
            raise ValueError("a position_map needs the thickness_map it belongs to")
        if any(given is not None and np.shape(given) != y.shape for given in (gmt, ppm)):
            raise ValueError(f"the thickness and position maps of {job['name']} must have "
                             f"the shape of the label image, {y.shape}")
        if verbose:
            print("\n{0}\n{1} (thickness {2:.2f} mm, {3:.2f} mm voxels)\n{0}".format(
                "=" * 80, job["name"], t, size))
            print(f"Recover the boundaries at {size / ss:.3f} mm and measure their distance.")

        # The distance map measures between sample centres.  Half a sample
        # less puts the surface on the face of the last object voxel instead
        # of on its centre, which is exact for a surface along the grid only.
        fine = _upsample(y, ss)
        face = 0.5 * size / ss
        outer = _sample_at_isolevel(fine, 1.5, _distance_to(fine >= 2.5, size / ss)) - face
        inner = _sample_at_isolevel(fine, 2.5, _distance_to(fine < 1.5, size / ss)) - face
        del fine
        buried = float(np.mean(inner > t + 0.5 * size)) if inner.size else float("nan")

        # The label mixes neighbouring tissues linearly, so the fraction of
        # a tissue is the triangular weight around its label
        volume = {}
        for tissue, value in LABELS.items():
            fraction = np.clip(1 - np.abs(y - value), 0, None)
            volume[tissue.lower()] = float(fraction.sum(dtype=np.float64)) * size ** 3 / 1000

        if gmt is None and pbt:
            if verbose:
                print("Run cat_surf.vol_thickness_pbt.")
            estimate = _pbt(y, size, pbt_options)
            if estimate is not None:
                gmt, ppm = estimate
                measured_by = "cat_surf.vol_thickness_pbt"
        measure = None
        if gmt is not None:
            gmt = np.asarray(gmt)
            measure = {"name": measured_by,
                       "gm": _stats(gmt[(gmt > 0) & (y > 1.5) & (y < 2.5)], t),
                       "central": None}
            if ppm is not None:
                ppm = np.asarray(ppm)
                # the central surface is where CAT samples the thickness for
                # the surface based analysis
                measure["central"] = _stats(gmt[(gmt > 0) & (ppm > 0.45) & (ppm < 0.55)], t)

        results.append({
            "name": job["name"],
            "thickness": t,
            "vx": size,
            "geom": {"outer": _stats(outer, t), "inner": _stats(inner, t), "buried": buried},
            "volume": {"label": volume, "ref": job["volume_ref"]},
            "measure": measure,
        })
        if verbose:
            print(_report(results[-1]))
        if fig is not None:
            # the samples are far too many to return them, they are only kept
            # for the histogram of the figure
            jobs.append({"label": y, "map": gmt, "outer": outer})

    if fig is not None:
        _figure(jobs, results, os.fspath(fig))
        if verbose:
            print(f"\nSave {os.fspath(fig)}")
    return results if several else results[0]


def _report(result: dict) -> str:
    """The result of one phantom as a table."""
    lines = ["", f"  Known thickness {result['thickness']:.3f} mm", "",
             "  {:<34} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
                 "measure", "n", "mean", "sd", "bias", "rmse", "p5-p95"),
             "  " + "-" * 78]

    def row(name, s):
        lines.append("  {:<34} {:6d} {:8.3f} {:8.3f} {:+8.3f} {:8.3f} {:6.2f}-{:.2f}".format(
            name, s["n"], s["mean"], s["sd"], s["bias"], s["rmse"], s["p5"], s["p95"]))

    row("geometry: WM surf -> GM/CSF surf", result["geom"]["outer"])
    row("geometry: GM/CSF surf -> WM surf", result["geom"]["inner"])
    measure = result["measure"]
    if measure is not None:
        short = measure["name"].split(".")[-1]
        row(f"{short}, all GM", measure["gm"])
        if measure["central"] is not None:
            row(f"{short}, central surf", measure["central"])
    lines += ["", "  {:.1f}% of the WM surface faces a sulcus without CSF (buried sulcus).".format(
        100 * result["geom"]["buried"])]

    volume, ref = result["volume"]["label"], result["volume"]["ref"]
    tissues = ("csf", "gm", "wm")
    lines += ["", "  {:<10} {:>10} {:>10} {:>10}".format("volume/ml", "CSF", "GM", "WM"),
              "  {:<10} {:10.2f} {:10.2f} {:10.2f}".format("label", *(volume[k] for k in tissues))]
    if ref is not None:
        lines += ["  {:<10} {:10.2f} {:10.2f} {:10.2f}".format(
                      "reference", *(ref[k] for k in tissues)),
                  "  {:<10} {:9.2f}% {:9.2f}% {:9.2f}%".format(
                      "error", *(100 * (volume[k] / ref[k] - 1) for k in tissues))]
    return "\n".join(lines)


def _figure(jobs: list, results: list, fname: str) -> None:
    """QC figure: per phantom the label, the thickness map and the histograms.

    A central slice of the label image, of the thickness map where there is
    one, and the histogram of the measured thickness against the known one.
    Drawn on an Agg canvas of its own, so the backend of the caller's
    session is left alone.

    Args:
        jobs: Per phantom the ``label``, the thickness ``map`` (or None) and
            the ``outer`` distances.
        results: The results of :func:`evaluate_phantom`.
        fname: Name of the PNG to write.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    # Crop to the bounding box of the largest phantom, so the slices of a
    # series stay comparable and no space is spent on background
    slices = [job["label"][:, :, _round(job["label"].shape[2] / 2) - 1] for job in jobs]
    boxes = []
    for image in slices:
        rows, cols = np.nonzero(image > 0)
        boxes.append((rows.min(), rows.max(), cols.min(), cols.max()) if rows.size
                     else (0, image.shape[0] - 1, 0, image.shape[1] - 1))
    if len({image.shape for image in slices}) == 1:
        boxes = [(min(b[0] for b in boxes), max(b[1] for b in boxes),
                  min(b[2] for b in boxes), max(b[3] for b in boxes))] * len(jobs)

    def crop(image, box):
        return np.rot90(image[max(0, box[0] - 2):box[1] + 3, max(0, box[2] - 2):box[3] + 3])

    n = len(jobs)
    figure = Figure(figsize=(10, 2.6 * n), dpi=110, facecolor="w")
    FigureCanvasAgg(figure)
    axes = figure.subplots(n, 3, squeeze=False)
    for job, result, image, box, (left, middle, right) in zip(jobs, results, slices, boxes, axes):
        t = result["thickness"]
        measure = result["measure"]

        left.imshow(crop(image, box), cmap="gray", vmin=0, vmax=3)
        left.set_title(f"label, known {t:.2f} mm", fontsize=9)
        left.set_axis_off()

        middle.set_axis_off()
        if job["map"] is None:
            middle.set_title("no thickness map", fontsize=9)
        else:
            # Only where the thickness was estimated, inside the GM band.
            # The background and the WM would otherwise dominate the colour
            # scale and hide the variation that matters.
            estimate = crop(job["map"][:, :, _round(job["map"].shape[2] / 2) - 1], box)
            band = crop(image, box)
            shown = np.ma.masked_where(~((band > 1.25) & (band < 2.75) & (estimate > 0)), estimate)
            mappable = middle.imshow(shown, cmap="jet", vmin=t - 0.4, vmax=t + 0.4)
            figure.colorbar(mappable, ax=middle)
            middle.set_title(measure["name"], fontsize=9)

        edges = t + np.linspace(-0.8, 0.8, 80)
        centres = (edges[:-1] + edges[1:]) / 2
        right.plot(centres, np.histogram(job["outer"], edges)[0] / max(1, job["outer"].size),
                   "k", linewidth=1, label="geometry")
        if job["map"] is not None:
            values = job["map"][(job["map"] > 0) & (job["label"] > 1.5) & (job["label"] < 2.5)]
            right.plot(centres, np.histogram(values, edges)[0] / max(1, values.size),
                       "b", linewidth=1, label=measure["name"].split(".")[-1])
        right.axvline(t, color="r", linestyle="--", label="known")
        right.set_xlim(edges[0], edges[-1])
        right.set_xlabel("thickness (mm)")
        right.set_ylabel("fraction")
        right.legend(loc="upper left", fontsize=7, frameon=False)
        title = f"geometry bias {result['geom']['outer']['bias']:+.3f} mm"
        if measure is not None:
            title = "geometry {:+.3f}, measure {:+.3f} mm".format(
                result["geom"]["outer"]["bias"], measure["gm"]["bias"])
        right.set_title(title, fontsize=9)

    figure.tight_layout()
    figure.savefig(fname, dpi=110)
