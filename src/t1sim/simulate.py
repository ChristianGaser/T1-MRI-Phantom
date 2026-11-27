"""The simulation: from a real scan to simulated images with their ground truth.

Three steps, which can be called on their own:

1. :func:`ground_truth` runs the front end once and makes the anatomy: the
   tissue fractions the images are rendered from, an optional change of the
   anatomy (cortical thickness, atrophy), the cleaned ground truth, and the
   image without artefacts.  This is the expensive step, about two minutes.
2. :func:`acquire` makes one simulated image of that anatomy: contrast, bias
   field, voxel size, movement, ringing and noise.  It takes seconds, and any
   number of acquisitions can share one ground truth.
3. :func:`write` stores an acquisition and its ground truth with names and
   sidecars after BIDS (:mod:`t1sim.bids`).

:func:`simulate` does all three for one scan.

**Two sets of fractions.**  The image is rendered from the fractions the
front end made of the scan (``image_fractions``), and the ground truth is the
same fractions after the cleanup (``fractions``, with their ``label``): a
vessel stays bright in the image and is CSF in the ground truth, as in
``matlab/mri_simulate.m``.  Without a cleanup the two are one array.  A
change of the anatomy changes both.

**The grid.**  The anatomy is made on the working grid of the front end, which
is affinely registered to MNI152NLin2009cAsym at 0.5 mm.  The registration
scales the head to the template, so every length and volume that defines the
anatomy (a thickness in mm, the ml of a sidecar) is one of that space and
not of the scan; the factor between the two volumes is in the sidecar
(``AffineVolumeScaling``).  An image is acquired on one of three grids:

- ``grid="scan"``, the default: the grid of the scan itself, with its
  voxels, its orientation and its field of view, in the world of the scan.
  The head has the size of the subject, the image lies on the scan voxel by
  voxel, and the ground truth is the mean of the fractions over every voxel
  of the scan.  What the scan shows beyond the working grid, the neck
  mostly, is taken from the scan.
- ``grid="working"``: the working grid, or one of another voxel size inside
  it (``resolution``).  The ground truth is the anatomy as it was made, with
  nothing averaged.
- ``grid="deepmriprep"``: the 336 x 384 x 336 voxels of the working grid the
  deepmriprep networks are trained on, so the output is training data as it
  is.  For that use the image can be acquired at the voxel size of a real
  scan and brought back onto the 0.5 mm grid (``resolution=1,
  upsample=True``): the ground truth keeps its 0.5 mm, and the image has the
  resolution and the correlated noise of a 1 mm scan that was resampled,
  which is what most real data on that grid is.

Wherever an image is resampled, it is by its k-space and a cubic B-spline,
and what swings below 0 next to an edge is set to 0; a ground truth is never
interpolated but averaged over the new voxels.
"""

from __future__ import annotations

import argparse
import sys
import zlib
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

from . import __version__, anatomy, artefacts, atlas, bids, cleanup, frontend, synthesis
from ._volume import grid_voxel_size

#: The part of the working grid (339 x 411 x 339) that deepmriprep's training
#: scripts use: 336 x 384 x 336 voxels.
DEEPMRIPREP_CROP = (slice(1, -2), slice(15, -12), slice(0, 336))

_CLEAN_DEFAULTS = {"vessels": 0.5, "dura": 4.0, "pve": True}


# ---------------------------------------------------------------- ground truth
def _clean_options(clean):
    """The options of the cleanup as a dict, or None for no cleanup."""
    if isinstance(clean, dict):
        unknown = set(clean) - set(_CLEAN_DEFAULTS)
        if unknown:
            raise ValueError(f"unknown cleanup option(s) {sorted(unknown)}; known: {sorted(_CLEAN_DEFAULTS)}")
        return {**_CLEAN_DEFAULTS, **clean}
    return dict(_CLEAN_DEFAULTS) if clean else None


def _atrophy_options(spec) -> list:
    """*spec* as a list of ``{"atlas", "region", "change", "mode", "tissues", "reach"}``, one per region."""
    if not spec:
        return []
    specs = [spec] if isinstance(spec, dict) else list(spec)
    out = []
    for item in specs:
        unknown = set(item) - {"atlas", "regions", "change", "mode", "tissues", "reach"}
        if unknown or "regions" not in item or "change" not in item:
            raise ValueError("an atrophy is a dict with 'regions' and 'change', and optionally "
                             "'atlas', 'mode', 'tissues' and 'reach'")
        regions_ = np.atleast_1d(item["regions"]).tolist()
        changes = np.atleast_1d(item["change"]).astype(float).tolist()
        if len(changes) == 1:
            changes *= len(regions_)
        if len(changes) != len(regions_):
            raise ValueError("an atrophy needs one change, or one per region")
        for region, change in zip(regions_, changes):
            out.append({"atlas": item.get("atlas", "Hammers"), "region": int(region),
                        "change": change, "mode": item.get("mode", "thinning"),
                        "tissues": tuple(item.get("tissues", ("GM",))), "reach": float(item.get("reach", 10.0))})
    return out


def _atlas_name(name) -> str:
    """The name of an atlas for a file name: that of the file, if it is one."""
    return bids.label(Path(str(name)).name.split(".")[0])


def _deform(atrophies: list, seg: dict, affine: np.ndarray, fractions: np.ndarray,
            image_fractions: np.ndarray, corrected: np.ndarray, verbose: bool) -> tuple:
    """Atrophy of one region after the other.

    The deformation is made from the ground truth, which knows what is no
    brain tissue, and applied to everything: the ground truth, the fractions
    of the image and the corrected scan.  From the second region on, the
    atlas has moved with the anatomy, and the Jacobian determinants are
    chained (:func:`t1sim.anatomy.compose`).

    Returns:
        ``(fractions, image_fractions, corrected, jacobian, done)``:
        the deformed volumes, the Jacobian determinant of all steps on the
        anatomy before the first, and what every step reached, for the
        sidecar.
    """
    cleaned = image_fractions is not fractions
    names = [item["atlas"] for item in atrophies]
    atlases = {name: atlas.to_grid(name, seg["iy"], seg["template_mat"]) for name in dict.fromkeys(names)}
    voxel_ml = float(np.prod(grid_voxel_size(affine))) / 1000
    jacobian = forward = None
    done = []
    for step, item in enumerate(atrophies):
        if verbose:
            print(f"Ground truth: atrophy of {item['atlas']} {item['region']}", flush=True)
        roi = anatomy.roi_from_atlas(atlases[item["atlas"]], [item["region"]])
        before = fractions
        result = anatomy.atrophy(fractions, corrected, affine, roi, item["change"], mode=item["mode"],
                                 tissues=item["tissues"], reach=item["reach"], verbose=verbose)
        fractions, corrected = result["fractions"], result["corrected"]
        if cleaned:
            image_fractions = anatomy.warp(image_fractions, result["displacement"])
        last = step == len(atrophies) - 1
        if jacobian is None:
            jacobian, forward = result["jacobian"], None if last else result["forward"]
        else:
            jacobian, forward = anatomy.compose(jacobian, forward, result)
        if not last:
            # the atlases the next regions come from move with the anatomy
            for name in set(names[step + 1:]):
                atlases[name] = anatomy.warp(atlases[name], result["displacement"], nearest=True)

        index = [frontend.TISSUES.index(tissue) for tissue in item["tissues"]]

        def in_region(volume):
            return sum(float(volume[..., i][roi].sum(dtype=np.float64)) for i in index)

        totals = (fractions.sum((0, 1, 2), dtype=np.float64) - before.sum((0, 1, 2), dtype=np.float64)) * voxel_ml
        done.append({"Atlas": str(item["atlas"]), "Region": item["region"], "Mode": item["mode"],
                     "Change": item["change"], "Tissues": list(item["tissues"]), "Reach": item["reach"],
                     "Achieved": round(result["achieved"], 4),
                     "ChangeInRegion": round(in_region(fractions) / in_region(before) - 1, 4),
                     "TissueChangeMl": {tissue: round(float(value), 3)
                                        for tissue, value in zip(frontend.TISSUES, totals)}})
    return fractions, (image_fractions if cleaned else fractions), corrected, jacobian, done


def ground_truth(scan, thickness=0, atrophy=None, clean=True, template_dir=None,
                 verbose: bool = False, seg: dict = None) -> dict:
    """The anatomy of a simulation.

    Args:
        scan: The T1w scan, a file.
        thickness: Give the cortex this thickness in mm; one value, or three
            for the mainly occipital parcels, the rest and the mainly frontal
            parcels.  0: keep the cortex.
        atrophy: A dict, or a list of dicts, each with ``regions`` (numbers
            of the atlas) and ``change`` (relative volume change of their
            grey matter, e.g. -0.1; one value or one per region), and
            optionally ``atlas`` (default ``"Hammers"``; a name of T1Prep's
            atlases or a file in template space), ``mode`` (``"thinning"``
            or ``"volume"``), ``tissues`` and ``reach``; see
            :func:`t1sim.anatomy.atrophy`.  Cannot be combined with
            *thickness*.
        clean: Clean the ground truth of vessels, dura and periventricular
            partial volume: True, False, or a dict with ``vessels``,
            ``dura`` and ``pve`` (:func:`t1sim.cleanup.detect`).  Where a
            *thickness* rebuilt the cortex nothing is cleaned: that cortex
            is its own ground truth.
        template_dir: Folder with the templates of the segmentation.
        verbose: Print the steps.
        seg: A result of ``frontend.segment(scan, deformation=True)`` to
            start from instead of running the front end.

    Returns:
        A dict:

        - ``fractions`` ((X, Y, Z, 3): CSF, GM, WM) and ``label``: the
          ground truth;
        - ``image_fractions``: what the image is rendered from, the ground
          truth before the cleanup (the same array without one);
        - ``corrected`` and ``ideal``: the bias-corrected scan and the image
          without artefacts, both with the changed anatomy;
        - ``means``: the tissue intensities; ``affine``: voxel-to-world
          matrix of the grid; ``iy`` and ``template_mat``: the deformation
          to the templates, of the anatomy before an atrophy;
        - ``jacobian``, with an atrophy: its Jacobian determinant, on the
          anatomy before the atrophy (see :func:`t1sim.anatomy.atrophy`);
        - ``scan``, ``tags`` (of the anatomy, for the file names),
          ``settings`` (what was asked for), ``cleaned`` and ``parameters``
          (for the sidecars);
        - ``scan_image`` and ``world``: the denoised scan on its own grid
          and the registration from the millimetres of the scan to those of
          the template, for an image on the grid of the scan.
    """
    atrophies = _atrophy_options(atrophy)
    thick = np.atleast_1d(np.asarray(thickness, dtype=float))
    if thick.any() and atrophies:
        raise ValueError("a simulated thickness and an atrophy cannot be combined")
    clean = _clean_options(clean)

    if seg is None:
        if verbose:
            print("Ground truth: segmentation of the scan", flush=True)
        seg = frontend.segment(scan, template_dir=template_dir, verbose=False, deformation=True)
    affine = np.asarray(seg["grid"].affine, dtype=np.float64)
    image_fractions, label, corrected = seg["fractions"], seg["label"], seg["corrected"]
    means = dict(seg["means"])
    voxel_ml = float(np.prod(grid_voxel_size(affine))) / 1000
    parameters, tags = {}, ""
    if "world" in seg:
        # the affine registration scales the head: volumes here, over volumes in the scan
        parameters["AffineVolumeScaling"] = round(float(abs(np.linalg.det(np.asarray(seg["world"])[:3, :3]))), 4)
    regions = None
    if thick.any() or clean:
        regions = atlas.to_grid("Neuromorphometrics", seg["iy"], seg["template_mat"])

    rebuilt = None
    if thick.any():
        if verbose:
            print("Ground truth: cortical thickness", flush=True)
        image_fractions, label, rebuilt = anatomy.thickness(image_fractions, label, regions, affine, thick,
                                                            verbose=verbose, return_rebuilt=True)
        tags += bids.thickness_tag(thick)
        parameters["Thickness"] = thick.tolist() if thick.size > 1 else float(thick[0])

    fractions = image_fractions
    if clean:
        if verbose:
            print("Ground truth: cleanup", flush=True)
        masks = cleanup.detect(label, seg["normalised"], affine, regions, **clean)
        if rebuilt is not None:
            # The rebuilt cortex holds no vessel and no dura, its mixtures
            # are those of its two neighbours, and it does not follow the
            # intensities of the scan any more, which the detection reads
            for mask in masks.values():
                mask &= ~rebuilt
        fractions, label = cleanup.apply(image_fractions, label, masks["vessels"] | masks["dura"],
                                         masks["pve"], means=means)
        parameters["LabelCleanup"] = {
            "BloodVessels": clean["vessels"], "DuraDistance": clean["dura"],
            "PeriventricularPVE": bool(clean["pve"]),
            "VolumeMl": {name: round(float(mask.sum()) * voxel_ml, 2) for name, mask in masks.items()},
            "Note": "the simulated image is rendered from the uncorrected fractions, "
                    "so it still shows what the ground truth rejects"}

    jacobian = None
    if atrophies:
        fractions, image_fractions, corrected, jacobian, done = _deform(
            atrophies, seg, affine, fractions, image_fractions, corrected, verbose)
        label = frontend.label_from_fractions(fractions)
        if len(done) == 1:
            kind = "Thin" if done[0]["Mode"] == "thinning" else "Vol"
            sign = "Plus" if done[0]["Change"] > 0 else ""
            tags += (f"{_atlas_name(done[0]['Atlas'])}Roi{done[0]['Region']}{kind}{sign}"
                     f"{bids.number(abs(100 * done[0]['Change']))}")
        else:
            tags += f"{_atlas_name(done[0]['Atlas'])}Multi"
        parameters["Atrophy"] = done

    parameters["TissueIntensities"] = {k: round(float(v), 3) for k, v in means.items()}
    settings = {"Thickness": thick.tolist() if thick.any() else None,
                "Atrophy": [{**item, "atlas": str(item["atlas"])} for item in atrophies],
                "Cleanup": clean}
    out = {
        "fractions": fractions, "label": label, "image_fractions": image_fractions,
        "corrected": corrected, "ideal": synthesis.render(image_fractions, means, corrected),
        "means": means, "affine": affine, "iy": seg["iy"], "template_mat": seg["template_mat"],
        "scan": str(scan) if scan is not None else "", "tags": bids.label(tags),
        "settings": settings, "cleaned": bool(clean), "parameters": parameters,
        "scan_image": seg.get("scan"), "world": seg.get("world"),
    }
    if jacobian is not None:
        out["jacobian"] = jacobian
    return out


# ----------------------------------------------------------------- acquisition
def _crop_affine(affine: np.ndarray, crop: tuple) -> np.ndarray:
    """The voxel-to-world matrix of a part of a grid."""
    shift = np.eye(4)
    shift[:3, 3] = [c.start or 0 for c in crop]
    return affine @ shift


def _seed_of(*parts) -> int:
    """A seed that follows from names: the same for the same file, another for another."""
    return zlib.crc32("_".join(str(p) for p in parts).encode()) & 0x7FFFFFFF


def _scan_id(scan) -> str:
    """What tells one scan from another for its seeds: its folder and its name.

    The name alone is not enough: ``subj01/T1.nii`` and ``subj02/T1.nii``
    would share their noise and their bias field.
    """
    if not scan:
        return "simulation"
    path = Path(str(scan)).resolve()
    return f"{path.parent.name}/{bids.prefix(path)}"


def _bias_field(truth: dict, shape, voxel, iy, percent: float, kind, scan_id: str) -> tuple:
    """The bias field of an acquisition, its tag and its entries in the sidecar."""
    sign = "Neg" if percent < 0 else ""
    parameters, seed_tag = {}, ""
    if isinstance(kind, str):
        bias = artefacts.template_bias_field(kind, iy, truth["template_mat"], percent)
        name = kind if len(kind) == 1 else "File"
    else:
        strength, seed = (list(kind) + [0])[:2]
        name = f"T{int(round(strength))}"
        if seed is None:
            seed = _seed_of(scan_id, "bias")            # a field of its own for every scan
        elif int(seed) != 0:
            seed_tag = f"S{int(seed)}"                  # another field, another file
        bias = artefacts.simulated_bias_field(shape, voxel, percent, strength, seed)
        parameters["BiasFieldSeed"] = int(seed)
    parameters["BiasFieldType"] = name
    return bias, f"Rf{sign}{bids.number(abs(percent))}{name}{seed_tag}", parameters


def _local_gain(reference: np.ndarray, raw: np.ndarray, valid: np.ndarray, voxel: np.ndarray,
                width: float = 5.0) -> np.ndarray:
    """The smooth factor that brings one image to the intensities of another.

    Where both are known it is the ratio of their local means, on a lattice
    of about 4 mm; beyond, it is continued with its nearest value.  The
    means are local enough for the factor to be right at the border of what
    is known, where the two images meet: a bias field changes within the
    width of a wider mean, and the images then met with a step.

    Args:
        reference: (X, Y, Z) image with the intensities that are wanted.
        raw: (X, Y, Z) image of the same object that is to be brought there.
        valid: (X, Y, Z) where *reference* is known.
        voxel: Voxel size in mm.
        width: Standard deviation in mm of the local means.

    Returns:
        (X, Y, Z) float32, the factor for *raw*.
    """
    step = np.maximum(1, np.round(4.0 / np.asarray(voxel))).astype(int)
    bright = raw > 0.1 * float(np.percentile(raw[::4, ::4, ::4], 99))
    usable = (valid & bright & (reference > 0)).astype(np.float32)

    def local_mean(volume):
        coarse = ndimage.uniform_filter(volume * usable, size=tuple(step), mode="nearest")
        coarse = coarse[tuple(slice(n // 2, None, n) for n in step)]
        return ndimage.gaussian_filter(coarse, width / (voxel * step), mode="constant")

    wanted, there = local_mean(reference), local_mean(raw)
    known = there > 1e-3 * float(there.max())
    if not known.any():
        return np.ones(raw.shape, dtype=np.float32)
    gain = np.where(known, wanted / np.maximum(there, 1e-12), 0).astype(np.float32)
    # a few voxels at the rim of what is known give any ratio: keep to what most of it has
    gain = np.clip(gain, *np.percentile(gain[known], [2, 98]))
    nearest = ndimage.distance_transform_edt(~known, return_distances=False, return_indices=True)
    gain = gain[tuple(nearest)]
    return ndimage.affine_transform(gain, 1.0 / step, offset=-(step // 2) / step, output_shape=raw.shape,
                                    order=1, mode="nearest")


def _on_the_scan(truth: dict, image: np.ndarray, affine: np.ndarray, resolution, contrast: float,
                 white: float, bias) -> tuple:
    """Acquire an image of the working grid on the grid of the scan.

    Inside the working grid the image is acquired as on any other grid: its
    k-space is cut to that of the scan, and it is read at the voxel centres
    of the scan.  A scan usually shows more than the working grid holds, the
    neck and more air.  There the image is the denoised scan itself, brought
    to the intensities of the bias-corrected scan next to it and given the
    contrast and the bias field of the simulation; the two are blended over
    3 mm at the border of the working grid.

    Returns:
        ``(image, acquisition, bias)``: the image on the grid of the scan,
        the :class:`t1sim.artefacts.ObliqueAcquisition` that made it, with
        the affine of the scan in its own world as ``scan_affine``, and the
        bias field on that grid (or None).
    """
    scan = truth.get("scan_image")
    if scan is None or truth.get("world") is None:
        raise ValueError("the grid of the scan needs the scan: make the ground truth from a result of "
                         "frontend.segment that holds 'scan' and 'world', or ask for grid='working'")
    shape, scan_affine = tuple(scan.shape[:3]), np.asarray(scan.affine, dtype=np.float64)
    rest = np.asarray(scan.dataobj, dtype=np.float32).reshape(shape)
    if resolution is not None:
        # another voxel size on the axes and in the field of view of the scan
        coarser = artefacts.Acquisition(shape, scan_affine, resolution)
        shape, scan_affine, rest = coarser.shape, coarser.affine, np.maximum(coarser.acquire(rest), 0)
    acquisition = artefacts.ObliqueAcquisition(image.shape, affine, shape, np.asarray(truth["world"]) @ scan_affine)
    acquisition.scan_affine = scan_affine
    # the cut of k-space and the spline both swing below 0 next to an edge: a magnitude image does not
    image = np.maximum(acquisition.acquire(image), 0)
    if bias is not None:
        bias = acquisition.sample(bias, extend=True)        # beyond the working grid: its nearest value

    weight = acquisition.coverage(ramp=6.0)
    if float(weight.min()) < 1:
        inside = weight >= 1
        gain = _local_gain(acquisition.sample(truth["corrected"]), rest, inside, grid_voxel_size(scan_affine))
        rest = np.maximum(rest, 0) * gain
        if contrast != 1:
            rest = artefacts.change_contrast(rest, contrast, white)
        if bias is not None:
            rest *= bias
        image = weight * image + (1 - weight) * rest
    return image, acquisition, bias


def acquire(truth: dict, snr: float = 40.0, noise_percent: float = 0.0, bias_percent: float = 30.0,
            bias_type=(2, 0), contrast: float = 1.0, resolution=None, upsample: bool = False,
            motion=0, ringing=0, seed: int = None, grid: str = "scan") -> dict:
    """One simulated image of a ground truth.

    The steps follow the acquisition:

    1. The contrast belongs to the tissue and the sequence: the image is
       mixed again from the tissue intensities of that contrast.
    2. The bias field belongs to the coil and multiplies that image.
    3. The voxel size and
    4. the movement and the ringing belong to the sampling of k-space.
    5. The noise comes last, on the grid of the acquisition.
    6. Optionally the image is brought back onto the grid of the ground
       truth.

    Args:
        truth: Result of :func:`ground_truth`.  It is not changed, except
            that the ground truth on the grid of the scan is kept in it for
            the next image on that grid.
        snr: Rician noise for this signal-to-noise ratio in white matter:
            its standard deviation in each of the two channels is the
            intensity of white matter divided by *snr*.  0: use
            *noise_percent*.
        noise_percent: Gaussian noise with this standard deviation in percent
            of the white matter intensity, if *snr* is 0.
        bias_percent: Peak-to-peak amplitude of the bias field in percent; 0
            for none, negative to invert the field.
        bias_type: ``(strength, seed)`` for a simulated field, or ``"A"``,
            ``"B"``, ``"C"`` or a file for a predefined one.  A seed of None
            gives every scan a field of its own; a number gives all scans
            the same.
        contrast: Exponent of the contrast change
            (:func:`t1sim.artefacts.change_contrast`): white matter keeps
            its intensity, and the ratio of two tissue intensities is raised
            to it.  1 leaves the contrast.
        resolution: Voxel size of the acquisition in mm, one value or three.
            None: that of the grid.  With ``grid="scan"`` the millimetres
            are those of the scan, along its axes and in its field of view.
        upsample: Bring the acquired image back onto the grid, so the image
            and the ground truth both have the voxel size of the grid.
            False: the ground truth is averaged to the voxels of the
            acquisition.  Not with ``grid="scan"``.
        motion, ringing: Movement artefacts and ringing: a severity or a
            dict (:mod:`t1sim.kspace`); 0 for none.
        seed: Seed of the noise and of the movement.  None: derived from
            the folder and the name of the scan and from the settings, so
            every output file has noise of its own that is the same in every
            run.  A number gives every scan the same noise pattern, which is
            useful for a comparison and harmful in training data.
        grid: ``"scan"`` for the grid of the scan, in the world of the scan
            (the default); ``"working"`` for the 0.5 mm grid in template
            space the anatomy is made on; or ``"deepmriprep"`` for the 336 x
            384 x 336 voxels of that grid its networks are trained on (see
            the module notes).

    Returns:
        A dict: ``image``; ``label`` and ``fractions`` ((X, Y, Z, 3): CSF,
        GM, WM), the ground truth on the grid of the image; ``affine`` of
        that grid; ``bias`` (the field, or None); ``jacobian`` (of an
        atrophy, on that grid, or None); ``tags`` (of the acquisition) and
        ``parameters`` (for the sidecar).
    """
    if grid not in ("working", "deepmriprep", "scan"):
        raise ValueError("grid is 'working', 'deepmriprep' or 'scan'")
    if grid == "scan" and upsample:
        raise ValueError("upsample brings the image back onto the 0.5 mm grid of the ground truth; "
                         "on the grid of the scan the ground truth has the voxels of the image")
    if snr < 0 or noise_percent < 0:
        raise ValueError("snr and noise_percent cannot be negative")
    if not contrast > 0:
        raise ValueError("contrast is the exponent of a power law and has to be positive")
    affine = truth["affine"]
    image, label, fractions = truth["ideal"], truth["label"], truth["fractions"]
    image_fractions, corrected = truth["image_fractions"], truth["corrected"]
    iy, jacobian = truth["iy"], truth.get("jacobian")
    if grid == "deepmriprep":
        if image.shape != (339, 411, 339):
            raise ValueError(f"the deepmriprep grid is cut from the 339 x 411 x 339 grid, not {image.shape}")
        crop = DEEPMRIPREP_CROP
        image, label, fractions, iy = image[crop], label[crop], fractions[crop], iy[crop]
        image_fractions, corrected = image_fractions[crop], corrected[crop]
        jacobian = None if jacobian is None else jacobian[crop]
        affine = _crop_affine(affine, crop)
    voxel = grid_voxel_size(affine)
    scan_id = _scan_id(truth["scan"])
    white = float(truth["means"]["WM"])
    parameters = {"Grid": grid, "Space": "scan" if grid == "scan" else bids.SPACE}
    # The parts of the desc label.  They are joined in the order of
    # mri_simulate.m, which is not the order of the steps.
    tag = dict.fromkeys(("noise", "bias", "acquisition", "contrast", "kspace", "seed"), "")

    if snr > 0:
        tag["noise"] = f"snr{bids.number(snr)}"
        parameters["SNR"] = float(snr)
    elif noise_percent > 0:
        tag["noise"] = f"pn{bids.number(noise_percent)}"
        parameters["NoiseFraction"] = noise_percent / 100

    if contrast != 1:
        # A voxel of two tissues stays the mixture of their intensities
        means = {tissue: artefacts.change_contrast(value, contrast, white)
                 for tissue, value in truth["means"].items()}
        image = synthesis.render(image_fractions, means, artefacts.change_contrast(corrected, contrast, white))
        tag["contrast"] = {0.5: "ConLow", 1.5: "ConHigh"}.get(contrast, f"Con{bids.number(contrast)}")
        parameters["ContrastChange"] = float(contrast)
        parameters["TissueIntensities"] = {tissue: round(value, 3) for tissue, value in means.items()}

    bias = None
    parameters["BiasFieldStrength"] = float(bias_percent)
    if bias_percent != 0:
        bias, tag["bias"], described = _bias_field(truth, image.shape, voxel, iy, bias_percent, bias_type, scan_id)
        parameters.update(described)
        image = image * bias                            # the field has a mean of 1

    acquisition = None
    if grid == "scan":
        image, acquisition, bias = _on_the_scan(truth, image, affine, resolution, contrast, white, bias)
        # The ground truth on that grid is the same for every image of it, and
        # the most expensive part: it is kept with the ground truth
        kept = truth.setdefault("_on_the_scan", {})
        key = None if resolution is None else tuple(np.atleast_1d(resolution).astype(float).tolist())
        # ... as long as it is still made from these arrays: a copy of the dict shares what is kept
        if key not in kept or kept[key][0] is not fractions or kept[key][1] is not jacobian:
            averaged = acquisition.average(fractions)
            kept[key] = (fractions, jacobian, averaged, frontend.label_from_fractions(averaged),
                         # nothing is deformed beyond the working grid
                         None if jacobian is None else acquisition.average(jacobian - 1) + 1)
        fractions, label, jacobian = kept[key][2:]
        if resolution is not None:
            parameters["AcquisitionVoxelSize"] = grid_voxel_size(acquisition.scan_affine).tolist()
    elif resolution is not None:
        acquisition = artefacts.Acquisition(image.shape, affine, resolution)
        image = np.maximum(acquisition.acquire(image), 0)     # as on the grid of the scan: nothing below 0
        parameters["AcquisitionVoxelSize"] = acquisition.voxel.tolist()
        if upsample:
            sizes = acquisition.voxel if np.ptp(acquisition.voxel) > 1e-6 else acquisition.voxel[:1]
            tag["acquisition"] = "Acq" + "x".join(bids.number(round(v, 2)) for v in sizes) + "mm"
        else:
            label, fractions = acquisition.average(label), acquisition.average(fractions)
            if bias is not None:
                bias = acquisition.average(bias)
            if jacobian is not None:
                jacobian = acquisition.average(jacobian)
    if grid == "scan":
        acquired_affine = acquisition.scan_affine           # from here on in the millimetres of the scan
    else:
        acquired_affine = acquisition.affine if acquisition is not None else affine

    derived = seed is None
    if derived:
        seed = _seed_of(scan_id, tag["noise"] + tag["bias"] + tag["acquisition"] + tag["contrast"],
                        truth["tags"], resolution, upsample, motion, ringing, grid)
    else:
        tag["seed"] = f"Seed{int(seed)}"                # other noise, other file
    parameters["Seed"] = int(seed)
    motion_seed, noise_seed = np.random.SeedSequence(int(seed)).spawn(2)
    if motion or ringing:
        from . import kspace

        ring = kspace.resolve_ringing(ringing, acquired_affine) if ringing else None
        if motion:
            move = kspace.resolve_motion(motion, acquired_affine, image.shape)
            image, parameters["Motion"] = kspace.add_motion(
                image, move, acquired_affine, int(motion_seed.generate_state(1)[0]), ringing=ring)
            tag["kspace"] = f"Motion{bids.number(move['severity'])}"
        else:
            image = kspace.add_ringing(image, ring)
        if ring is not None:
            parameters["Ringing"] = dict(ring)
            tag["kspace"] += f"Ringing{bids.number(ring['strength'])}" + ("Gibbs" if ring["type"] == "gibbs" else "")

    sigma = white / snr if snr > 0 else noise_percent / 100 * white
    parameters["NoiseSD"] = round(float(sigma), 4)
    image = artefacts.add_noise(image, sigma, rician=snr > 0, seed=np.random.default_rng(noise_seed))

    out_affine = acquired_affine
    if acquisition is not None and upsample:
        image = np.maximum(acquisition.to_fine(image), 0)
        out_affine = affine
    parameters["VoxelSize"] = [round(float(v), 4) for v in grid_voxel_size(out_affine)]
    return {"image": image.astype(np.float32, copy=False), "label": label, "fractions": fractions,
            "affine": out_affine, "bias": bias, "jacobian": jacobian,
            "tags": bids.label("".join(tag.values())), "parameters": parameters}


# ----------------------------------------------------------------------- files
def _refuse_another_ground_truth(path: Path, settings: dict) -> None:
    """Stop before the ground truth of another simulation is replaced.

    Images that differ in the acquisition alone share the files of one
    ground truth, and its name does not hold every setting.  Replacing it
    with one of another grid or another cleanup would leave the images that
    are already there with a ground truth that is not theirs.
    """
    if not path.is_file():
        return
    known = (bids.read_sidecar(path) or {}).get("GroundTruthSettings")
    if known is None:                                   # written without this record: the grid at least
        if tuple(nib.load(str(path)).shape) != tuple(settings["Shape"]):
            raise FileExistsError(f"{path} exists on another grid; write to another folder")
        return
    new = bids.as_json(settings)
    differing = sorted(key for key in set(known) | set(new) if known.get(key) != new.get(key))
    if differing:
        raise FileExistsError(
            f"{path} is the ground truth of another simulation of this scan (it differs in "
            f"{', '.join(differing)}), and the images that belong to it would lose it: "
            "write to another folder (out_dir)")


def write(truth: dict, acquired: dict, out_dir=None, derivatives: bool = True,
          save_bias: bool = False, save_jacobian: bool = True, compress: bool = True) -> dict:
    """Store an acquisition and its ground truth.

    Args:
        truth, acquired: Results of :func:`ground_truth` and :func:`acquire`.
        out_dir: Folder of the outputs; default: see
            :func:`t1sim.bids.output_folder`.
        derivatives: Used when *out_dir* is None.
        save_bias: Also write the bias field, as ``..._desc-<tags of the
            image>Biasfield_T1w``.
        save_jacobian: Also write the Jacobian determinant of an atrophy.
        compress: Write ``.nii.gz``; False: ``.nii``.

    Returns:
        The files that were written: ``{"T1w": ..., "dseg": ..., "probseg":
        {"CSF": ..., "GM": ..., "WM": ...}, ...}``.

    Raises:
        FileExistsError: The folder holds a ground truth of the same name
            that was made with other settings or on another grid.
    """
    scan = truth["scan"] or "simulation_T1w.nii.gz"
    if out_dir is None:
        folder, pipeline, source = bids.output_folder(scan, derivatives)
    else:
        folder, pipeline, source = Path(out_dir), None, ""
    ext = ".nii.gz" if compress else ".nii"
    affine = acquired["affine"]
    on_scan = acquired["parameters"]["Grid"] == "scan"
    if on_scan:
        # as mri_simulate.m names what stays on the grid of the scan: no space, and a
        # resolution only if it is not that of the scan
        base = bids.prefix(scan)
        if "AcquisitionVoxelSize" in acquired["parameters"]:
            base += f"_{bids.resolution_entity(grid_voxel_size(affine))}"
        code = 1 if truth.get("scan_image") is None else int(truth["scan_image"].header["sform_code"]) or 1
    else:
        base = f"{bids.prefix(scan)}_space-{bids.SPACE}_{bids.resolution_entity(grid_voxel_size(affine))}"
        code = bids.TEMPLATE_CODE
    anat = truth["tags"]
    image_desc = acquired["tags"] + anat or "simu"
    truth_desc = anat + ("Clean" if truth["cleaned"] else "")
    truth_base = base + (f"_desc-{truth_desc}" if truth_desc else "")
    settings = {**truth["settings"], "Grid": acquired["parameters"]["Grid"],
                "Shape": list(acquired["label"].shape)}
    _refuse_another_ground_truth(folder / f"{truth_base}_dseg{ext}", settings)
    if pipeline is not None:
        bids.write_dataset_description(pipeline)
    common = {"GeneratedBy": [bids.GENERATED_BY]}
    if source:
        common["Sources"] = [source]
    files = {}

    files["T1w"] = bids.save(acquired["image"], affine, folder / f"{base}_desc-{image_desc}_T1w{ext}", code=code)
    bids.write_sidecar(files["T1w"], {**common, "SimulationParameters": {**truth["parameters"],
                                                                         **acquired["parameters"]}})
    # the label map as CAT12 stores it: 16 bit in steps of 1/765
    files["dseg"] = bids.save(acquired["label"], affine, folder / f"{truth_base}_dseg{ext}",
                              scale=1 / 765, dtype=np.int16, code=code)
    bids.write_sidecar(files["dseg"], {
        **common, "Manual": False,
        "Description": "Ground truth partial volume label image. Values are continuous and interpolate "
                       "between the tissue labels, e.g. 2.5 is an equal mixture of GM and WM. The tissue "
                       "fractions it is a view of are the probseg files.",
        "Labels": {"CSF": 1, "GM": 2, "WM": 3}, "SimulationParameters": truth["parameters"],
        "GroundTruthSettings": settings})
    files["probseg"] = {}
    for index, tissue in enumerate(frontend.TISSUES):
        files["probseg"][tissue] = bids.save(
            np.clip(acquired["fractions"][..., index], 0, 1), affine,
            folder / f"{truth_base}_label-{tissue}_probseg{ext}", scale=1 / 255, dtype=np.uint8, code=code)
    if save_bias and acquired["bias"] is not None and acquired["bias"].shape == acquired["image"].shape:
        files["bias"] = bids.save(acquired["bias"], affine,
                                  folder / f"{base}_desc-{acquired['tags'] + anat}Biasfield_T1w{ext}", code=code)
    if save_jacobian and acquired.get("jacobian") is not None:
        files["jacobian"] = bids.save(acquired["jacobian"], affine, folder / f"{truth_base}_jacobian{ext}",
                                      code=code)
        bids.write_sidecar(files["jacobian"], {
            **common,
            "Description": "Jacobian determinant of the deformation that makes the atrophy: the factor by "
                           "which the tissue at a voxel changes its volume, below 1 for a loss. It lies on "
                           "the anatomy BEFORE the atrophy, i.e. on the simulation of this scan without "
                           "atrophy, and not on the label map of the same name, in which that tissue has "
                           "moved by a fraction of a mm.",
            "SimulationParameters": {"Atrophy": truth["parameters"].get("Atrophy")}})
    return files


def simulate(scan, out_dir=None, derivatives: bool = True, thickness=0, atrophy=None, clean=True,
             template_dir=None, truth: dict = None, save_bias: bool = False, save_jacobian: bool = True,
             compress: bool = True, verbose: bool = False, **acquisition) -> dict:
    """Simulate one image of a scan and write it with its ground truth.

    Args:
        scan: The T1w scan, a file.
        out_dir, derivatives, save_bias, save_jacobian, compress: See
            :func:`write`.
        thickness, atrophy, clean, template_dir: See :func:`ground_truth`.
        truth: A result of :func:`ground_truth` for this scan, to make
            another acquisition of the same anatomy without running the
            front end again.  *thickness*, *atrophy* and *clean* are then
            not used.
        verbose: Print the steps.
        **acquisition: The arguments of :func:`acquire` (``snr``,
            ``bias_percent``, ``resolution``, ``motion``, ``grid``, ...).

    Returns:
        ``{"files": ..., "truth": ..., "acquired": ...}``.
    """
    if truth is None:
        truth = ground_truth(scan, thickness=thickness, atrophy=atrophy, clean=clean,
                             template_dir=template_dir, verbose=verbose)
    if verbose:
        print("Acquisition", flush=True)
    acquired = acquire(truth, **acquisition)
    files = write(truth, acquired, out_dir=out_dir, derivatives=derivatives, save_bias=save_bias,
                  save_jacobian=save_jacobian, compress=compress)
    if verbose:
        print(f"Wrote {files['T1w']}\n      {files['dseg']}", flush=True)
    return {"files": files, "truth": truth, "acquired": acquired}


# ---------------------------------------------------------------- command line
def _bias_type_argument(values: list):
    """``--bias-type``: strength and seed of a simulated field, or the name or file of a stored one."""
    try:
        strength = float(values[0])
    except ValueError:
        if len(values) > 1:
            raise SystemExit("t1sim: --bias-type takes a strength and a seed, or one name or file") from None
        return values[0]
    if len(values) > 2:
        raise SystemExit("t1sim: --bias-type takes a strength and a seed, or one name or file")
    seed = values[1] if len(values) > 1 else "0"
    return strength, None if seed == "scan" else int(seed)


def main(argv=None) -> int:
    """Command line: ``t1sim scan.nii.gz [options]``."""
    parser = argparse.ArgumentParser(
        prog="t1sim", description="Simulate T1-weighted images with an exact ground truth from a real scan.")
    parser.add_argument("scans", nargs="+", help="T1w scan(s), NIfTI")
    parser.add_argument("--version", action="version", version=f"t1sim {__version__}")
    out = parser.add_argument_group("output")
    out.add_argument("--out-dir", help="folder of the outputs (default: derivatives/t1sim-<version>)")
    out.add_argument("--no-derivatives", action="store_true", help="write next to the scan")
    out.add_argument("--grid", choices=("scan", "working", "deepmriprep"), default="scan",
                     help="output grid: that of the scan (default); the 0.5 mm working grid in template "
                          "space; or the 336x384x336 part of it deepmriprep trains on")
    out.add_argument("--save-bias", action="store_true", help="also write the bias field")
    anat = parser.add_argument_group("anatomy")
    anat.add_argument("--thickness", type=float, nargs="+", default=[0],
                      help="cortical thickness in mm: one value, or three (occipital, rest, frontal)")
    anat.add_argument("--atrophy", nargs=3, action="append", metavar=("ATLAS", "REGION", "CHANGE"),
                      help="atrophy of a region, e.g. Hammers 28 -0.1; can be given several times")
    anat.add_argument("--atrophy-mode", choices=("thinning", "volume"), default="thinning")
    anat.add_argument("--no-clean", action="store_true", help="keep vessels and dura in the ground truth")
    acq = parser.add_argument_group("acquisition; several values of --snr or --seed give several images of one ground truth")
    acq.add_argument("--snr", type=float, nargs="+", default=[40.0], help="SNR in white matter, Rician noise (40)")
    acq.add_argument("--noise-percent", type=float, default=0.0, help="Gaussian noise in %% of WM, with --snr 0")
    acq.add_argument("--bias", type=float, default=30.0, help="bias field amplitude in %% (30; 0 for none)")
    acq.add_argument("--bias-type", nargs="+", default=["2", "0"],
                     help="strength and seed of a simulated field (2 0; 'scan' as seed gives every "
                          "scan its own field), or A, B, C or the file of a field in template space")
    acq.add_argument("--contrast", type=float, default=1.0, help="exponent of the contrast change (1)")
    acq.add_argument("--resolution", type=float, nargs="+", help="voxel size of the acquisition in mm")
    acq.add_argument("--upsample", action="store_true",
                     help="bring the acquired image back onto the grid (for training data)")
    acq.add_argument("--motion", type=float, default=0, help="severity of the movement artefacts (0)")
    acq.add_argument("--ringing", type=float, default=0, help="strength of the ringing (0)")
    acq.add_argument("--seed", type=int, nargs="+", default=[None],
                     help="seed of noise and movement (default: every output file its own)")
    parser.add_argument("--verbose", action="store_true")
    if argv is None and len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        return 1
    args = parser.parse_args(argv)

    bias_type = _bias_type_argument(args.bias_type)
    atrophy = [{"atlas": a, "regions": [int(r)], "change": float(c), "mode": args.atrophy_mode}
               for a, r, c in (args.atrophy or [])]
    resolution = None if not args.resolution else (args.resolution[0] if len(args.resolution) == 1
                                                   else tuple(args.resolution))
    thickness = args.thickness[0] if len(args.thickness) == 1 else tuple(args.thickness)
    for scan in args.scans:
        truth = ground_truth(scan, thickness=thickness, atrophy=atrophy, clean=not args.no_clean,
                             verbose=args.verbose)
        for snr in args.snr:
            for seed in args.seed:
                result = simulate(scan, out_dir=args.out_dir, derivatives=not args.no_derivatives,
                                  truth=truth, save_bias=args.save_bias, verbose=args.verbose,
                                  snr=snr, noise_percent=args.noise_percent, bias_percent=args.bias,
                                  bias_type=bias_type, contrast=args.contrast, resolution=resolution,
                                  upsample=args.upsample, motion=args.motion, ringing=args.ringing,
                                  seed=seed, grid=args.grid)
                print(result["files"]["T1w"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
