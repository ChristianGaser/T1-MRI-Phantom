"""Maps of the template space on the working grid.

The cleanup of the label, the thickness simulation, the atrophy of a region
and the predefined bias fields all need something that is defined in the
template space: an atlas, a prior, a field.  The unified segmentation
deforms its templates onto the scan, and its inverse deformation ``iy``
says for every voxel of the working grid where it lies in the templates
(``frontend.segment(..., deformation=True)``).  :func:`to_grid` follows it,
so a map arrives through the same nonlinear deformation the segmentation
found and not through the affine registration alone.

Atlases are taken from T1Prep's template folder, where every atlas has a
table of its regions next to it (``<name>.csv``).  The region numbers of the
Neuromorphometrics atlas that t1sim relies on are collected here.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates

#: Region numbers of the Neuromorphometrics atlas.  Every cortical parcel has
#: a number from 100 on, with the right hemisphere even and the left odd;
#: everything else lies below 100.
CORTEX_FROM = 100
VENTRICLES = (4, 11, 49, 50, 51, 52)
CEREBRAL_WM = (44, 45)                       # right, left
BASAL_GANGLIA = (23, 30, 36, 37, 55, 56, 57, 58)
THALAMUS = (59, 60, 61, 62)                  # with the ventral diencephalon
HIPPOCAMPUS = (31, 32, 47, 48)               # with the amygdala
BRAINSTEM = (35,)
CEREBELLUM = (38, 39, 40, 41, 71, 72, 73)
#: Cortical parcels that are mainly occipital and mainly frontal: the regions
#: of a thickness that differs between the lobes.
OCCIPITAL = (108, 109, 114, 115, 128, 129, 134, 135, 144, 145, 148, 149, 156, 157, 158, 159,
             160, 161, 170, 171, 176, 177, 196, 197)
FRONTAL = (102, 103, 104, 105, 120, 121, 132, 133, 136, 137, 146, 147, 152, 153, 154, 155,
           162, 163, 164, 165, 172, 173, 178, 179, 190, 191, 202, 203, 204, 205)


def template_folder() -> Path:
    """T1Prep's folder of templates and atlases."""
    try:
        from t1prep.utils import TEMPLATE_PATH_T1PREP
    except ImportError as exc:
        raise RuntimeError(f"The atlases are taken from T1Prep, which is not installed: {exc}") from exc
    return Path(TEMPLATE_PATH_T1PREP)


def find(name) -> Path:
    """The file of an atlas or template map.

    Args:
        name: A file, or the name of a map in T1Prep's template folder
            without its extension; upper and lower case do not matter there.
    """
    path = Path(os.fspath(name))
    if path.is_file():
        return path
    folder = template_folder()
    wanted = {f"{path.name.lower()}{ext}" for ext in (".nii.gz", ".nii")}
    for entry in sorted(folder.iterdir()):
        if entry.name.lower() in wanted:
            return entry
    raise FileNotFoundError(f"'{name}' is neither a file nor a map in {folder}")


def regions(name) -> dict:
    """The regions of an atlas: number -> name, from the table next to it."""
    path = find(name)
    table = path.with_name(path.name.split(".")[0] + ".csv")
    if not table.is_file():
        raise FileNotFoundError(f"no table of regions next to {path}")
    with table.open(newline="") as handle:
        return {int(row["ROIid"]): row["ROIname"] for row in csv.DictReader(handle, delimiter=";")}


def to_grid(image, iy: np.ndarray, template_mat: np.ndarray, order: int = 0,
            fill: float = 0.0) -> np.ndarray:
    """A map of the template space, on the working grid.

    Args:
        image: The map: a nibabel image, a file, or a name :func:`find`
            knows.  It may lie on any grid of the template space.
        iy: (X, Y, Z, 3), for every voxel of the working grid the voxel of
            the segmentation's templates it maps to (0-based, NaN where
            undefined).
        template_mat: Voxel-to-world matrix of those templates.
        order: 0 for an atlas (the nearest voxel), 1 for a continuous map
            (linear interpolation).
        fill: Value where the deformation is undefined or the map ends.

    Returns:
        (X, Y, Z) array of the map's type (float32 for ``order=1``).
    """
    if order not in (0, 1):
        raise ValueError("order must be 0 (nearest) or 1 (linear)")
    if not isinstance(image, nib.spatialimages.SpatialImage):
        image = nib.load(find(image))
    data = np.asarray(image.dataobj)
    if data.ndim != 3:
        raise ValueError(f"the map must be a 3-D volume, not {data.shape}")
    # template voxel -> world -> voxel of the map
    mat = np.linalg.inv(image.affine) @ np.asarray(template_mat, dtype=np.float64)
    rot, shift = mat[:3, :3].astype(np.float32), mat[:3, 3].astype(np.float32)
    shape = np.array(data.shape)
    out = np.full(iy.shape[:3], fill, dtype=data.dtype if order == 0 else np.float32)
    step = 24                                   # planes at a time: bounds the temporaries
    for x0 in range(0, iy.shape[0], step):
        coords = iy[x0:x0 + step] @ rot.T + shift
        defined = np.isfinite(coords).all(-1)
        if order == 0:
            index = np.rint(np.nan_to_num(coords)).astype(np.int64)
            defined &= ((index >= 0) & (index < shape)).all(-1)
            index = index[defined]
            out[x0:x0 + step][defined] = data[index[:, 0], index[:, 1], index[:, 2]]
        else:
            values = map_coordinates(data.astype(np.float32, copy=False), coords[defined].T,
                                     order=1, mode="constant", cval=np.nan)
            values[~np.isfinite(values)] = fill
            out[x0:x0 + step][defined] = values
    return out


def fill_from_nearest(atlas: np.ndarray, vx, hole: np.ndarray = None) -> np.ndarray:
    """Give every voxel without a region that of the nearest voxel with one.

    Args:
        atlas: (X, Y, Z) region numbers, 0 where there is none.
        vx: Voxel size in mm.
        hole: Voxels whose region is dropped first, so that they are filled
            from their surroundings as well.
    """
    from ._volume import nearest_index

    atlas = atlas.copy()
    if hole is not None:
        atlas[hole] = 0
    if not (atlas > 0).any():
        return atlas
    return atlas[nearest_index(atlas > 0, vx)]
