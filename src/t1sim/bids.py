"""Names, folders and sidecars of the outputs, after BIDS.

The scheme is that of ``matlab/mri_simulate.m``, so that tools written for
its output read the output of t1sim as well (``deepmriprep-train`` pairs
every simulated image with its ground truth by these names):

- A file name is the entities of the input, ``space`` and ``res`` of the
  output grid, one ``desc`` label and a suffix.
- The simulated image is a ``T1w``.  Its ``desc`` holds the tags of the
  acquisition (noise, bias field, contrast, movement, ringing) followed by
  the tags of the anatomy (atrophy, thickness).
- The ground truth is a ``dseg`` (the label map; continuous, which its
  sidecar says) and one ``probseg`` per tissue (the fractions).  Their
  ``desc`` holds the tags of the anatomy only, plus ``Clean``: runs that
  differ in the acquisition alone share one ground truth.

BIDS labels are alphanumeric, so a decimal point becomes ``p`` and
everything else that is not a letter or a digit is dropped.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import nibabel as nib
import numpy as np

from . import __version__

GENERATED_BY = {"Name": "t1sim", "Version": __version__}
SPACE = "MNI152NLin2009cAsym"


def label(text) -> str:
    """*text* as a BIDS label: ``p`` for the decimal point, letters and digits only."""
    return re.sub(r"[^a-zA-Z0-9]", "", str(text).replace(".", "p"))


def number(value: float) -> str:
    """A number as short as MATLAB's ``%g`` prints it, as a BIDS label."""
    return label(f"{value:g}")


def strip_extension(name: str) -> str:
    for ext in (".nii.gz", ".nii"):
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return name


def prefix(scan) -> str:
    """The entities of the input that the outputs keep.

    A ``_T1w`` suffix and a ``desc`` entity of the input are dropped; the
    whole name of an input that is not named after BIDS is kept.
    """
    name = strip_extension(Path(str(scan)).name)
    if name.endswith("_T1w"):
        return re.sub(r"_desc-[a-zA-Z0-9]+$", "", name[:-4])
    return name


def resolution_entity(voxel) -> str:
    """``res-0p5mm``, or ``res-0p5x0p5x1p5mm`` for voxels that are not isotropic."""
    voxel = [float(f"{v:.2f}") for v in voxel]
    if max(voxel) - min(voxel) < 1e-6:
        return f"res-{number(voxel[0])}mm"
    return "res-" + "x".join(number(v) for v in voxel) + "mm"


def output_folder(scan, derivatives: bool = True) -> tuple:
    """Where the outputs of a scan go.

    Args:
        scan: The input file.
        derivatives: Write into ``derivatives/t1sim-<version>`` of the
            dataset: at the root of a BIDS dataset (the folder above the
            first ``sub-*`` folder), with the folders of subject and session
            kept below it, or next to the folder of the scan otherwise.
            False: next to the scan.

    Returns:
        ``(folder, pipeline_folder, source)``: the folder of the outputs, the
        folder of the pipeline (None without *derivatives*) and the BIDS URI
        of the scan (empty if it is not inside a BIDS dataset).
    """
    scan = Path(str(scan)).resolve()
    if not derivatives:
        return scan.parent, None, ""
    pipeline = f"t1sim-{__version__}"
    parts = scan.parent.parts
    subject = next((i for i, part in enumerate(parts) if part.startswith("sub-")), None)
    if subject:
        root = Path(*parts[:subject])
        below = parts[subject:]
        return (root / "derivatives" / pipeline).joinpath(*below), root / "derivatives" / pipeline, \
            "bids::" + "/".join(below + (scan.name,))
    pipeline_folder = scan.parent.parent / "derivatives" / pipeline
    return pipeline_folder, pipeline_folder, ""


def write_dataset_description(pipeline_folder) -> None:
    """The ``dataset_description.json`` a BIDS derivative needs, unless it is there."""
    path = Path(pipeline_folder) / "dataset_description.json"
    if path.is_file():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"Name": Path(pipeline_folder).name, "BIDSVersion": "1.8.0",
                                "DatasetType": "derivative", "GeneratedBy": [GENERATED_BY]}, indent=2))


def _plain(value):
    """*value* with numpy types replaced by what json writes."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_sidecar(image_path, content: dict) -> Path:
    """The JSON sidecar of an image."""
    path = Path(strip_extension(str(image_path)) + ".json")
    path.write_text(json.dumps(_plain(content), indent=2))
    return path


def save(volume: np.ndarray, affine: np.ndarray, path, scale: float = None, dtype=np.float32) -> Path:
    """Write a volume as NIfTI.

    Args:
        volume: The volume.
        affine: Its voxel-to-world matrix.
        path: The file (``.nii.gz`` or ``.nii``).
        scale: Store the volume as integers of *dtype* that are multiplied by
            *scale* when read: 1/255 with uint8 for a fraction, as CAT12
            stores its tissue maps.  None: as *dtype* itself.
        dtype: Type on disk.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if scale is None:
        data = np.asarray(volume, dtype=dtype)
    else:
        info = np.iinfo(dtype)
        data = np.clip(np.rint(np.asarray(volume, dtype=np.float64) / scale), info.min, info.max).astype(dtype)
    image = nib.Nifti1Image(data, affine)
    image.header.set_data_dtype(dtype)
    if scale is not None:
        image.header.set_slope_inter(scale, 0.0)
    image.header.set_xyzt_units("mm")
    image.set_qform(affine, code=4)          # MNI152
    image.set_sform(affine, code=4)
    nib.save(image, str(path))
    return path
