"""The artefacts t1sim makes in k-space: head movement and ringing.

The port of ``add_motion_artefacts`` and ``add_ringing`` of
``matlab/mri_simulate.m`` (documented in ``matlab/docs/artefacts.md``).  Both
run on the output grid, after the contrast of the image is final and before
the noise is added, and neither changes the ground truth.

**Movement.**  Every phase-encoding line of k-space is sampled at a different
time, so a line acquired after the head has moved belongs to a displaced
object, while the reconstruction treats all lines as one.  The lines are
therefore put into the order in which they are acquired (``ordering``) and
cut into ``blocks``, the temporal grid of the motion, and k-space is
assembled from the blocks of rigidly moved copies of the image.  The mismatch
between the blocks is the ringing and ghosting along the phase-encoding axis
that is typical for motion.

- A translation multiplies k-space by ``exp(-2 pi i k dr)`` (Fourier shift
  theorem): an exact phase ramp without any interpolation.
- A rotation rotates k-space itself and cannot be written as a phase, so the
  volume is resampled before its transform.  One resampling serves all blocks
  that share a rotation, thus the runtime follows the number of events and
  not the number of blocks.
- The pose time course is a number of ``events``, nods with pitch dominating:
  an excursion of one block plus a residual offset of 30% of it that lasts to
  the end of the scan.  Between the events a slow drift and two slow
  oscillations give every block its own pose (``continuous``); they are
  translations only, which are free.
- Only the differences between the poses make the artefact.  The mean of the
  poses weighted by the k-space energy of their blocks is where the complex
  image lies, so it is subtracted.
- The assembled k-space is no longer conjugate symmetric, so the image is the
  magnitude of a complex volume, as it is in a scanner.
- The magnitude is not linear, so the magnitude image does not lie where the
  complex image does: up to two voxels away for a severe motion.  Its
  displacement from the unmoved volume is therefore estimated and removed
  (``recentre``, see :func:`apply_motion`), so that the image stays where its
  label map is.

**Ringing.**  A gain on k-space, of two types that are different artefacts:
``notch`` damps or inverts a narrow band at ``|k| = k0`` (gain
``1 - strength``), which lays one spatial frequency over the whole image, the
pronounced regular ripples of mriaug; ``gibbs`` cuts k-space off at
``0.7 - 0.15 * strength`` of the Nyquist frequency, which is what a short
acquisition matrix does, and overshoots an edge by 9% whatever the fraction.
Along one axis (the phase encoding, by default) or along all three.  With
motion the gain is applied to the k-space the motion assembled, before the
magnitude is taken: the two share one acquisition, and that is not the same
as one after the other.

What is random (:func:`sample_motion`) is separate from what is not
(:func:`apply_motion`), so that the latter can be compared with MATLAB: fed
with the pose table of ``mri_simulate.m`` and with ``recentre=False``, it
gives the image of ``mri_simulate.m``.

Volumes are 3-D numpy arrays or torch tensors indexed ``[x, y, z]`` and come
back as what they were, a tensor on its device.  float64 is computed in
float64, everything else in float32 (complex64), which is what MPS has.  An
array is processed on the CPU; hand in a tensor to use a GPU.

What differs from the MATLAB version:

- Voxel axes and blocks are counted from 0: ``pe``, ``PhaseEncodingAxis``,
  ``EventBlocks`` and ``CentreBlock`` are one less than in MATLAB.
- The random numbers are numpy's and not MATLAB's.  What is drawn, in which
  order and from which distribution is the same.
- An option that does not exist is an error and is not ignored.
- The magnitude image is brought back onto the unmoved volume (``recentre``).
  ``mri_simulate.m`` stops at the mean pose, which left the image 0.4 to
  1.6 mm from its label map on a real head.

The rotation interpolates as ``rotate_volume`` of ``mri_simulate.m`` does,
with SPM's sinc interpolation (``spm_slice_vol`` with a hold of -5): a sinc
under a Hanning window of six points per axis, whose weights are normalised
to a sum of 1.  Points beyond the border of the volume are left out of the
kernel, and what is rotated in from more than 0.05 voxels outside is 0.
"""

from __future__ import annotations

import math

import numpy as np
import torch

__all__ = [
    "resolve_motion", "resolve_ringing", "sample_motion", "apply_motion",
    "add_motion", "add_ringing",
]

#: The motion draws from a stream of its own, so that the noise of a simulation
#: is the same with and without motion.  ``mri_simulate.m`` adds this number to
#: the seed; here it is the second word of the seed sequence.
MOTION_STREAM = 104729

#: Share of an excursion that stays as an offset: the head does not return
#: exactly to where it was.
RESIDUAL = 0.3

# Nodding, a pitch rotation around the left-right axis, dominates instructed
# head motion; the translation that goes with it is mostly anterior-posterior
# and inferior-superior.  Weights of (x, y, z).
_ROTATION_WEIGHTS = np.array([1.0, 0.35, 0.35])
_TRANSLATION_WEIGHTS = np.array([0.35, 1.0, 1.0])

# SPM's sinc interpolation: a hold of -5 is a kernel of 5 + 1 points per axis,
# which reaches from 2 points below the coordinate to 3 above.  A coordinate
# closer than _ON_GRID to a grid point takes that point alone, and one more
# than _FOV voxels outside the volume gives the background.
_POINTS = 6
_FIRST_TAP = -2
_ON_GRID = 1e-5
_FOV = 5e-2

# Voxels the resampling handles at a time.  Measured on a 256^3 volume: the CPU
# and MPS are both fastest around this size, and smaller is slower.
_CHUNK = 1 << 17

# Recentring: the image counts as aligned with the object once the displacement
# that is still found is below _ALIGNED voxels, and after _ROUNDS estimates the
# last one is removed without another look.
_ALIGNED = 1e-3
_ROUNDS = 8

_MOTION_OPTIONS = ("severity", "events", "translation", "rotation", "blocks",
                   "continuous", "pe", "ordering", "centre")
_RINGING_OPTIONS = ("strength", "type", "pe", "k0")


# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------

def _round(value: float) -> int:
    """MATLAB's ``round``: halves away from zero, where Python rounds them to even."""
    return int(math.copysign(math.floor(abs(value) + 0.5), value))


def _scalar(value, name: str) -> float:
    """A finite real number, or an error that names the option."""
    if isinstance(value, (str, bytes)) or np.ndim(value) != 0:
        raise ValueError(f"{name} has to be a scalar.")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} has to be a scalar.") from None
    if not math.isfinite(number):
        raise ValueError(f"{name} has to be finite.")
    return number


def _options(options, shortcut: str, known: tuple, name: str) -> dict:
    """The options as a dict: a scalar is the shortcut, ``None`` values count as not given."""
    if not isinstance(options, dict):
        if isinstance(options, (str, bytes)) or np.ndim(options) != 0:
            raise ValueError(f"{name} has to be a scalar {shortcut} or a dict.")
        return {shortcut: options}
    unknown = sorted(set(options) - set(known))
    if unknown:
        raise ValueError(f"{name} has no option {', '.join(map(repr, unknown))}; "
                         f"the options are {', '.join(known)}.")
    return {key: value for key, value in options.items() if value is not None}


def _resolved(options, key: str, name: str) -> dict:
    """The options, if they are what ``resolve_motion`` or ``resolve_ringing`` returns."""
    if not isinstance(options, dict) or key not in options:
        raise TypeError(f"{name} has to be what resolve_{name} returns, not {options!r}")
    return options


def _pe_axis(pe, affine, name: str) -> int:
    """The voxel axis a phase-encoding direction runs along.

    k-space is spanned by the voxel axes, so the result has to be one of
    them, but a world axis describes the same anatomical direction whatever
    orientation the image is stored in.  An oblique image is served by the
    closest voxel axis.

    Args:
        pe: 0, 1 or 2 for a voxel axis, or 'x' (left-right), 'y'
            (anterior-posterior) or 'z' (inferior-superior) for a world axis.
        affine: (4, 4) voxel-to-world matrix of the grid.
        name: Name of the option, for the error message.

    Returns:
        The voxel axis, 0, 1 or 2.
    """
    error = ValueError(f"{name} has to be 0, 1, 2 or 'x', 'y', 'z'.")
    if isinstance(pe, str):
        if len(pe) != 1 or pe.lower() not in "xyz":
            raise error
        directions = np.asarray(affine, dtype=np.float64)[:3, :3]
        # world direction of every voxel axis, and how much of the wanted one it has
        directions = directions / np.sqrt((directions ** 2).sum(axis=0))
        return int(np.argmax(np.abs(directions["xyz".index(pe.lower())])))
    if isinstance(pe, (int, np.integer)) and not isinstance(pe, bool) and 0 <= pe <= 2:
        return int(pe)
    raise error


def resolve_motion(options, affine, shape) -> dict | None:
    """The full set of movement parameters, checked.

    A scalar is the severity: 1, 2 and 3 are mild, moderate and severe, and
    values in between and beyond are allowed.  The number of events and both
    amplitudes grow with ``severity ** 1.5``, which maps 1, 2 and 3 onto 1, 3
    and 5 events of 1, 2.8 and 5.2 mm and degrees, the ranges of the severity
    classes of datasets like MR-ART.  A dict overrides single values:

    - ``severity``: the base of all defaults below (1 if it is missing: a
      dict asks for motion);
    - ``events``: number of motion events (``round(severity ** 1.5)``, at
      least 1);
    - ``translation``, ``rotation``: maximum amplitude of an event in mm and
      in degrees (``severity ** 1.5``);
    - ``blocks``: number of k-space blocks, the temporal resolution of the
      motion and so the duration of one excursion (32, i.e. 5-10 s of a 5 min
      MPRAGE, the time scale of a nod);
    - ``continuous``: amplitude of the drift and tremor between the events,
      as a fraction of ``translation`` (0.2; 0 leaves the pose piecewise
      constant);
    - ``pe``: phase-encoding direction, a voxel axis (0, 1, 2) or a world
      axis ('x', 'y', 'z'), which is mapped to the closest voxel axis ('y',
      the usual in-plane direction of a sagittal MPRAGE);
    - ``ordering``: 'linear' fills k-space from -kmax to +kmax, so its centre
      is sampled in the middle of the scan; 'centric' starts at the centre
      ('linear');
    - ``centre``: place the first event in the block that samples the centre
      of k-space (True).  The artefact depends far more on whether the centre
      is hit than on the amplitude, so this keeps the severities comparable
      between images.

    Args:
        options: The severity, a dict of the options above, or None.
        affine: (4, 4) voxel-to-world matrix of the output grid (nibabel's).
        shape: Shape of the volume the motion is made for.

    Returns:
        A dict with all options above, ``pe`` as a voxel axis (from 0),
        ``blocks`` limited to the number of phase-encoding lines, and that
        number as ``lines``; or None when no motion is asked for (None, or a
        severity of 0 or less).
    """
    if options is None:
        return None
    options = _options(options, "severity", _MOTION_OPTIONS, "motion")
    severity = _scalar(options.get("severity", 1), "motion['severity']")
    if severity <= 0:
        return None

    amplitude = severity ** 1.5
    motion = {
        "severity": severity,
        "events": max(1, _round(amplitude)),
        "translation": amplitude,   # mm
        "rotation": amplitude,      # degrees
        "blocks": 32,
        "continuous": 0.2,
        "pe": "y",
        "ordering": "linear",
        "centre": True,
    }
    motion.update(options, severity=severity)

    for key in ("events", "translation", "rotation", "blocks", "continuous"):
        motion[key] = _scalar(motion[key], f"motion['{key}']")
    motion["pe"] = _pe_axis(motion["pe"], affine, "motion['pe']")
    motion["ordering"] = str(motion["ordering"]).lower()
    if motion["ordering"] not in ("linear", "centric"):
        raise ValueError("motion['ordering'] has to be 'linear' or 'centric'.")
    if motion["events"] < 1:
        raise ValueError("motion['events'] has to be at least 1.")
    if motion["blocks"] < 2:
        raise ValueError("motion['blocks'] has to be at least 2, otherwise all of k-space "
                         "shares one pose and no artefact is created.")
    if motion["translation"] < 0 or motion["rotation"] < 0:
        raise ValueError("motion['translation'] and motion['rotation'] must not be negative.")
    if motion["translation"] == 0 and motion["rotation"] == 0:
        raise ValueError("motion['translation'] and motion['rotation'] are both 0, thus "
                         "nothing would move.")
    if motion["continuous"] < 0:
        raise ValueError("motion['continuous'] must not be negative.")

    if len(shape) != 3:
        raise ValueError(f"the volume must be 3-D, not {tuple(shape)}")
    lines = int(shape[motion["pe"]])
    motion["events"] = _round(motion["events"])
    # more blocks than lines would leave empty blocks behind
    motion["blocks"] = max(2, min(_round(motion["blocks"]), lines))
    motion["centre"] = bool(motion["centre"])
    motion["lines"] = lines
    return motion


def resolve_ringing(options, affine) -> dict | None:
    """The shape of the k-space gain that makes the ringing, checked.

    A scalar is the strength: 1, 2 and 3 are mild, moderate and severe, any
    positive value is allowed.  A dict selects and configures the type:

    - ``strength``: as above (1 if it is missing);
    - ``type``: 'notch' (default) changes the gain of a narrow band at
      ``|k| = k0`` to ``1 - strength``, i.e. removes the band for 1 and
      inverts it beyond.  It looks like ringing and is not what a scanner
      does.  'gibbs' keeps ``max(0.15, 0.7 - 0.15 * strength)`` of k-space,
      as a short acquisition matrix does: the overshoot is 9% whatever the
      fraction, only the period of the ripples grows as it shrinks (about
      ``2 / fraction`` voxels), and the image is blurred by as much;
    - ``pe``: the axis the ringing runs along, as for the motion, or 'all'
      ('y').  One axis is a flat band or a shorter phase encoding, 'all' a
      spherical shell or a lower resolution in every direction, which for
      'gibbs' reads as smoothing rather than as ringing;
    - ``k0``: centre of the band of 'notch', in units of the Nyquist
      frequency (0.6, where mriaug puts it).

    Args:
        options: The strength, a dict of the options above, or None.
        affine: (4, 4) voxel-to-world matrix of the output grid (nibabel's).

    Returns:
        A dict with ``strength``, ``type``, ``axes`` (voxel axes, from 0) and
        ``fraction`` for 'gibbs' or ``k0`` and ``gain`` for 'notch'; or None
        when no ringing is asked for (None, or a strength of 0).
    """
    if options is None:
        return None
    options = _options(options, "strength", _RINGING_OPTIONS, "ringing")
    strength = _scalar(options.get("strength", 1), "ringing['strength']")
    if strength < 0:
        raise ValueError("ringing['strength'] must not be negative.")
    if strength == 0:
        return None

    kind = str(options.get("type", "notch")).lower()
    if kind not in ("notch", "gibbs"):
        raise ValueError("ringing['type'] has to be 'notch' or 'gibbs'.")
    pe = options.get("pe", "y")
    if isinstance(pe, str) and pe.lower() == "all":
        axes = [0, 1, 2]
    else:
        axes = [_pe_axis(pe, affine, "ringing['pe']")]

    ringing = {"strength": strength, "type": kind, "axes": axes}
    if kind == "gibbs":
        ringing["fraction"] = max(0.15, 0.7 - 0.15 * strength)
    else:
        k0 = _scalar(options.get("k0", 0.6), "ringing['k0']")
        if not 0 < k0 < 1:
            raise ValueError("ringing['k0'] has to be between 0 and 1, in units of the "
                             "Nyquist frequency.")
        ringing["k0"] = k0
        # 0 removes the band, a negative gain inverts it, which is what makes
        # the ripples strong enough to see
        ringing["gain"] = 1 - strength
    return ringing


# ---------------------------------------------------------------------------
# volumes and k-space
# ---------------------------------------------------------------------------

def _working(volume):
    """The volume as the tensor the computation runs on, and the way back.

    Returns:
        The tensor, float64 for a float64 input and float32 otherwise, and a
        function that turns a result into what the input was: an array for an
        array, a tensor on the same device for a tensor.
    """
    if isinstance(volume, torch.Tensor):
        dtype = torch.float64 if volume.dtype == torch.float64 else torch.float32
        work = volume.detach().to(dtype)

        def restore(result):
            return result.contiguous()
    else:
        array = np.asarray(volume)
        double = array.dtype.kind == "f" and array.dtype.itemsize == 8
        array = np.ascontiguousarray(array, dtype=np.float64 if double else np.float32)
        if not array.flags.writeable:
            array = array.copy()        # torch does not take a read-only array
        work = torch.from_numpy(array)

        def restore(result):
            return result.contiguous().cpu().numpy()
    if work.ndim != 3:
        raise ValueError(f"the volume must be 3-D, not {tuple(work.shape)}")
    return work, restore


def _frequencies(size: int) -> np.ndarray:
    """Signed frequency index of every element of an FFT axis: 0, 1, ..., -2, -1."""
    return (np.arange(size) + size // 2) % size - size // 2


def _along(values: np.ndarray, *axes: int) -> np.ndarray:
    """Make ``values`` 3-D, with its dimensions on ``axes`` and a length of 1 elsewhere."""
    values = np.transpose(values, np.argsort(axes))
    shape = [1, 1, 1]
    for axis, length in zip(sorted(axes), values.shape):
        shape[axis] = length
    return values.reshape(shape)


def _window(shape, ringing: dict, like: torch.Tensor) -> list:
    """The factors of the gain that is multiplied onto k-space to make the ringing.

    'gibbs' is a plain rectangle per axis and not a smooth filter, which
    would damp the very ringing that is wanted.  'notch' is a band whose gain
    is changed: a flat one for a single axis, a spherical shell for all three.
    Every shape depends on ``|k|`` alone, and the one Nyquist bin without a
    counterpart is never inside it, so the gain is symmetric and the
    reconstruction of a real image stays real.

    Args:
        shape: Shape of the volume.
        ringing: The options :func:`resolve_ringing` returns.
        like: A real tensor with the dtype and the device of the result.

    Returns:
        Real tensors that broadcast against k-space.  The gain is separable
        but for the shell, so it is kept as its factors.
    """
    # distance from the centre of k-space along each axis, 1 at Nyquist; in
    # float64, so that a bin on the edge of the band falls where it does in MATLAB
    distance = [_along(np.abs(_frequencies(size)) / (size / 2), axis)
                for axis, size in enumerate(shape)]
    axes = ringing["axes"]
    if ringing["type"] == "gibbs":
        factors = [distance[axis] <= ringing["fraction"] for axis in axes]
    else:
        low, high = 0.95 * ringing["k0"], 1.05 * ringing["k0"]
        if len(axes) == 1:
            radius = distance[axes[0]]
        else:
            radius = np.sqrt(distance[0] ** 2 + distance[1] ** 2 + distance[2] ** 2)
        factors = [1 + (ringing["gain"] - 1) * ((radius > low) & (radius < high))]
    return [torch.as_tensor(factor, dtype=like.dtype, device=like.device) for factor in factors]


def _acquisition(lines: int, blocks: int, ordering: str):
    """When every phase-encoding line is acquired.

    Args:
        lines: Number of phase-encoding lines.
        blocks: Number of blocks the acquisition is cut into.
        ordering: 'linear' runs from -kmax to +kmax, 'centric' starts at the
            centre and works outwards, alternating between +k and -k.

    Returns:
        ``line``, the index into the FFT axis of every acquisition step,
        ``block``, the block of every step, and the block that samples the
        centre of k-space.
    """
    k = np.arange(lines) - lines // 2
    if ordering == "linear":
        order = np.arange(lines)
    else:
        # the offset for negative k only breaks the tie between +k and -k
        order = np.argsort(np.abs(k) + 0.25 * (k < 0), kind="stable")
    # the FFT puts k = 0 first and the negative frequencies into the upper half
    line = k[order] % lines
    # round(linspace(0, lines, blocks + 1)) of MATLAB, in integers
    edges = (2 * np.arange(blocks + 1) * lines + blocks) // (2 * blocks)
    block = np.repeat(np.arange(blocks), np.diff(edges))
    centre = int(block[np.flatnonzero(order == lines // 2)[0]])
    return line, block, centre


def _block_weights(work: torch.Tensor, pe: int, line: np.ndarray, block: np.ndarray,
                   blocks: int) -> np.ndarray:
    """Share of the k-space energy every block samples.

    Only a transform along the phase-encoding axis is needed: by Parseval the
    energy of a line does not depend on how it is distributed along the two
    other axes.
    """
    other = tuple(axis for axis in range(3) if axis != pe)
    energy = torch.fft.fft(work, dim=pe).abs().square().sum(dim=other)
    energy = energy.cpu().numpy().astype(np.float64)
    weights = np.bincount(block, weights=energy[line], minlength=blocks)
    total = weights.sum()
    if not total > 0:
        return np.full(blocks, 1.0 / blocks)    # an empty volume has no pose to prefer
    return weights / total


# ---------------------------------------------------------------------------
# the rotation
# ---------------------------------------------------------------------------

def _rotation(degrees, affine: np.ndarray, shape) -> np.ndarray:
    """Where in the volume every voxel of the rotated volume comes from.

    The rotation is around the centre of the volume and defined in world
    coordinates, so that the same angles describe the same anatomical
    rotation whatever orientation the image is stored in.  The angles are
    those of ``spm_matrix``: R = Rx Ry Rz.

    Args:
        degrees: Rotation around the world axes x, y and z.
        affine: (4, 4) voxel-to-world matrix of the grid.
        shape: Shape of the volume.

    Returns:
        (4, 4) matrix from a voxel of the rotated volume to its position in
        the volume, both counted from 0.
    """
    rx, ry, rz = np.deg2rad(np.asarray(degrees, dtype=np.float64))
    about_x = np.array([[1, 0, 0], [0, np.cos(rx), np.sin(rx)], [0, -np.sin(rx), np.cos(rx)]])
    about_y = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    about_z = np.array([[np.cos(rz), np.sin(rz), 0], [-np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    rotation = about_x @ about_y @ about_z

    centre = affine[:3, :3] @ ((np.asarray(shape, dtype=np.float64) - 1) / 2) + affine[:3, 3]
    world = np.eye(4)
    world[:3, :3] = rotation
    world[:3, 3] = centre - rotation @ centre
    # the same transform between voxels moves the object; resampling asks for
    # the way back
    return np.linalg.inv(np.linalg.solve(affine, world @ affine))


def _sinc_weights(index: torch.Tensor, shift: torch.Tensor, size: int):
    """SPM's interpolation kernel along one axis.

    Args:
        index: (N,) integer position of the voxels along this axis.
        shift: (N,) how far from it their value is taken: the coordinate is
            ``index + shift``.  Kept apart because the shift is small, so its
            fraction survives float32, which that of the coordinate would not.
        size: Length of the axis.

    Returns:
        ``first``, (N,) position of the first of the six points in a volume
        padded by three voxels, ``weights``, (N, 6), and ``inside``, (N,),
        False where the coordinate lies outside the field of view.
    """
    below = torch.floor(shift)
    nearest = torch.floor(shift + 0.5)
    taps = torch.arange(_FIRST_TAP, _FIRST_TAP + _POINTS, device=shift.device)
    x = math.pi * ((shift - below)[:, None] - taps)
    weights = torch.sin(x) / x * (0.5 * (1 + torch.cos(2 * x / _POINTS)))
    # on a grid point the kernel is that point alone (and x is 0 there)
    on_grid = ((shift - nearest).abs() < _ON_GRID)[:, None]
    alone = (taps == (nearest - below)[:, None]).to(shift.dtype)
    weights = torch.where(on_grid, alone, weights)

    # points beyond the border are left out, and the rest shares the weight
    position = index + below.long()
    points = position[:, None] + taps
    weights = weights * ((points >= 0) & (points < size))
    total = weights.sum(dim=1, keepdim=True)
    weights = weights / torch.where(total == 0, 1, total)

    coordinate = index.to(shift.dtype) + shift
    inside = (coordinate >= -_FOV) & (coordinate <= size - 1 + _FOV)
    # outside the field of view the result is discarded, the points only have
    # to stay within the padded volume
    first = position.clamp(-1, size - 1) + (_FIRST_TAP + _POINTS // 2)
    return first, weights, inside


def _resample(work: torch.Tensor, mat: np.ndarray) -> torch.Tensor:
    """Sinc interpolation of a volume at the positions an affine map gives.

    ``spm_slice_vol`` with a hold of -5, slice by slice, as ``rotate_volume``
    of ``mri_simulate.m`` uses it: the kernel of :func:`_sinc_weights` along
    each axis, 216 points per voxel, and 0 for what comes from outside the
    field of view.

    Args:
        work: (X, Y, Z) real tensor.
        mat: (4, 4) matrix from a voxel of the result to its position in
            ``work``, both counted from 0.

    Returns:
        The resampled volume, same shape, dtype and device.
    """
    nx, ny, nz = work.shape
    dtype, device = work.dtype, work.device
    # Zeros around the volume: the six points of a voxel are then at fixed
    # offsets from its first one, also at the border, where they have no weight
    pad = _POINTS // 2
    flat = torch.nn.functional.pad(work, (pad,) * 6).reshape(-1)
    step_y = nz + 2 * pad
    step_x = (ny + 2 * pad) * step_y
    taps = torch.arange(_POINTS, device=device)
    plane = (taps[:, None] * step_y + taps).reshape(1, -1)   # the 36 points of a y-z plane

    # the map as a displacement, see _sinc_weights
    linear = torch.as_tensor(mat[:3, :3] - np.eye(3), dtype=dtype, device=device)
    offset = torch.as_tensor(mat[:3, 3].copy(), dtype=dtype, device=device)

    result = torch.empty(nx * ny * nz, dtype=dtype, device=device)
    # in chunks, because every voxel needs 36 indices and values at a time
    for start in range(0, result.numel(), _CHUNK):
        voxel = torch.arange(start, min(start + _CHUNK, result.numel()), device=device)
        index = (torch.div(voxel, ny * nz, rounding_mode="floor"),
                 torch.div(voxel, nz, rounding_mode="floor") % ny,
                 voxel % nz)
        real = [i.to(dtype) for i in index]

        first, weights, inside = [], [], None
        for axis in range(3):
            shift = (linear[axis, 0] * real[0] + linear[axis, 1] * real[1]
                     + linear[axis, 2] * real[2] + offset[axis])
            begin, weight, within = _sinc_weights(index[axis], shift, work.shape[axis])
            first.append(begin)
            weights.append(weight)
            inside = within if inside is None else inside & within

        points = ((first[0] * step_x + first[1] * step_y + first[2])[:, None] + plane).reshape(-1)
        in_plane = (weights[1][:, :, None] * weights[2][:, None, :]).reshape(len(voxel), -1)
        value = torch.zeros(len(voxel), dtype=dtype, device=device)
        for tap in range(_POINTS):
            # the plane one step further along x is the same indices in a
            # volume that starts one step later
            taken = flat[tap * step_x:][points].view(len(voxel), -1)
            value += weights[0][:, tap] * (taken * in_plane).sum(dim=1)
        result[start:start + len(voxel)] = torch.where(inside, value, 0)
    return result.view(nx, ny, nz)


# ---------------------------------------------------------------------------
# the place of the image
# ---------------------------------------------------------------------------

def _shift(kspace: torch.Tensor, voxels) -> None:
    """Move the image of a k-space by a translation, in place: an exact phase ramp."""
    for axis, size in enumerate(kspace.shape):
        ramp = np.exp(-2j * np.pi * _frequencies(size) / size * voxels[axis])
        kspace *= torch.as_tensor(_along(ramp, axis), dtype=kspace.dtype, device=kspace.device)


def _displacement(image: torch.Tensor, reference: torch.Tensor, coarse: bool) -> np.ndarray:
    """How far the object in an image lies from where the unmoved volume has it.

    The translation that fits the image to the volume best is the lag at
    which their cross-correlation is largest (least squares, and for a
    translation the same as the largest correlation coefficient).  The volume
    is close to 0 outside the object, so ghosts in the background of the
    image do not take part.  The peak is looked for among the lags of whole
    voxels first, and from there by a Newton step on the correlation between
    the voxels, whose slope and curvature at lag 0 are moments of the cross
    power spectrum: they are sums and not differences of nearly equal
    correlations, which float32 would not resolve.

    Args:
        image: (X, Y, Z) real tensor, the magnitude image.
        reference: Conjugate of the real Fourier transform (``rfftn``) of the
            unmoved volume.
        coarse: Look among the whole voxels; False when the peak is known to
            be the one at lag 0.

    Returns:
        (3,) float64, the displacement of the image in voxels: whole voxels
        if the peak is not at lag 0, and the Newton step otherwise, which is
        0 when the image is where the volume is.
    """
    shape = np.array(image.shape)
    cross = torch.fft.rfftn(image)
    cross *= reference
    if coarse:
        peak = int(torch.argmax(torch.fft.irfftn(cross, s=image.shape)))
        lag = (np.array(np.unravel_index(peak, image.shape)) + shape // 2) % shape - shape // 2
        if lag.any():           # the upper half of an axis are the negative lags
            return lag.astype(np.float64)

    # c(d) = sum_k C(k) exp(i w d) / N, so its slope at 0 is -sum w Im(C) / N
    # and its curvature -sum w w' Re(C) / N.  Summed along one axis first, and
    # along the last one with the weight of the half that rfftn leaves out.
    omega = [2 * np.pi * _frequencies(size) / size for size in shape]
    omega[2] = np.abs(omega[2][:cross.shape[2]])
    twice = np.full(cross.shape[2], 2.0)
    twice[0] = 1                                    # k = 0 and Nyquist are their own mirror
    if shape[2] % 2 == 0:
        twice[-1] = 1
    planes = {(0, 1): torch.matmul(cross, torch.as_tensor(twice, dtype=cross.dtype, device=cross.device)),
              (0, 2): cross.sum(dim=1), (1, 2): cross.sum(dim=0)}
    planes = {axes: plane.cpu().numpy().astype(np.complex128) for axes, plane in planes.items()}
    planes[0, 2] *= twice
    planes[1, 2] *= twice
    lines = [planes[0, 1].sum(axis=1), planes[0, 1].sum(axis=0), planes[0, 2].sum(axis=0)]
    slope = np.array([-(omega[axis] * lines[axis].imag).sum() for axis in range(3)])
    curvature = np.diag([-(omega[axis] ** 2 * lines[axis].real).sum() for axis in range(3)])
    for (a, b), plane in planes.items():
        curvature[a, b] = curvature[b, a] = -(omega[a] @ plane.real @ omega[b])
    if not np.linalg.eigvalsh(curvature).max() < 0:
        return np.zeros(3)      # no peak here, e.g. an empty volume
    # the peak is within the voxel; a step beyond it is not one to trust
    return np.clip(-np.linalg.solve(curvature, slope), -1, 1)


def _recentre(kspace: torch.Tensor, work: torch.Tensor):
    """The magnitude image of a k-space, moved to where the unmoved volume is.

    The mean pose that :func:`apply_motion` subtracts is where the complex
    image lies.  The magnitude image lies elsewhere, because the magnitude is
    not linear, and the label map belongs to the magnitude image.  Its
    displacement is estimated (:func:`_displacement`) and removed with a
    phase ramp on k-space, which is exact, and this is repeated, since the
    magnitude of a shifted complex image is not the shifted magnitude, until
    nothing is left to remove.

    Args:
        kspace: (X, Y, Z) complex tensor, the assembled k-space.  Changed in
            place.
        work: (X, Y, Z) real tensor, the unmoved volume.

    Returns:
        The magnitude image, and the displacement in voxels that was removed
        from it, (3,) float64.
    """
    reference = torch.fft.rfftn(work).conj()
    removed = np.zeros(3)
    coarse = True
    for _ in range(_ROUNDS):
        image = torch.fft.ifftn(kspace).abs()
        found = _displacement(image, reference, coarse)
        if np.abs(found).max() < _ALIGNED:
            return image, removed
        coarse = coarse and np.abs(found).max() >= 1    # whole voxels until lag 0 is the peak
        removed += found
        _shift(kspace, -found)
    return torch.fft.ifftn(kspace).abs(), removed


# ---------------------------------------------------------------------------
# movement
# ---------------------------------------------------------------------------

def sample_motion(motion: dict, seed):
    """Draw the pose time course of a scan: the random part of the motion.

    Every event is a nod of a random amplitude between 30% and 100% of the
    maximum, in a random direction in which pitch and the translations that
    go with it dominate: an excursion during its block and
    :data:`RESIDUAL` of it from then on.  Between the events, and only as a
    translation, a drift of half the continuous amplitude over the scan and
    two oscillations of 1.5 to 2.5 and of 4 to 6 periods per scan.

    Args:
        motion: The options :func:`resolve_motion` returns.
        seed: Seed of the simulation; the motion gets a stream of its own
            from it (see :data:`MOTION_STREAM`).  None draws another motion
            every time, and a ``numpy.random.Generator`` is used as it is.

    Returns:
        ``pose``, (blocks, 6) float64, one row per k-space block in
        acquisition order: the translation in mm and the rotation in degrees
        around x, y and z, in world coordinates; and ``events``, the block of
        every event (from 0).  The pose is as drawn: :func:`apply_motion`
        subtracts its mean.
    """
    motion = _resolved(motion, "lines", "motion")
    if hasattr(seed, "standard_normal"):
        rng = seed
    elif seed is None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng([int(seed), MOTION_STREAM])

    blocks = motion["blocks"]
    events = np.asarray(rng.integers(0, blocks, size=motion["events"]))
    if motion["centre"]:
        events[0] = _acquisition(motion["lines"], blocks, motion["ordering"])[2]

    eps = np.finfo(np.float64).eps
    pose = np.zeros((blocks, 6))
    for at in events:
        amplitude = 0.3 + 0.7 * rng.random()
        turn = rng.standard_normal(3) * _ROTATION_WEIGHTS
        move = rng.standard_normal(3) * _TRANSLATION_WEIGHTS
        nod = np.concatenate([
            amplitude * motion["translation"] * move / max(np.linalg.norm(move), eps),
            amplitude * motion["rotation"] * turn / max(np.linalg.norm(turn), eps)])
        # the excursion of the nod itself, and the offset that stays
        pose[at] += nod
        pose[at + 1:] += RESIDUAL * nod

    # A pose that is piecewise constant with a handful of steps gives a handful
    # of discontinuities in k-space and so ringing along the edges alone; a
    # real motion corrupted image carries the dense ripples of many small ones.
    if motion["continuous"] > 0:
        time = np.arange(blocks) / blocks
        amplitude = motion["continuous"] * motion["translation"]
        for axis in range(3):
            slow = 0.3 * np.sin(2 * np.pi * (1.5 + rng.random()) * time + 2 * np.pi * rng.random())
            fast = 0.2 * np.sin(2 * np.pi * (4 + 2 * rng.random()) * time + 2 * np.pi * rng.random())
            pose[:, axis] += amplitude * (0.5 * time + slow + fast)
    return pose, [int(at) for at in events]


def apply_motion(volume, pose, motion: dict, affine, ringing: dict | None = None, *,
                 events=None, subtract_mean: bool = True, recentre: bool = True):
    """Corrupt a volume with a given pose time course: the motion without its chance.

    k-space is assembled from the phase-encoding lines of moved copies of the
    volume, every block of lines from the copy in the pose of that block, and
    the magnitude of its inverse transform is the image.

    Args:
        volume: (X, Y, Z) array or tensor, the image on the output grid
            before the noise.
        pose: (blocks, 6), one row per k-space block in acquisition order:
            translation in mm and rotation in degrees around x, y and z, in
            world coordinates (:func:`sample_motion`).
        motion: The options :func:`resolve_motion` returns, for a volume of
            this shape.
        affine: (4, 4) voxel-to-world matrix of the output grid (nibabel's).
        ringing: The options :func:`resolve_ringing` returns, or None.  The
            gain is applied to the assembled k-space before the magnitude is
            taken, because ringing and motion belong to one acquisition.
        events: The blocks of the events, for the record; they are not used.
        subtract_mean: Subtract the mean of the poses weighted by the k-space
            energy of their blocks, which is where the complex image lies.
            False takes the pose as it is, e.g. the ``Pose`` of a sidecar.
        recentre: Bring the magnitude image back onto the unmoved object.
            The mean pose does not do that: the magnitude is not linear, and
            on a head-sized object the image was left 0.5, 1.4 and 1.9
            voxels away for a severity of 1, 2 and 3 (mean of four seeds; on
            a real head at 1 mm: 0.4, 0.9 and 1.6 mm).  A label map that far
            from its image counts as an error of the segmentation and is
            noise in training data.  The displacement is estimated from the
            cross-correlation with the unmoved volume and removed by a phase
            ramp on the assembled k-space, which is exact, until less than
            0.001 voxels are found; afterwards no shift fits better than
            none (0.005, 0.03 and 0.04 voxels left by an independent
            measure).  Only the translation is removed.  False gives the
            image as ``mri_simulate.m`` makes it.

    Returns:
        The image with the movement artefacts, the type of ``volume``, and
        the realised motion as ``mri_simulate.m`` writes it to the sidecar of
        the image, which can be dumped as JSON: ``Severity``, ``Events``,
        ``EventBlocks``, ``MaxTranslation``, ``MaxRotation``, ``Blocks``,
        ``PhaseEncodingAxis``, ``Ordering``, ``CentreBlock``, ``Continuous``,
        ``Recentred`` (the displacement in mm of the world that the
        recentring removed; None without it) and ``Pose``, the movement
        relative to the image that is returned, i.e. after the subtraction
        and the recentring.  Axes and blocks count from 0.
    """
    motion = _resolved(motion, "lines", "motion")
    if ringing is not None:
        ringing = _resolved(ringing, "axes", "ringing")
    work, restore = _working(volume)
    affine = np.asarray(affine, dtype=np.float64)
    pose = np.array(pose, dtype=np.float64)
    pe, blocks = motion["pe"], motion["blocks"]
    lines = work.shape[pe]
    if lines != motion["lines"]:
        raise ValueError(f"the motion was resolved for {motion['lines']} phase-encoding "
                         f"lines, the volume has {lines}")
    if pose.shape != (blocks, 6):
        raise ValueError(f"the pose must have the shape ({blocks}, 6), not {pose.shape}")

    line, block, centre = _acquisition(lines, blocks, motion["ordering"])
    if subtract_mean:
        pose = pose - _block_weights(work, pe, line, block, blocks) @ pose

    realised = {
        "Severity": motion["severity"],
        "Events": motion["events"],
        "EventBlocks": None if events is None else [int(at) for at in events],
        "MaxTranslation": motion["translation"],
        "MaxRotation": motion["rotation"],
        "Blocks": blocks,
        "PhaseEncodingAxis": pe,
        "Ordering": motion["ordering"],
        "CentreBlock": centre,
        "Continuous": motion["continuous"],
        "Recentred": [0.0, 0.0, 0.0] if recentre else None,
        "Pose": pose.tolist(),
    }

    # every block ended up with the same pose: nothing to splice
    if not pose.any():
        return restore(_ring(work, ringing)), realised

    # Phase ramp of the translation of every line, as its factors along the
    # axes: each is 2-D, an axis by the lines, because every line has the
    # translation of its block.  The full ramp is never built.
    block_of_line = np.empty(lines, dtype=np.intp)
    block_of_line[line] = block
    shift = np.linalg.solve(affine[:3, :3], pose[:, :3].T).T[block_of_line]   # in voxels
    complex_dtype = torch.complex128 if work.dtype == torch.float64 else torch.complex64
    ramp = []
    for axis in range(3):
        k = _frequencies(work.shape[axis]) / work.shape[axis]    # cycles per voxel
        if axis == pe:
            factor = _along(np.exp(-2j * np.pi * k * shift[:, axis]), pe)
        else:
            factor = _along(np.exp(-2j * np.pi * np.outer(k, shift[:, axis])), axis, pe)
        ramp.append(factor)
    along_pe = ramp.pop(pe)
    ramp = [torch.as_tensor(factor, dtype=complex_dtype, device=work.device) for factor in ramp]

    # One resampling per distinct rotation, and the translations inside it as
    # phase ramps: a rotation costs a resampling of the whole volume while a
    # translation is free, so the runtime follows the number of events.
    groups = {}
    for index, rotation in enumerate(pose[:, 3:] + 0.0):    # + 0.0: -0 is no rotation of its own
        groups.setdefault(tuple(rotation.tolist()), []).append(index)

    kspace = None
    for rotation, members in groups.items():
        moved = _resample(work, _rotation(rotation, affine, work.shape)) if any(rotation) else work
        part = torch.fft.fftn(moved)
        del moved
        # only the lines of these blocks are taken, the others get a gain of 0
        select = along_pe * _along(np.isin(block_of_line, members), pe)
        part *= torch.as_tensor(select, dtype=complex_dtype, device=work.device)
        for factor in ramp:
            part *= factor
        if kspace is None:
            kspace = part
        else:
            kspace += part
        del part

    if ringing is not None:
        for factor in _window(work.shape, ringing, work):
            kspace *= factor
    if not recentre:
        return restore(torch.fft.ifftn(kspace).abs()), realised

    image, removed = _recentre(kspace, work)
    removed = affine[:3, :3] @ removed          # in mm
    pose[:, :3] -= removed
    realised.update(Recentred=removed.tolist(), Pose=pose.tolist())
    return restore(image), realised


def add_motion(volume, motion: dict | None, affine, seed, ringing: dict | None = None, *,
               recentre: bool = True):
    """Add movement artefacts, and the ringing that shares their k-space.

    :func:`sample_motion` and :func:`apply_motion` in one call.

    Args:
        volume: (X, Y, Z) array or tensor, the image on the output grid
            before the noise.
        motion: The options :func:`resolve_motion` returns.  None adds the
            ringing alone.
        affine: (4, 4) voxel-to-world matrix of the output grid (nibabel's).
        seed: Seed of the simulation (see :func:`sample_motion`).
        ringing: The options :func:`resolve_ringing` returns, or None.
        recentre: Keep the image on the unmoved object (see
            :func:`apply_motion`).

    Returns:
        The image with the artefacts, the type of ``volume``, and the
        realised motion (see :func:`apply_motion`), None without motion.
    """
    if motion is None:
        return add_ringing(volume, ringing), None
    pose, events = sample_motion(motion, seed)
    return apply_motion(volume, pose, motion, affine, ringing, events=events, recentre=recentre)


# ---------------------------------------------------------------------------
# ringing
# ---------------------------------------------------------------------------

def _ring(work: torch.Tensor, ringing: dict | None) -> torch.Tensor:
    """The ringing of a volume that did not move."""
    if ringing is None:
        return work.clone()
    kspace = torch.fft.fftn(work)
    for factor in _window(work.shape, _resolved(ringing, "axes", "ringing"), work):
        kspace *= factor
    # the gain is symmetric, so the result is real but for rounding
    return torch.fft.ifftn(kspace).real


def add_ringing(volume, ringing: dict | None):
    """Add the ringing alone, to an image without movement artefacts.

    With motion the gain belongs onto the k-space the motion assembled:
    :func:`apply_motion` takes the ringing along.

    Args:
        volume: (X, Y, Z) array or tensor, the image on the output grid
            before the noise.
        ringing: The options :func:`resolve_ringing` returns; None returns a
            copy of the volume.

    Returns:
        The image with the ringing, the type of ``volume``.  It is the real
        part and not the magnitude, so an undershoot below 0 stays one.
    """
    work, restore = _working(volume)
    return restore(_ring(work, ringing))
