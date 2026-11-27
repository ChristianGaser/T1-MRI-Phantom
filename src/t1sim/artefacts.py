"""What an acquisition does to the image: bias field, resolution, contrast, noise.

Every function takes the image and explicit parameters and returns a new
image; the ground truth is never touched.  Movement artefacts and ringing,
which are made in k-space, are in :mod:`t1sim.kspace`.

What is and is not a port of ``matlab/mri_simulate.m``:

- The bias field is: a simulated field is the very field MATLAB makes for a
  seed.
- The noise is the same model without its upper limit (see
  :func:`add_noise`).
- The contrast change is the same power law, but applied to the intensities
  of the tissues before the image is mixed from them and before the bias
  field (see :func:`change_contrast`); MATLAB applies it to the finished
  image.
- The resolution is not.  ``mri_simulate.m`` interpolates the image at the
  new voxel centres, which aliases, and interpolates the label linearly,
  which samples it at a point.  Here the image is acquired: its k-space is
  cut at the Nyquist frequency of the new voxel size, which is what a smaller
  acquisition matrix does, and the tissue fractions are averaged over the box
  of every new voxel, which is what a partial volume is.
"""

from __future__ import annotations

import itertools
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import fft as sp_fft
from scipy import ndimage

from ._volume import smooth, voxel_size

#: Smoothness of a simulated bias field: FWHM in mm of its final smoothing.
BIAS_FWHM_MM = 15.0


# ------------------------------------------------------------------ bias field
def _linear_weights(n_in: int, positions: np.ndarray) -> np.ndarray:
    """Matrix that interpolates n_in samples linearly at *positions* (in samples)."""
    positions = np.clip(positions, 0, n_in - 1)
    low = np.minimum(np.floor(positions).astype(int), max(n_in - 2, 0))
    frac = positions - low
    weights = np.zeros((positions.size, n_in), dtype=np.float64)
    rows = np.arange(positions.size)
    weights[rows, low] += 1 - frac
    weights[rows, np.minimum(low + 1, n_in - 1)] += frac
    return weights


def _along_axes(volume: np.ndarray, matrices) -> np.ndarray:
    """Apply one matrix per axis: out = M0 x M1 x M2 applied to the volume."""
    out = volume
    for axis, matrix in enumerate(matrices):
        out = np.moveaxis(np.tensordot(matrix, out, axes=(1, axis)), 0, axis)
    return out


def _resize_linear(volume: np.ndarray, shape) -> np.ndarray:
    """Trilinear interpolation onto a grid whose corner voxels are those of the volume."""
    return _along_axes(volume, [_linear_weights(n, np.linspace(0, n - 1, m))
                                for n, m in zip(volume.shape, shape)])


def simulated_bias_field(shape, vx, percent: float = 30.0, strength: int = 2,
                         seed: int = 0) -> np.ndarray:
    """A smooth random bias field with a mean of 1.

    A random field of ``2**strength`` voxels per axis, low-pass filtered,
    brought to the size of the image and smoothed with 15 mm: a larger
    *strength* gives a field with more local detail.  For one seed it is the
    field ``mri_simulate.m`` makes (``rf.type = [strength, seed]``).

    Args:
        shape: Size of the image.
        vx: Its voxel size in mm.
        percent: Peak-to-peak amplitude in percent, over the whole volume;
            inside the brain it is a little less (27 of 30 on the scan this
            was measured on).  A negative value swaps the bright and the
            dark areas.
        strength: 1 to 4.
        seed: Seed of the field.

    Returns:
        ``shape`` float32, the field the image is multiplied with.
    """
    shape = tuple(int(n) for n in shape)
    vx = voxel_size(vx)
    if round(strength) < 1:
        raise ValueError("the strength of a simulated bias field is at least 1")
    n = 2 ** int(round(strength))
    # MATLAB's rand(n, n, n) for rng(seed, 'twister'): the same generator,
    # filled column by column.  MATLAB starts its seed 0 from 5489, the
    # generator's own default.
    state = np.random.RandomState(int(seed) if int(seed) != 0 else 5489)
    field = state.random_sample(n ** 3).reshape((n, n, n), order="F")

    axis = np.arange(-n // 2, n // 2)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    lowpass = np.exp(-np.sqrt(x ** 2 + y ** 2 + z ** 2) / (n / 8))
    field = np.real(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(field)) * lowpass)))

    # The field is smooth, so it is built and smoothed on a coarse grid of at
    # least six samples per FWHM and only then brought to the image
    coarse_vx = np.maximum(vx, BIAS_FWHM_MM / 6)
    coarse = np.maximum(4, np.ceil(np.array(shape) * vx / coarse_vx)).astype(int)
    fwhm = BIAS_FWHM_MM / (np.array(shape) * vx / coarse)        # in coarse voxels
    pad = int(np.ceil(3 * fwhm.max()))

    field = _resize_linear(field, coarse)
    field -= field.min()
    field /= field.max()
    padded = np.full(coarse + 2 * pad, 0.5)      # 0.5 around it: no dark rim from the smoothing
    padded[pad:-pad, pad:-pad, pad:-pad] = field
    field = smooth(padded, fwhm, 1.0)[pad:-pad, pad:-pad, pad:-pad].astype(np.float64)

    field = _resize_linear(field, shape)
    field -= field.min()
    field /= field.max()
    field *= percent / 100.0                     # the sign of percent inverts the field
    return (1 + field - field.mean()).astype(np.float32)


def template_bias_field(name, iy: np.ndarray, template_mat: np.ndarray,
                        percent: float = 30.0, folder=None) -> np.ndarray:
    """One of the predefined bias fields ``A``, ``B`` or ``C``, on the working grid.

    The fields are measured ones in template space (``rf100_<name>.nii`` of
    the MATLAB tools); they follow the deformation of the segmentation.

    Args:
        name: ``"A"``, ``"B"`` or ``"C"``, or the file of a field in template
            space.
        iy, template_mat: The deformation (``frontend.segment(...,
            deformation=True)``).
        percent: Percent of the stored field.  With 30 the three fields
            spanned 30 to 31% over the brain of the scan this was measured
            on (26 to 28% over its central 98%), like a simulated field of
            30%.  Negative inverts the field.
        folder: Where ``rf100_<name>.nii`` lies; default: ``matlab/`` of the
            repository.

    Returns:
        (X, Y, Z) float32 with a mean of 1.  Where the deformation is not
        defined or the stored field ends, it is continued with its nearest
        value.
    """
    import nibabel as nib

    from . import atlas
    from ._volume import nearest_index

    path = Path(str(name))
    if not path.is_file():
        folder = Path(folder) if folder is not None else Path(__file__).resolve().parents[2] / "matlab"
        path = folder / f"rf100_{name}.nii"
        if not path.is_file():
            raise FileNotFoundError(
                f"The bias field '{name}' was not found: {path} does not exist.  The predefined "
                "fields are part of the repository (matlab/rf100_A.nii, _B, _C).")
    stored = nib.load(path)
    stored = nib.Nifti1Image(_repair_edges(np.asarray(stored.dataobj, dtype=np.float32)), stored.affine)
    field = atlas.to_grid(stored, iy, template_mat, order=1, fill=np.nan)
    defined = np.isfinite(field)
    if not defined.any():
        raise ValueError(f"the bias field {path} does not cover the grid")
    mean = float(field[defined].mean(dtype=np.float64))
    if not defined.all():
        field = field[nearest_index(defined)]        # a step to 1 would be an edge in the image
    return (1 + (field - mean) * (percent / 100.0)).astype(np.float32)


def _repair_edges(field: np.ndarray) -> np.ndarray:
    """Replace an outermost plane that holds about half of the plane next to it.

    Such a plane is what interpolating a field against the zeros beyond its
    border leaves behind.  The first axial plane of ``rf100_A``, ``_B`` and
    ``_C`` is one: at 30% it darkened the lowest 2 mm of the grid to a fifth.
    """
    field = field.copy()
    for axis in range(3):
        for edge, inner in ((0, 1), (-1, -2)):
            outer, next_to = np.take(field, edge, axis), np.take(field, inner, axis)
            valid = next_to != 0
            if valid.any() and 0.4 < float(np.median(outer[valid] / next_to[valid])) < 0.6:
                index = [slice(None)] * 3
                index[axis] = edge
                field[tuple(index)] = next_to
    return field


# -------------------------------------------------------------------- contrast
def change_contrast(intensity, exponent: float, white: float):
    """Power-law change of the contrast around the intensity of white matter.

    ``white * (intensity / white) ** exponent``: white matter keeps its
    intensity, and the ratio of any two intensities is raised to *exponent*.
    0.5 lowers the contrast between the tissues (a ratio of GM to WM of 0.76
    becomes 0.87), 1.5 raises it (0.66).

    This is a map for the intensities of pure tissues, not for a finished
    image.  A voxel that holds two tissues is the mixture of their
    intensities, also at another contrast, and a power law of the mixture is
    not the mixture of the power laws: applied to the image, an exponent of
    0.5 or 1.5 moved the apparent pial boundary by 0.07 mm.  Applied after a
    bias field, it also halved the field or made it 1.6 times as strong.
    :func:`t1sim.simulate.acquire` therefore maps the tissue intensities and
    what is not brain, and mixes the image again
    (:func:`t1sim.synthesis.render`).

    Args:
        intensity: Intensities on the scale of *white*: the three tissue
            intensities, or an image of what is not brain.  Values below 0
            count as 0.
        exponent: Greater than 0; 1 changes nothing.
        white: The intensity of white matter.

    Returns:
        The mapped intensities: an array of float32, or a float for a number.
    """
    if not exponent > 0 or not white > 0:
        raise ValueError("the exponent of a contrast change and the intensity of white matter are positive")
    if np.ndim(intensity) == 0:
        return float(white * (max(float(intensity), 0.0) / white) ** exponent)
    scaled = np.maximum(np.asarray(intensity, dtype=np.float32), 0) / np.float32(white)
    return (scaled ** np.float32(exponent) * np.float32(white)).astype(np.float32)


# ----------------------------------------------------------------------- noise
def add_noise(image: np.ndarray, sigma: float, rician: bool = True, seed=0,
              ceiling: float = None) -> np.ndarray:
    """Noise of a given standard deviation.

    Rician noise is what a magnitude image has: Gaussian noise of *sigma* in
    the real and in the imaginary part of the signal, ``sqrt((s + n1)**2 +
    n2**2)``.  At a high SNR it is close to Gaussian; in the dark background
    it has the Rayleigh floor of a real scan.

    Args:
        image: The image.
        sigma: Standard deviation of the noise, on the scale of the image.
            For a target SNR in white matter: its intensity divided by the
            SNR.
        rician: Rician noise; False: Gaussian noise added to the image.
        seed: Seed, or a ``numpy.random.Generator``.
        ceiling: Upper limit of the result; None for none.  A limit at the
            brightest voxel of the image, which ``mri_simulate.m`` applies,
            cuts off the upper half of the noise wherever the image is that
            bright; measured in its output without a bias field, white
            matter had 61% of the noise it was given and was 1.3% darker.

    Returns:
        The noisy image, float32, not below 0.
    """
    image = np.asarray(image, dtype=np.float32)
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    if sigma <= 0:
        return np.clip(image, 0, ceiling)
    out = image + rng.standard_normal(image.shape, dtype=np.float32) * np.float32(sigma)
    if rician:
        np.square(out, out=out)
        imaginary = rng.standard_normal(image.shape, dtype=np.float32) * np.float32(sigma)
        out += np.square(imaginary, out=imaginary)
        np.sqrt(out, out=out)
    return np.clip(out, 0, ceiling, out=out)


# ------------------------------------------------------------------ resolution
#: Voxels by which an image is continued beyond every face before its k-space is cut.
MIRROR = 16


def _low_pass(image: np.ndarray, keep, axes=(0, 1, 2)) -> np.ndarray:
    """Remove the waves of an image that an acquisition does not hold.

    The Fourier transform takes the image as periodic: what touches one end
    of the grid, as the neck does, would ring into the other end.  The image
    is therefore continued with its mirror image for :data:`MIRROR` voxels
    beyond the faces of *axes* first, which keeps the two ends apart and
    leaves no step at a face.

    Args:
        image: (X, Y, Z) float32.
        keep: Function of the three frequency axes (cycles per voxel, each
            with a length of 1 along the other two axes) that returns the
            waves that stay, as a boolean array.
        axes: The axes along which waves are removed.

    Returns:
        The image without the other waves, float32; the image itself if all
        stay.
    """
    pad = []
    for axis, n in enumerate(image.shape):
        if axis in axes and n > 2:
            total = sp_fft.next_fast_len(n + 2 * min(MIRROR, n - 1), real=True) - n
            total = min(total, 2 * (n - 1))
            pad.append((total // 2, total - total // 2))
        else:
            pad.append((0, 0))
    size = [n + before + after for n, (before, after) in zip(image.shape, pad)]
    freq = [(sp_fft.rfftfreq(n) if axis == 2 else sp_fft.fftfreq(n)).astype(np.float32)
            .reshape([-1 if a == axis else 1 for a in range(3)]) for axis, n in enumerate(size)]
    mask = keep(freq)
    if mask.all():
        return image
    spectrum = sp_fft.rfftn(np.pad(image, pad, mode="reflect"), workers=-1)
    spectrum *= mask
    inner = tuple(slice(before, before + n) for n, (before, _) in zip(image.shape, pad))
    return sp_fft.irfftn(spectrum, s=size, workers=-1)[inner].astype(np.float32)


class Acquisition:
    """The grid of an acquisition at another voxel size, inside the same field of view.

    The new grid has as many whole voxels as fit, centred in the field of
    view, so its voxel size is exactly the one asked for.

    Args:
        shape: Size of the fine grid.
        affine: Its voxel-to-world matrix (axes along the world axes or not).
        voxel: The voxel size of the acquisition in mm, a scalar or one value
            per axis of the fine grid.

    Attributes:
        shape, affine: The grid of the acquisition.
        voxel, fine_voxel: Voxel sizes in mm.
    """

    def __init__(self, shape, affine, voxel):
        self.fine_shape = tuple(int(n) for n in shape)
        self.fine_affine = np.asarray(affine, dtype=np.float64)
        self.fine_voxel = np.sqrt((self.fine_affine[:3, :3] ** 2).sum(0))
        self.voxel = voxel_size(voxel)
        self.ratio = self.voxel / self.fine_voxel
        extent = np.array(self.fine_shape) * self.fine_voxel
        self.shape = tuple(int(n) for n in np.maximum(np.floor(extent / self.voxel + 1e-6), 1))
        #: position of the centre of the first new voxel, in fine voxels
        self.first = ((extent - np.array(self.shape) * self.voxel) / 2 + self.voxel / 2) / self.fine_voxel - 0.5
        to_fine = np.eye(4)
        to_fine[:3, :3] = np.diag(self.ratio)
        to_fine[:3, 3] = self.first
        self.affine = self.fine_affine @ to_fine

    def _overlap(self, axis: int) -> np.ndarray:
        """Share of every fine voxel in every new voxel along one axis (rows sum to 1)."""
        n, m, r = self.fine_shape[axis], self.shape[axis], self.ratio[axis]
        low = self.first[axis] - r / 2 + r * np.arange(m)            # edges, in fine voxels
        edges = np.arange(n + 1) - 0.5
        overlap = (np.minimum(low[:, None] + r, edges[None, 1:])
                   - np.maximum(low[:, None], edges[None, :-1]))
        return np.clip(overlap, 0, None) / r

    def acquire(self, image: np.ndarray) -> np.ndarray:
        """The image at the voxel size of the acquisition.

        k-space is cut at the Nyquist frequency of the new voxel size along
        every axis that gets coarser (a rectangular window: the point spread
        function of a Cartesian acquisition, Gibbs ringing included), and the
        band-limited image is read at the new voxel centres.
        """
        image = np.asarray(image, dtype=np.float32)
        if image.shape != self.fine_shape:
            raise ValueError(f"the image must be {self.fine_shape}, not {image.shape}")
        coarser = [int(axis) for axis in np.where(self.ratio > 1 + 1e-6)[0]]
        if coarser:
            def keep(freq):                             # in cycles per fine voxel
                mask = np.ones((), dtype=bool)
                for axis in coarser:
                    mask = mask & (np.abs(freq[axis]) <= 0.5 / self.ratio[axis] + 1e-6)
                return mask

            image = _low_pass(image, keep, axes=coarser)
        if np.allclose(self.ratio, 1) and np.allclose(self.first, 0):
            return image
        return ndimage.affine_transform(image, self.ratio, offset=self.first,
                                        output_shape=self.shape, order=3, mode="nearest")

    def average(self, volume: np.ndarray) -> np.ndarray:
        """A map of fractions (or a label) at the voxel size of the acquisition.

        Every new voxel is the mean over the box it covers, which is the
        partial volume of that voxel.  A last axis of classes is kept.
        """
        volume = np.asarray(volume, dtype=np.float32)
        if volume.shape[:3] != self.fine_shape:
            raise ValueError(f"the volume must be {self.fine_shape}, not {volume.shape[:3]}")
        return _along_axes(volume, [self._overlap(a) for a in range(3)]).astype(np.float32)

    def to_fine(self, image: np.ndarray) -> np.ndarray:
        """An acquired image back on the fine grid, by cubic B-spline interpolation.

        This is what happens to a real scan that is resampled onto a finer
        grid: its resolution stays that of the acquisition.
        """
        image = np.asarray(image, dtype=np.float32)
        if image.shape != self.shape:
            raise ValueError(f"the image must be {self.shape}, not {image.shape}")
        return ndimage.affine_transform(image, 1 / self.ratio, offset=-self.first / self.ratio,
                                        output_shape=self.fine_shape, order=3, mode="nearest")


def _box_counts(linear: np.ndarray) -> np.ndarray:
    """Into how many slices a new voxel is cut along each of its axes, for :meth:`ObliqueAcquisition.average`.

    One box on a grid along the axes of the fine one, where it is exact;
    otherwise 2 or 3 slices along every axis, the more the further the
    grids are turned against each other, and along a longer axis of the
    voxel as many more as it is longer, so that the slices stay near cubes.

    Args:
        linear: (3, 3) the new grid in the fine one: column j is axis j of
            a new voxel, in fine voxels.
    """
    rows = np.sqrt((linear ** 2).sum(1))
    # how far a fine axis is from the new axis next to it: 0 on a grid along the axes
    turned = float(1 - (np.abs(linear).max(1) / rows).min())
    if turned < 1e-9:
        return np.ones(3, dtype=int)
    per_axis = 2 if turned < 0.06 else 3
    lengths = np.sqrt((linear ** 2).sum(0))
    return per_axis * np.maximum(1, np.round(lengths / lengths.min())).astype(int)


#: How far, in fine voxels, a new voxel centre may lie beyond a face of the
#: fine grid and still count as inside: the rounding of two affines.
FACE_TOLERANCE = 1e-3


class ObliqueAcquisition:
    """The grid of an acquisition that lies anywhere in the world of the image.

    Other voxel sizes, another orientation and another field of view: the
    grid of the scan the working grid was made from.  It does what
    :class:`Acquisition` does for a grid along the axes of the image: the
    image is acquired (its k-space cut at the Nyquist frequencies of the new
    lattice, then read at the new voxel centres), and a map of fractions is
    averaged over the box of every new voxel.  The new voxels are
    parallelepipeds in the image, so neither is separable here.

    Args:
        shape: Size of the fine grid the image lies on.
        affine: Its voxel-to-world matrix.
        target_shape: Size of the grid of the acquisition.
        target_affine: Its voxel-to-world matrix, into the same world.

    Attributes:
        shape, affine: The grid of the acquisition.
        voxel: Its voxel size, in mm of that world.
        linear, offset: Where a voxel of the acquisition lies in the fine
            grid: ``linear @ index + offset``, in fine voxels.
    """

    def __init__(self, shape, affine, target_shape, target_affine):
        self.fine_shape = tuple(int(n) for n in shape)
        self.fine_affine = np.asarray(affine, dtype=np.float64)
        self.shape = tuple(int(n) for n in target_shape)
        self.affine = np.asarray(target_affine, dtype=np.float64)
        if len(self.shape) != 3 or len(self.fine_shape) != 3:
            raise ValueError("both grids have three dimensions")
        self.voxel = np.sqrt((self.affine[:3, :3] ** 2).sum(0))
        to_fine = np.linalg.solve(self.fine_affine, self.affine)
        self.linear, self.offset = to_fine[:3, :3].copy(), to_fine[:3, 3].copy()

    def _check(self, volume: np.ndarray) -> np.ndarray:
        volume = np.asarray(volume, dtype=np.float32)
        if volume.shape[:3] != self.fine_shape:
            raise ValueError(f"the volume must be {self.fine_shape}, not {volume.shape[:3]}")
        return volume

    def acquire(self, image: np.ndarray) -> np.ndarray:
        """The image on the grid of the acquisition.

        A wave of the image with ``f`` cycles per fine voxel has ``linear.T @
        f`` cycles per new voxel.  What exceeds half a cycle along one of the
        three axes of the new lattice lies outside its k-space and is
        removed; the band-limited image is then read at the new voxel centres
        by cubic B-spline interpolation.  A new voxel whose centre lies
        beyond a face of the fine grid is 0; within :meth:`coverage` of such a
        face the image is not to be trusted.
        """
        image = self._check(image)
        if image.ndim != 3:
            raise ValueError("the image is one volume")

        def keep(freq):
            mask = np.ones((), dtype=bool)
            for new_axis in range(3):
                cycles = sum(np.float32(self.linear[axis, new_axis]) * freq[axis] for axis in range(3))
                mask = mask & (np.abs(cycles) <= 0.5 + 1e-6)
            return mask

        image = ndimage.affine_transform(_low_pass(image, keep), self.linear, offset=self.offset,
                                         output_shape=self.shape, order=3, mode="nearest")
        image[self.coverage(ramp=1.0) <= 0] = 0.0
        return image

    def sample(self, volume: np.ndarray, extend: bool = False) -> np.ndarray:
        """A smooth map at the centres of the new voxels, by linear interpolation.

        For a map that changes little within a voxel, a bias field, where an
        average over the voxel would be the same.

        Args:
            volume: (X, Y, Z) map on the fine grid.
            extend: Continue the map beyond the fine grid with its nearest
                value; False: 0 there.
        """
        return ndimage.affine_transform(self._check(volume), self.linear, offset=self.offset,
                                        output_shape=self.shape, order=1,
                                        mode="nearest" if extend else "constant", cval=0.0)

    def average(self, volume: np.ndarray) -> np.ndarray:
        """A map of fractions on the grid of the acquisition.

        Every new voxel is the mean over the parallelepiped it covers in the
        fine grid, which is the partial volume of that voxel.  A fine voxel
        counts with the share of it that lies inside, as in
        :meth:`Acquisition.average`: the map is taken as constant within a
        fine voxel, which is what a fraction of that voxel is.

        The parallelepiped is cut into slices along its own axes, and every
        slice is covered with a box along the axes of the fine grid that has
        its spread along each of them; the mean over a box is exact: the
        integral of a map that is constant within its voxels is the linear
        interpolation of its running sums.  On a grid along the axes the one
        box is the new voxel and the result is that of
        :meth:`Acquisition.average`.  Otherwise there are 2 or 3 slices along
        every axis, the more the further the grids are turned against each
        other, and along a longer axis of the voxel as many more as it is
        longer (:func:`_box_counts`).  The boxes fit the slices only roughly,
        which was measured on a shell whose boundaries are one fine voxel
        wide, against the same boxes with 16 to 20 slices along every axis,
        which have converged to 0.0005 (6 to 8 slices for the long voxels,
        to 0.0015):

        =======================  =========  =====  =============  =======
        New voxels               Turned by  Boxes  Largest error  RMS at
                                                                  boundary
        =======================  =========  =====  =============  =======
        1.7 fine voxels          1 degree   8      0.0004         0.0001
        1.7 fine voxels          6          8      0.003          0.0006
        1.7 fine voxels          11         8      0.005          0.0011
        1.7 fine voxels          17         8      0.009          0.0021
        1.7 fine voxels          25 to 29   27     0.010          0.0024
        2 x 2 x 6 fine voxels    25         81     0.002          0.0007
        1 x 1 x 4 fine voxels    25 and 15  108    0.001          0.0003
        =======================  =========  =====  =============  =======

        With the same number of slices along every axis, the last two were
        0.019 and 0.037 off.  Counting the fine voxel at points inside the
        new voxel, the obvious way, needs 40 points per axis to come as close
        for the grid turned by 1 degree: the points fall in step with the
        fine grid.  A last axis of classes is kept; beyond the fine grid the
        map is 0.
        """
        volume = self._check(volume)
        classes = volume.ndim == 4
        stack = np.moveaxis(volume, -1, 0) if classes else volume[None]
        out = np.zeros((stack.shape[0],) + self.shape, dtype=np.float32)
        filled = np.argwhere(stack.any(0))
        if not filled.size:
            return np.ascontiguousarray(np.moveaxis(out, 0, -1)) if classes else out[0]

        # the part of the fine grid that is not 0, and the new voxels that reach it
        first, last = filled.min(0), filled.max(0) + 1
        corners = np.array(list(itertools.product(*zip(first - 1.0, last + 0.0))))
        index = np.linalg.solve(self.linear, (corners - self.offset).T).T
        low = np.clip(np.floor(index.min(0)).astype(int) - 1, 0, self.shape)
        high = np.clip(np.ceil(index.max(0)).astype(int) + 2, 0, self.shape)
        size = tuple(int(b - a) for a, b in zip(low, high))
        if min(size) < 1:
            return np.ascontiguousarray(np.moveaxis(out, 0, -1)) if classes else out[0]
        box = (slice(None),) + tuple(slice(a, b) for a, b in zip(low, high))

        # the new voxel cut into counts[j] slices along its axis j, each covered by a box along the
        # fine axes with the spread of the slice along every fine axis (half-widths in fine voxels)
        counts = _box_counts(self.linear)
        half = 0.5 * np.sqrt(((self.linear / counts) ** 2).sum(1))
        steps = [(np.arange(n) + 0.5) / n - 0.5 for n in counts]
        calls = [(self.offset + self.linear @ (low + np.array(step)) - first + 0.5 + np.array(signs) * half,
                  float(np.prod(signs)))
                 for step in itertools.product(*steps) for signs in itertools.product((-1, 1), repeat=3)]

        for channel, values in enumerate(stack):
            # running sums along the three axes, at the corners of the fine voxels
            table = np.zeros(tuple(int(n) + 1 for n in last - first), dtype=np.float64)
            table[1:, 1:, 1:] = values[tuple(slice(a, b) for a, b in zip(first, last))]
            for axis in range(3):
                np.cumsum(table, axis=axis, out=table)

            def part(some, table=table):
                total = np.zeros(size, dtype=np.float64)
                for corner, sign in some:
                    total += sign * ndimage.affine_transform(table, self.linear, offset=corner, output_shape=size,
                                                             order=1, mode="nearest", output=np.float64)
                return total

            workers = max(1, min(len(calls), os.cpu_count() or 1, 16))
            total = np.zeros(size, dtype=np.float64)
            with ThreadPoolExecutor(workers) as pool:
                for result in pool.map(part, [calls[n::workers] for n in range(workers)]):
                    total += result
            mean = total / (np.prod(counts) * np.prod(2 * half))
            mean[np.abs(mean) < 1e-7] = 0.0             # what the rounding of the sums leaves of nothing
            out[(channel,) + box[1:]] = mean
        return np.ascontiguousarray(np.moveaxis(out, 0, -1)) if classes else out[0]

    def coverage(self, ramp: float = 6.0) -> np.ndarray:
        """How far every new voxel lies inside the fine grid.

        It falls only towards the faces of the fine grid that the centre of
        a new voxel lies beyond, by more than a thousandth of a fine voxel:
        where the new grid ends at a face of the fine one or before it, as at
        the faces of the field of view of a scan that the grid of the
        anatomy reaches, there is nothing to make room for.  A face is the
        outer face of the outermost fine voxels, half a voxel beyond their
        centres.

        Args:
            ramp: Width in fine voxels of the border in which the result
                falls from 1 to 0.

        Returns:
            (X, Y, Z) float32 on the grid of the acquisition: 1 inside the
            fine grid, 0 at a face with new voxels beyond it, and beyond.
        """
        index = [np.arange(n, dtype=np.float32).reshape([-1 if a == axis else 1 for a in range(3)])
                 for axis, n in enumerate(self.shape)]
        weight = np.ones(self.shape, dtype=np.float32)
        for axis, size in enumerate(self.fine_shape):
            position = sum(np.float32(self.linear[axis, new]) * index[new] for new in range(3))
            position = position + np.float32(self.offset[axis])
            # the faces of the fine grid; a face is only a border where new voxels lie beyond it
            lower, upper = np.float32(-0.5), np.float32(size - 0.5)
            below = position - lower if float(position.min()) < lower - FACE_TOLERANCE else np.float32(np.inf)
            above = upper - position if float(position.max()) > upper + FACE_TOLERANCE else np.float32(np.inf)
            inside = np.minimum(below, above) / np.float32(ramp)
            np.minimum(weight, np.clip(inside, 0, 1), out=weight)
        return weight
