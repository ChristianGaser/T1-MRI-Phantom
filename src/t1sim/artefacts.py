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
        coarser = self.ratio > 1 + 1e-6
        if coarser.any():
            spectrum = sp_fft.rfftn(image, workers=-1)
            for axis in np.where(coarser)[0]:
                n = self.fine_shape[axis]
                freq = sp_fft.rfftfreq(n) if axis == 2 else sp_fft.fftfreq(n)   # cycles per fine voxel
                keep = np.abs(freq) <= 0.5 / self.ratio[axis] + 1e-9
                spectrum *= keep.reshape([-1 if a == axis else 1 for a in range(3)])
            image = sp_fft.irfftn(spectrum, s=image.shape, workers=-1).astype(np.float32)
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
