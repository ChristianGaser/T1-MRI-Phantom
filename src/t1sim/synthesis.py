"""Synthesis: the T1-weighted image that belongs to tissue fractions.

Inside the brain the image is the mixture of the tissue intensities, weighted
by the fractions of CSF, GM and WM: the linear partial volume model, and the
reason the fractions are the exact ground truth of the image.  Everything
that is not brain keeps the intensity of the bias-corrected scan, so skull,
scalp and background stay as they are, and where the fractions sum to less
than one the two blend.

The tissue intensities are those of the mixture of Gaussians the
segmentation estimated (:func:`t1sim.frontend.tissue_means`); nothing is
fitted here.
"""

from __future__ import annotations

import numpy as np

from ._volume import percentile
from .frontend import TISSUES


def render(fractions: np.ndarray, means, corrected: np.ndarray,
           clip_percentile: float = 99.99) -> np.ndarray:
    """The image of tissue fractions.

    ``(sum_k mean_k * fraction_k + corrected * (s - P)) / s`` with ``P`` the
    sum of the fractions and ``s = max(P, 1)``, as ``mri_simulate.m``
    synthesises it.

    Args:
        fractions: (X, Y, Z, 3) fractions of CSF, GM and WM.
        means: Their intensities: a dict with the keys ``CSF``, ``GM`` and
            ``WM``, or three values in that order.
        corrected: (X, Y, Z) bias-corrected scan, which fills what is not
            brain.  Voxels that are 0 or not finite there carry no data and
            get a value just above 0.  Negative values, which the resampling
            of the scan leaves in the background, become 0: a magnitude image
            has none.
        clip_percentile: The image is limited to this percentile of all its
            voxels, because the bias correction can leave single very bright
            voxels in the noisy background.  None: no limit.
            ``mri_simulate.m`` limits at the 98th percentile, which on a
            whole head at 0.5 mm is the intensity of white matter: more than
            2% of the voxels are pure white matter, so the scalp fat, 1.4 to
            3 times as bright in the scans this was measured on, was cut down
            to it.

    Returns:
        (X, Y, Z) float32 image on the intensity scale of *corrected*.
    """
    fractions = np.asarray(fractions)
    if fractions.ndim != 4 or fractions.shape[-1] != 3:
        raise ValueError(f"fractions must be (X, Y, Z, 3), not {fractions.shape}")
    if isinstance(means, dict):
        means = [means[t] for t in TISSUES]
    means = np.asarray(means, dtype=np.float32)
    if means.shape != (3,):
        raise ValueError("means must hold the intensity of CSF, GM and WM")
    corrected = np.asarray(corrected, dtype=np.float32)
    if corrected.shape != fractions.shape[:3]:
        raise ValueError("corrected must be on the grid of the fractions")

    brain = fractions.sum(-1, dtype=np.float32)
    scale = np.maximum(brain, 1.0)
    image = fractions @ means
    no_data = (corrected == 0) | ~np.isfinite(corrected)
    image += np.where(no_data, 0.0, np.maximum(corrected, 0)) * (scale - brain)
    image /= scale
    image[no_data] = 1e-3
    if clip_percentile is not None:
        np.minimum(image, percentile(image, clip_percentile), out=image)
    return image
