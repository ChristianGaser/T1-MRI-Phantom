"""t1sim -- simulated T1-weighted MR images with a known ground truth.

The Python counterpart of ``matlab/mri_simulate.m``, for validation phantoms
and for training data.  A real scan is turned into tissue fractions
(:mod:`t1sim.frontend`), whose anatomy can be changed (:mod:`t1sim.anatomy`).
The image is rendered from them (:mod:`t1sim.synthesis`) and acquired with a
contrast, a bias field, a voxel size, noise (:mod:`t1sim.artefacts`),
movement and ringing (:mod:`t1sim.kspace`); its ground truth is the same
fractions after a cleanup of what only looks like grey matter
(:mod:`t1sim.cleanup`).  :mod:`t1sim.simulate` runs these steps and writes
the results after BIDS (:mod:`t1sim.bids`); :mod:`t1sim.atlas` brings maps of
the template space onto the grid, and :mod:`t1sim.models` fetches the
templates of the segmentation.  :mod:`t1sim.phantom` builds a folded phantom
whose cortical thickness is known by construction.
"""

__version__ = "0.1.0.dev0"
