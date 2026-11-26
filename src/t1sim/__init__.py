"""t1sim -- simulated T1-weighted MR images with a known ground truth.

The Python counterpart of ``matlab/mri_simulate.m``, for validation phantoms
and for training data.  A real scan is turned into tissue fractions
(:mod:`t1sim.frontend`), which can be cleaned (:mod:`t1sim.cleanup`) and
changed (:mod:`t1sim.anatomy`); the image is rendered from them
(:mod:`t1sim.synthesis`) and acquired with a bias field, a voxel size, a
contrast, noise (:mod:`t1sim.artefacts`), movement and ringing
(:mod:`t1sim.kspace`).  :mod:`t1sim.simulate` runs these steps and writes
the results after BIDS (:mod:`t1sim.bids`).  :mod:`t1sim.phantom` builds a
folded phantom whose cortical thickness is exact by construction.
"""

__version__ = "0.1.0.dev0"
