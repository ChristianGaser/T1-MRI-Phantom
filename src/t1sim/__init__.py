"""t1sim -- simulated T1-weighted MR images with a known ground truth.

The Python counterpart of ``matlab/mri_simulate.m``, for validation phantoms
and for training data.  In development: so far it holds the download of the
model files (:mod:`t1sim.models`) and the front end that turns a real scan
into a label map and tissue fractions (:mod:`t1sim.frontend`).
"""

__version__ = "0.1.0.dev0"
