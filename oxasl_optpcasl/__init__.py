"""
OCASL_OPTPCASL

Python library for optimizing multi-PLD pCASL acquisitions
Ported to Python from MATLAB
"""
__author__ = "Joseph Woods"
__copyright__ = "Copyright 2019 University of Oxford"
__credits__ = ["Joseph Wood, Martin Craig"]
__maintainer__ = "Martin Craig"
__email__ = "martin.craig@eng.ox.ac.uk"

from ._version import __version__
from .optimize import optimize, LOptimal, DOptimal, TRWeightingOrNAveFloor
from .structures import ASLParams, BATDist, Scan, Limits
