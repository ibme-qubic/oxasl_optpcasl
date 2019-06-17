"""
OCASL_OPTPCASL

Python library for optimizing multi-PLD pCASL acquisitions

@author Joseph Woods, FMRIB Oxford July 2017

Ported to Python from the MATLAB code by Martin Craig

Copyright 2019 University of Oxford
"""
from .optimize import optimize, LOptimal, DOptimal, TRWeightingOrNAveFloor
from .structures import ASLParams, BATDist, Scan, Limits
