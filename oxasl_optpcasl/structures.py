"""
OXASL_OPTPCASL - Basic data structures used

Copyright 2019 University of Oxford
"""
import numpy as np

VAR_MULTI_PCASL = 'var_multi_pCASL'
VAR_TE_PCASL = 'var_te_pCASL'
LOOK_LOCKER = 'look_locker'
VAR_TE_PCASL_NPLD = 'var_te_pCASL_nPLD'

class ASLScan(object):
    """
    Details of the desired scan to optimize for
    """
    def __init__(self, asltype, duration, npld, slices=1, slicedt=0.0, readout=0.5):
        self.asltype = asltype
        self.duration = duration
        self.npld = npld
        self.slices = slices
        self.slicedt = slicedt
        self.readout = readout

    def __str__(self):
        if self.slices > 1:
            return "%is 2D scan with %i slices (time per slice=%.5fs) and readout time %.2fs" % (self.duration, self.slices, self.slicedt, self.readout)
        else:
            return "%is 3D scan with readout time %fs" % (self.duration, self.readout)

class ASLParams(object):
    """
    Define the physiological parameters to use
    
    (as per Buxton et al. MRM 1998)
    """
    def __init__(self, **kwargs):
        self.f = kwargs.get("f", 50.0/6000)
        self.bat = kwargs.get("bat", 1.0)
        self.m0b = kwargs.get("m0b", 1.0)
        self.t1b = kwargs.get("t1b", 1.65)
        self.t1t = kwargs.get("t1t", 1.445)
        self.lam = kwargs.get("lam", 0.9)
        self.alpha = kwargs.get("alpha", 0.85)
        self.tau = kwargs.get("tau", 1.4)
        self.noise = kwargs.get("noise", 0.002)

class ATTDist(object):
    """
    ATT (BAT) distribution

    The 'taper' parameter (seconds) causes the weighting to decay from 1 to 0.5
    at the beginning and end of the distribution
    """
    def __init__(self, start, end, step, taper=0):
        total_points = int(1 + np.ceil((end-start)/step))
        taper_points = int(np.floor(taper / step))
        self.start, self.end, self.step, self.taper = start, end, step, taper
        self.dist = np.linspace(start, end, total_points)
        self.exclude_taper = np.linspace(start+taper, end-taper, total_points - 2*taper_points)
        self.weight = np.concatenate((
            np.linspace(0.5, 1.0, taper_points),
            np.ones(total_points - 2*taper_points),
            np.linspace(1.0, 0.5, taper_points),
        ))
        self.length = len(self.dist)

    def __str__(self):
        return "BAT distribution: %i values between %.2fs and %fs (weight taper=%.2fs)" % (self.length, self.start, self.end, self.taper)

class Limits(object):
    """
    PLD limits and step size to search over
    """
    def __init__(self, lb, ub, step):
        self.lb, self.ub, self.step = lb, ub, step

    def __str__(self):
        return "PLDs between %.2fs and %.2fs in steps of %.5fs" % (self.lb, self.ub, self.step)
