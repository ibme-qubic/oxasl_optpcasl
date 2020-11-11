"""
OXASL_OPTPCASL - Basic data structures used

Copyright 2019 University of Oxford
"""
import numpy as np

class ScanParams(object):
    """
    Parameters of the scan to optimize for
    """
    def __init__(self, duration, npld, nslices=1, slicedt=0.0, readout=0.638, ld=1.8, noise=0.0013, plds=None, had_size=8):
        self.duration = duration
        self.npld = npld
        self.nslices = nslices
        self.slicedt = slicedt
        self.readout = readout
        if isinstance(ld, (float, int)):
            self.ld = [float(ld)]
        else:
            self.ld = ld
        self.plds = plds
        self.noise = noise
        self.had_size=had_size
        if self.plds is not None and len(self.plds) != self.npld:
            raise ValueError("Supplied initial PLDs %s is inconsistent with stated number of PLDs: %i" % (plds, npld))

    def __str__(self):
        if self.nslices > 1:
            return "%is 2D scan with %.3fs label duration, %i slices (time per slice=%.5fs) and readout time %.3fs" % (self.duration, self.ld, self.nslices, self.slicedt, self.readout)
        else:
            return "%is 3D scan with  %.3fs label duration and readout time %.3fs" % (self.duration, self.ld, self.readout)

class PhysParams(object):
    """
    Define the physiological parameters to use

    (as per Buxton et al. MRM 1998)
    """
    def __init__(self, **kwargs):
        self.t1b = kwargs.get("t1b", 1.65)
        self.t1t = kwargs.get("t1t", 1.445)
        self.alpha = kwargs.get("alpha", 0.85)
        self.lam = kwargs.get("lam", 0.9)
        self.f = kwargs.get("f", 50.0/6000)
        self.m0b = kwargs.get("m0b", 1.0)

class ATTDist(object):
    """
    ATT (BAT) distribution

    The 'taper' parameter (seconds) causes the weighting to decay from 1 to 0.5
    at the beginning and end of the distribution
    """
    def __init__(self, start=0.5, end=2.0, step=0.025, taper=0):
        total_points = int(1 + np.ceil((end-start)/step))
        taper_points = int(np.ceil(taper / step)) - 1 # FIXME I don't agree with this but it is compatible with Joe's code
        self.start, self.end, self.step, self.taper = start, end, step, taper
        self.atts = np.linspace(start, end, total_points)
        self.exclude_taper = np.linspace(start+taper, end-taper, total_points - 2*taper_points)
        self.weight = np.concatenate((
            np.linspace(0.5, 1.0, taper_points+1),
            np.ones(total_points - 2*(taper_points+1)),
            np.linspace(1.0, 0.5, taper_points+1),
        ))
        self.length = len(self.atts)

    def __str__(self):
        return "ATT distribution: %i values between %.3fs and %.3fs (weight taper=%.3fs)" % (self.length, self.start, self.end, self.taper)

class Limits(object):
    """
    Limits and step size to search over for PLDs / labelling durations
    """
    def __init__(self, lb, ub, step, name="PLD"):
        self.lb, self.ub, self.step, self.name = lb, ub, step, name

    def __str__(self):
        return "%ss between %.3fs and %.3fs in steps of %.5fs" % (self.name, self.lb, self.ub, self.step)
