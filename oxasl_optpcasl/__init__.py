"""
OCASL_OPTPCASL

Python library for optimizing multi-PLD pCASL acquisitions

@author Joseph Woods, FMRIB Oxford July 2017

Ported to Python from the MATLAB code by Martin Craig

Copyright 2019 University of Oxford
"""
import copy

import numpy as np

from OED_CASL_OptPLD_acrossSlices_analytical_seq_minDist import optimize, TRWeightingOrNAveFloor

class ASLParams:
    """
    Define the ASL parameters (as per Buxton et al. MRM 1998)
    """
    def __init__(self, filename, f, **kwargs):
        self.filename, self.f = filename, f
        self.bat = kwargs.get("bat", 1.0)
        self.m0b = kwargs.get("m0b", 1.0)
        self.t1b = kwargs.get("t1b", 1.65)
        self.t1t = kwargs.get("t1t", 1.445)
        self.lam = kwargs.get("lam", 0.9)
        self.alpha = kwargs.get("alpha", 0.85)
        self.tau = kwargs.get("tau", 1.4)
        self.noise = kwargs.get("noise", 0.002)

class BATDist:
    """
    ATT (BAT) distribution

    The 'taper' parameter (seconds) causes the weighting to decay from 1 to 0.5
    at the beginning and end of the distribution
    """
    def __init__(self, start, end, step, taper=0):
        self.start, self.end, self.step, self.taper = start, end, step, taper
        self.dist = np.arange(start, end+step, step)
        self.exclude_taper = np.arange(start+taper, end-taper+step, step)
        self.weight = np.concatenate((
            np.linspace(0.5, 1.0, np.floor(taper / step)),
            np.ones(int(np.ceil((end - start - 2*taper) / step))),
            np.linspace(1.0, 0.5, np.floor(taper / step)),
        ))
        self.length = len(self.dist)

    def __str__(self):
        return "BAT distribution: %i values between %.2fs and %fs (weight taper=%.2fs)" % (self.length, self.start, self.end, self.taper)

class Scan:
    """
    Details of the desired scan to optimize for
    """
    def __init__(self, duration, npld, slices=1, slicedt=0.0, readout=0.5):
        self.duration, self.npld, self.slices, self.slicedt, self.readout = duration, npld, slices, slicedt, readout

    def __str__(self):
      if self.slices > 1:
          return "%is 2D scan with %i slices (time per slice=%.5fs) and readout time %.2fs" % (self.duration, self.slices, self.slicedt, self.readout)
      else:
          return "%is 3D scan with readout time %fs" % (self.duration, self.readout)

class Limits:
    """
    PLD limits and step size to search over
    """
    def __init__(self, lb, ub, step):
        self.lb, self.ub, self.step = lb, ub, step

    def __str__(self):
        return "PLDs between %.2fs and %.2fs in steps of %.5fs" % (self.lb, self.ub, self.step)

class OptimizationMethod:

    def resize_inputs(self, params):
        """
        Returns a modified version of params with t, tau, f and BAT
        all set to a common shape
        """
        params = copy.copy(params)
        params.tau = np.full(params.t.shape, params.tau)
        
        A_f = params.f
        A_BAT = params.bat

        # Tile t and tau so third dimension has size of BAT
        params.t = np.squeeze(np.tile(params.t[..., np.newaxis], [1, 1, params.bat.shape[0]]))
        params.tau = np.squeeze(np.tile(params.tau[..., np.newaxis], [1, 1, params.bat.shape[0]]))

        # Make f and BAT the same dimensions as t
        # f is constant, BAT is 1D (which should be preserved in the 3rd dimension)
        params.f = np.full(params.t.shape, params.f)
        if len(params.t.shape) == 2:
            params.bat = np.tile(params.bat[np.newaxis, :], list(params.t.shape)[:-1] + [1])
        else:
            params.bat = np.tile(params.bat[np.newaxis, np.newaxis, :], list(params.t.shape)[:-1] + [1])
        #print("Updated shapes in HessianParseInputs")
        #print("t shape ", params.t.shape)
        #print("tau shape ", params.tau.shape)
        #print("f shape ", params.f.shape)
        #print("BAT shape ", params.bat.shape)

        # FIXME deleted MATLAB error handling here for 1 PLD - test?
        return params

    def sensitivity(self, params):
        """
        This function calculates the sensitivty functions of the Buxtion CASL
        model (Buxton et al. MRM 1998) and are given in Woods et al. MRM 2019.
        The BAT sensitivty function assumes a fixed outflow in T1prime in order
        to simplify the equations.
        """

        # This sensitivty function includes an assumed fixed T1prime, so we fix the
        # f in T1prime to a sensible value.
        fixedF = 50.0/6000

        T1prime = 1.0/((1.0/params.t1t) + (fixedF/params.lam))
        #print(params.m0b, params.alpha, T1prime, params.t1b)
        #print(np.mean(2*params.m0b * params.alpha * T1prime))
        #print(params.bat[..., 0])
        M = 2*params.m0b * params.alpha * T1prime * np.exp(-params.bat/params.t1b)
        #print("BATo ", params.bat.shape, np.mean(params.bat))
        #print("M ", M.shape, np.mean(M))

        # Initialise
        df = np.zeros(params.t.shape, dtype=np.float32)
        dBAT = np.zeros(params.t.shape, dtype=np.float32)

        # for t between deltaT and tau plus deltaT
        tRegion = np.logical_and(params.t > params.bat, params.t <= (params.tau + params.bat))
        t = params.t[tRegion]
        # if sum(size(params.bat)>1) > 1 # Check whether multiple values are being used
        if np.any(np.array(params.bat.shape) > 1):
            BAT = params.bat[tRegion]
            f = params.f[tRegion]
            M_temp = M[tRegion]
        else:
            BAT = params.bat
            f = params.f
            M_temp = M

        T1prime_temp = T1prime

        df[tRegion] = M_temp * (1 - np.exp((BAT-t) / T1prime_temp))
        dBAT[tRegion] = M_temp * f * ((-1.0/params.t1b) - np.exp((BAT-t) / T1prime_temp) * ((1.0/T1prime_temp) - (1.0/params.t1b)))
        #print("BAT ", BAT.shape, np.mean(BAT))
        #print("M_temp ", M_temp.shape, np.mean(M_temp))
        #print("f ", f.shape, np.mean(f))
        #print("dBAT ", dBAT.shape, np.mean(dBAT))
        #print("df ", df.shape, np.mean(df))

        # for t greater than tau plus deltaT
        tRegion = params.t > (params.tau+params.bat)
        t = params.t[tRegion]
        if np.any(np.array(params.tau.shape) >1):
            tau = params.tau[tRegion]
        else:
            tau = params.tau

        # if sum(size(params.bat)>1) > 1
        if np.any(np.array(params.bat.shape) > 1):
            BAT = params.bat[tRegion]
            f = params.f[tRegion]
            M_temp = M[tRegion]
        else:
            BAT = params.bat
            f = params.f
            M_temp = M
        T1prime_temp = T1prime

        df[tRegion] = M_temp * np.exp((-t+tau+BAT) / T1prime_temp) * (1-np.exp(-tau/T1prime_temp))
        dBAT[tRegion] = M_temp * f * (1-np.exp(-tau/T1prime_temp)) * np.exp((BAT+tau-t)/T1prime_temp) * ((1.0/T1prime_temp)-(1.0/params.t1b))

        return df, dBAT

    def Hessian(self, params, scan, slice, tDim=0, **kwargs):
        """
        This function calculates an approximation of the Hessian for optimising
        an experimental design. param is the struct containing the variables,
        
        Inputs:
                param    - A struct containing all of the parameters needed for the
                        Buxton standard CASL model (Buxton et al. MRM 1998)
                scanTime - The scan duration avalaible for the ASL protocol (units: seconds)
                slice    - The number of slices to optimise for (units: seconds)
                scan.slicedt  - The time to acquire each slice (the inter-slice
                        timing) (units: seconds)
        """
        # Get parameters into the correct matrix sizes
        params = self.resize_inputs(params)

        # Work out how many average we can fit in and divide by the noise SD
        TRWeight = kwargs.get("nAverage", None)
        if TRWeight is None:
            TRWeight, _ = TRWeightingOrNAveFloor(params, scan, tDim, slice)

        TRWeight = TRWeight/(params.noise**2)
        #print("TRWeight shape ", TRWeight.shape, np.mean(TRWeight))

        df, dBAT = self.sensitivity(params)
        #print("df shape ", df.shape, np.mean(df))
        #print("dBAT shape ", dBAT.shape, np.mean(dBAT))

        # Form  the Hessian 
        # Hessian dimensions can be: 2 x 2 x PLD x Tau x BAT x f
        H = np.zeros(list(df.shape)[1:]+ [2, 2])
        H[..., 0, 0] = TRWeight * np.squeeze(np.sum(df*df, tDim))
        H[..., 0, 1] = TRWeight * np.squeeze(np.sum(df*dBAT, tDim))
        H[..., 1, 0] = TRWeight * np.squeeze(np.sum(dBAT*df, tDim))
        H[..., 1, 1] = TRWeight * np.squeeze(np.sum(dBAT*dBAT, tDim))
        return H

    def covOEDNAveFloor(self, params, scan, slice, **kwargs):
        """
        This function calculates an approximation of the Hessian for optimising
        an experimental design. param is the struct containing the variables,
        """
        H = self.Hessian(params, scan, slice, **kwargs)

        # Numpy determinant function batches over leading dimensions
        covMatrix = np.linalg.inv(H)

        # Change into (ml/100g/min)
        covMatrix[..., 0, 0] = covMatrix[..., 0, 0] * 6000 * 6000
        return covMatrix, H

class LOptimal(OptimizationMethod):
    """
    Choose CBF variance to minimise
    """
    def __init__(self, A):
        self.A = A
        self.name = 'L-optimal'

    def HessianCov(self, params, scan, slice):
        H = self.Hessian(params, scan, slice)
        
        # Numpy matrix ops batch over leading dimensions
        cov = np.abs(np.matmul(self.A, np.linalg.inv(H)))

        # Correct for inf*0 errors in A*inverse
        cov[np.isnan(cov)] = np.inf

        # Force trace function to batch across leading dimensions
        r = np.trace(cov, axis1=-1, axis2=-2)
        return r

class DOptimal(OptimizationMethod):
    """
    """
    def __init__(self):
        self.name = 'D-optimal'
        
    def HessianCov(self, params, scan, slice):
        H = self.Hessian(params, scan, slice)
        
        # Numpy determinant function batches over leading dimensions
        detH = np.linalg.det(H)
        #detH = np.square(TRWeight) * detH
        #print("detH shape ", detH.shape)

        return 1.0/np.abs(detH)
    