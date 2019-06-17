"""
Functions for optimizing the PLDs of a multi-PLD pCASL protocol

@author Joseph Woods, FMRIB, Oxford, July 2017
"""
import sys
import copy

import numpy as np

from ._version import __version__

def optimize(opttype, params, att_dist, scan, lims, log=sys.stdout):
    """
    Optimise the PLDs of a multi-PLD PCASL protocol
    
    :param opttype: Optimization type object (L-optimal or D-optimal)
    :param params: ASL parameters
    :param att_dist: ATT time distribution and weighting
    :param scan: Scan details to optimize for
    :param lims: PLD limiting values (max/min)
    :param log: Stream-like object for logging output
    """
    params.bat = att_dist.dist
    log.write("OPTPCASL v%s\n\n" % __version__)
    log.write("Optimizing PLDs for %s\n" % scan)
    log.write("PLD search limits: %s\n" % lims)
    log.write("Optimizing for %i PLDs\n" % scan.npld)
    log.write("%s\n" % str(att_dist))
    log.write("Optimization method: %s\n" % opttype.name)

    # Initialise
    bestminVariance = 1e99
    count = 0
    allCost = []
    allPLD = []

    # Initialise the PLDs and check that the number of avearges > 0. Reduce
    # the PLD range if not
    factor = 1.0/lims.step
    PLDSubtract = 0
    TRWeight = 0
    while TRWeight < 1:

        # Initial sequence of PLDs spaced evenly between upper and lower bound
        PLD = np.linspace(lims.lb, lims.ub - PLDSubtract, scan.npld)
        PLD = np.round(PLD*factor) / factor
        params.PLD = PLD
        #print("PLDs", params.PLD)
        # Sequence of TIs corresponding to the PLDs. FIXME tau could be sequence
        params.t = params.tau + params.PLD
        #print("TIs", params.t)
        TRWeight, _ = TRWeightingOrNAveFloor(params, scan, 0, 0)
        #print("TRWeight", TRWeight)

        if TRWeight < 1:
            PLDSubtract += 0.1

    # Optimisation loop
    while 1:
        # Each time through all of the samples, save them
        oldPLD = np.copy(PLD)

        for pld_idx in range(scan.npld):

            #print("\nDoing PLD %i" % pld_idx)

            # Work out the set of available PLDs to search over. This is a series of
            # values between the PLD before the current PLD and the PLD after, with a 
            # step size of lims.step however special case for first and last PLDs
            # where we use the upper/lower bounds instead of the previous/next PLD
            if pld_idx == 0:
                start, stop = lims.lb, PLD[pld_idx+1]
            elif pld_idx == scan.npld-1:
                start, stop = PLD[pld_idx-1], lims.ub
            else:
                start, stop = PLD[pld_idx-1], PLD[pld_idx+1]

            PLDTry = np.round(np.arange(start, stop+0.001, lims.step), 5)
            #print("PLDTry", PLDTry)

            PLDTryL = len(PLDTry)
            if PLDTryL > 0:

                # Current PLDs list is tiled with 1 column per trial PLD
                params.PLD = np.tile(PLD[:, np.newaxis], (1, PLDTryL))

                # Distweight sequence is tiled with 1 row per trial PLD
                distWeightCurr = np.tile(att_dist.weight, (PLDTryL, 1))

                # Flag matrix shape with 1 row per dist weight and 1 column per slice
                distInd = np.zeros((att_dist.length, scan.slices), dtype=np.bool)

                # Variance matrix shape (#trial PLDs, #distweights, #slices)
                variance = np.zeros((PLDTryL, att_dist.length, scan.slices), dtype=np.float)

                for slice in range(scan.slices): # Loop over the slices

                    # Can only calculate the cost function for BAT >
                    # shortest PLD, otherwise we don't see the inflow and
                    # determinant blows up.
                    #
                    # NB Matlab uses min() here on lims.lb suggesting it could be a sequence?
                    minPLDPossible = lims.lb + slice*scan.slicedt

                    # For this slice distInd = True if corresponding BAT is above the minimum PLD
                    distInd[:, slice] = att_dist.dist > minPLDPossible
                    params.bat = att_dist.dist[distInd[:, slice]]

                    # TRWeightingOrNAveFloor only adjusts the readout time.
                    # Adjust the PLDs (and t) here. Add slice*scan.slicedt to
                    # the PLDs. If first slice, add 0s. If 12th slice, add
                    # 0.4972s.

                    # By specifying the slice, this calculation will ensure
                    # that the shortest PLD allows for the preceding
                    # slices. The PLDs given out will be for the first
                    # slice.
                    otherInd = np.concatenate((np.arange(pld_idx), np.arange(pld_idx+1, scan.npld)))
                    params.PLD[otherInd, :] = np.tile(PLD[otherInd, np.newaxis], (1, PLDTryL)) + (slice*scan.slicedt)
                    params.PLD[pld_idx, :] = PLDTry + (slice*scan.slicedt)
                    params.t = params.tau + params.PLD
                    #print(params.PLD)
                    #print("PLDs shape", params.PLD.shape, np.mean(params.PLD))
                    #print("t shape", params.t.shape, np.mean(params.t))
                    #print("distWeightCurr shape", distWeightCurr.shape, np.mean(distWeightCurr))
                    #print("distInd shape", distInd.shape, np.mean(distInd))
                    
                    # FIXME refactor to class methods
                    if opttype.name == 'D-optimal':
                        # Calculate the determinant of the hessian
                        FIMCurr, _ = Hessian_DOptimal_analytical(params, scan, slice)

                        # Correct for small rounding errors leading to negative values
                        FIMCurr = np.abs(FIMCurr)
                        #print("FIMCUrr shape ", FIMCurr.shape, np.mean(FIMCurr))

                        # Take reciprocal to get the generalised variance
                        # Important to do this before taking the mean across slices
                        variance[:, distInd[:, slice], slice] = (1.0/FIMCurr)

                    elif opttype.name == 'L-optimal':
                        # Calculate the Trace of the covariance matrix
                        variance[:, distInd[:, slice], slice] = Hessian_LOptimal_analytical(params, opttype.A, scan, slice) # Here A is a weighting for the Trace
                    
                variance = variance * 6000 * 6000  # Change into (ml/100g/min)^2

                # Take mean of generalised variance across slices
                varianceMean = np.zeros((PLDTryL, att_dist.length))
                for ll in range(att_dist.length):
                    varianceMean[:, ll] = np.mean(variance[:, ll, distInd[ll, :]], 1)
                
                # Correct for artificial zeros (probably only happened when
                # I was using the pseudo-inverse, but leave here for
                # robustness)
                if np.any(varianceMean == 0):
                    print('WARNING: Some variances are zero. Setting to inf...')
                    varianceMean[varianceMean == 0] = np.inf

                #print("variance shape ", variance.shape, np.mean(variance))
                #print("varianceMean shape ", varianceMean.shape, np.mean(varianceMean))
                # Take mean of generalised variance across the BAT distribution
                cost = distWeightCurr * varianceMean                       # Weight by the ATT distribution
                cost[distWeightCurr == 0] = 0                              # To correct for 0*nan in att_dist.weight * CovCurr
                cost[np.isnan(cost)] = np.inf                                 # To correct for 0*inf in att_dist.weight * CovCurr
                #print("cost shape ", cost.shape, np.mean(cost))
                
                #sliceWeighting = sum(distInd,2)./allSliceL;               % Weighting based on how many slices contribute
                
                sliceWeighting = np.sum(distInd, 1)/np.arange(1, scan.slices+1)  # Weighting based on how many slices contribute
                #print("sliceWeighting shape ", sliceWeighting.shape, np.mean(sliceWeighting))

                costWeighted = cost * np.tile(sliceWeighting[np.newaxis, :], (PLDTryL, 1))   # Apply the slice weighting
                #print("costWeighted shape ", costWeighted.shape)
                costMean = np.mean(costWeighted, 1)                        # Weighted mean
                #print("costMean shape ", costMean.shape, np.mean(costMean))
                
                # Find the PLD that leads to the minimum generalised variance
                minVariance = np.min(costMean)
                jj = np.argmin(costMean)
                #print("Best is %i" % jj)

                PLD[pld_idx] = PLDTry[jj]
                #print("PLDs ", PLD)
                #print("oldPLD ", oldPLD)

                # Save the results from each step for visualising
                allCost.append(minVariance)
                allPLD.append(PLD)

                count += 1
                
        # If the PLDs stay constant, then exit
        #print("PLDs ", PLD)
        #print("oldPLD ", oldPLD)
        #print("minVar ", minVariance)
        if not np.any(PLD != oldPLD):
            log.write("\nFinished optimization after %i iters - PLDs unchanged\n" % count)
            break

    if (bestminVariance-minVariance) / bestminVariance > 1e-12:

        bestPLD = sorted(PLD)
        log.write('Optimal PLDs: %s\n' % bestPLD)

        bestminVariance = minVariance

        params.t = params.tau + PLD
        numAv, TotalTR = TRWeightingOrNAveFloor(params, scan, 0, 0)
        log.write('numAv = %i\n' % numAv)
        log.write('Scan time = %f\n' % (numAv*TotalTR))

    return bestPLD, numAv, bestminVariance

def TRWeightingOrNAveFloor(params, scan, tDim=0, slice=0):
    # I round the TR since there are occaisionally computational rounding
    # errors. I have used 5 decimal places to all dense BAT sampling, but I am
    # very unlikely to go finer than 0.001 density.
    fname = params.filename
    if fname == 'var_te_pCASL':            
        TEInd = max(params.t[:params.num_enc, 0, 0, 0, 0], [], 0)
        TR = np.squeeze(params.t[TEInd, :, :, :, :]) + scan.readout - (slice*scan.slicedt)
        # Multiply by the number of images that have to be acquired
        TotalTR = round(TR*(params.num_enc+1), 5)

    elif fname == 'var_te_pCASL_nPLD':
        TEInd = max(params.t[:params.num_enc, 0, 0, 0, 0], [], 0)
        TR = 0
        for ii in range(params.multiPLD):
            ind = TEInd + ii*params.num_enc
            TR = TR + np.squeeze(params.t[ind, :, :, :, :]) + scan.readout - (slice*scan.slicedt)

        # Multiply by the number of images that have to be acquired
        TotalTR = round(TR*(params.num_enc+1), 5)
    elif fname == 'var_multi_pCASL':
        TR = params.t + scan.readout - (slice*scan.slicedt)
        #print("TR shape ", TR.shape)
        TotalTR = np.round((2*(np.sum(TR, tDim))), 5)

    elif fname == 'look-locker':
        TR = params.t(end, 1) + (scan.readout/2) - (slice*scan.slicedt)
        TotalTR = np.round(TR*2, 5)

    #print(TotalTR.shape)
    return np.floor(np.round(np.squeeze(scan.duration/TotalTR), 5)), TotalTR

def Hessian_DOptimal_analytical(params, scan, slice):
    """
    function [detH, H] = Hessian_DOptimal_analytical(param, scanTime, slice, scan.slicedt, tReadout)
    
    This function calculates an approximation of the Hessian for optimising
    an experimental design. param is the struct containing the variables,
    
    Inputs:
            param    - A struct containing all of the parameters needed for the
                       Buxton standard CASL model (Buxton et al. MRM 1998)
            scanTime - The scan duration avalaible for the ASL protocol (units: seconds)
            slice    - The number of slices to optimise for (units: seconds)
            scan.slicedt  - The time to acquire each slice (the inter-slice
                       timing) (units: seconds)

    Written by Joseph Woods, FMRIB, Oxford, July 2017
    """
    # Get parameters into the correct matrix sizes
    params = HessianParseInputs(params)

    # Time should be put into the 1st dimension
    tDim = 0

    # Work out how many average we can fit in and divide by the noise SD
    TRWeight, _ = TRWeightingOrNAveFloor(params, scan, tDim, slice)
    TRWeight = TRWeight/(params.noise**2)
    #print("TRWeight shape ", TRWeight.shape, np.mean(TRWeight))

    df, dBAT = sensitivity(params)
    #print("df shape ", df.shape, np.mean(df))
    #print("dBAT shape ", dBAT.shape, np.mean(dBAT))

    # Form  the Hessian 
    # Hessian dimensions can be: 2 x 2 x PLD x Tau x BAT x f
    H = np.zeros(list(df.shape)[1:]+ [2, 2])
    H[..., 0, 0] = np.squeeze(np.sum(df*df, tDim))
    H[..., 0, 1] = np.squeeze(np.sum(df*dBAT, tDim))
    H[..., 1, 0] = np.squeeze(np.sum(dBAT*df, tDim))
    H[..., 1, 1] = np.squeeze(np.sum(dBAT*dBAT, tDim))
    #print("H shape ", H.shape)

    # Numpy determinant function batches over leading dimensions
    detH = np.linalg.det(H)
    detH = np.square(TRWeight) * detH
    #print("detH shape ", detH.shape)
    return detH, H

def HessianParseInputs(params):
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

def sensitivity(params):
    """
    function [df, dBAT] = sensitivity(param)
    
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

def covOEDNAveFloor(params, scan, slice, **kwargs):
    """
    function [covMatrix,H] = covOEDNAveFloor(param,scanTime,slice,inversion,slicedt,nAverage)
    %
    % This function calculates an approximation of the Hessian for optimising
    % an experimental design. param is the struct containing the variables,
    %
    % Written by Joseph Woods, FMRIB, Oxford, May 2017
    """
    #print("BAT shape1 ", params.bat.shape)
    tDim = 0
    params = HessianParseInputs(params)
    #print("BAT shape2 ", params.bat.shape)

    # Multiply TR by weights to get the correct weighted TotalTR

    TRWeight = kwargs.get("nAverage", None)
    if TRWeight is None:
        TRWeight, _ = TRWeightingOrNAveFloor(params, scan,tDim, slice)

    TRWeight = TRWeight/(params.noise**2)
    #print("TRWeight ", TRWeight)

    df, dBAT = sensitivity(params)
    #print("df shape ", df.shape, np.mean(df))
    #print("dBAT shape ", dBAT.shape, np.mean(dBAT))

    # Form  the Hessian 
    # Hessian dimensions can be: 2 x 2 x PLD x Tau x BAT x f
    H = np.zeros(list(df.shape)[1:]+ [2, 2])
    H[..., 0, 0] = TRWeight * (np.sum(df*df, tDim))
    H[..., 0, 1] = TRWeight * (np.sum(df*dBAT, tDim))
    H[..., 1, 0] = TRWeight * (np.sum(dBAT*df, tDim))
    H[..., 1, 1] = TRWeight * (np.sum(dBAT*dBAT, tDim))
    #print("H shape ", H.shape, np.mean(H[..., 0, 0]))
    #print("H shape ", H.shape, np.mean(H[..., 0, 1]))
    #print("H shape ", H.shape, np.mean(H[..., 1, 0]))
    #print("H shape ", H.shape, np.mean(H[..., 1, 1]))

    # Numpy determinant function batches over leading dimensions
    detH = np.linalg.det(H)
    sing = H[detH == 0]
    #print(np.count_nonzero(detH))
    #print(detH.size)
    invH = np.linalg.inv(H)
    #print("invH shape ", invH.shape)

    #covMatrixReshape = cell2mat(inverse);
    #covMatrix = reshape(covMatrixReshape,2,2,prod(Hsize(3:end)))
    covMatrix = invH

    # Change into (ml/100g/min)^
    covMatrix[..., 0, 0] = covMatrix[..., 0, 0] * 6000 * 6000;
    return covMatrix, H
