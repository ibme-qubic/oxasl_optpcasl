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
                    
                    variance[:, distInd[:, slice], slice] = opttype.HessianCov(params, scan, slice)

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
