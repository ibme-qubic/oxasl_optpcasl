#!/usr/bin/env python
"""
PCASL optimisation script for Logan (Xin Zhang)

This script dmonstrates how the optimisation functions are used, mainly
how to define the paramter inputs and how the plot the predicted CRLBs
SDs.
"""
import numpy as np
import matplotlib.pyplot as plt

import oxasl_optpcasl as opt

def main():
    """
    Simple example of optimizing the PLDs of a standard pCASL acquisition
    as per Buxton et al. MRM 1998
    """
    # Define the ASL parameters
    params = opt.ASLParams(opt.VAR_MULTI_PCASL, f=50.0/6000)

    # ATT (BAT) distribution
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)

    # Details of the desired scan to optimize for
    scan = opt.Scan(duration=300, npld=6, readout=0.5, slices=3, slicedt=0.0452)

    # PLD limits and step size to search over
    lims = opt.Limits(0.1, 3.0, 0.025)

    # Type of optimisation
    # Note: the output best_min_variance is not comparable between D-optimal and L-optimal
    #opttype = opt.LOptimal([[1, 0],  [0, 0]])
    opttype = opt.DOptimal()

    # Run the optimisation
    best_plds, num_av, best_min_variance = opt.optimize(opttype, params, att_dist, scan, lims)

    ## Example plotting of the CRLB SDs
    params.pld = best_plds
    params.bat = att_dist.exclude_taper

    cov_optimized = np.zeros((scan.slices, len(params.bat), 2, 2))
    for slice_idx in range(scan.slices):
        params.t = params.tau + np.array(params.pld) + scan.slicedt*slice_idx
        cov_optimized[slice_idx, ...], _ = opttype.cov_oedn_ave_floor(params, scan, slice_idx)

    # Plot optimized protocol
    plt.subplot(1, 2, 1)
    plt.title("CBFSD")
    plt.ylabel('SD (ml/100g/min)')
    plt.xlabel("ATT (s)")
    plt.ylim(0, 9)
    for slice_idx in range(scan.slices):
        plt.plot(params.bat, np.squeeze(np.sqrt(np.abs(cov_optimized[slice_idx, :, 0, 0]))),
                 label="Optimised Protocol - slice %i" % slice_idx)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.title("BATSD")
    plt.ylabel("SD (s)")
    plt.xlabel("ATT (s)")
    plt.ylim(0, 0.25)
    for slice_idx in range(scan.slices):
        plt.plot(params.bat, np.squeeze(np.sqrt(np.abs(cov_optimized[slice_idx, :, 1, 1]))),
                 label="Optimised Protocol - slice %i" % slice_idx)
    plt.legend()

    plt.show()

if __name__ == "__main__":
    main()
