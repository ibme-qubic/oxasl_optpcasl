"""
PCASL optimisation script for Logan (Xin Zhang)

Written by Joseph Woods, FMRIB, Oxford, June 2019

This script dmonstrates how the optimisation functions are used, mainly
how to define the paramter inputs and how the plot the predicted CRLBs
SDs.
"""
import numpy as np
import matplotlib.pyplot as plt

import oxasl_optpcasl as opt

# Define the ASL parameters (as per Buxton et al. MRM 1998)
params = opt.ASLParams('var_multi_pCASL', 50.0/6000)

# ATT (BAT) distribution
att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)

# Details of the desired scan to optimize for
scan = opt.Scan(duration=300, readout=0.5)
        
# PLD limits and step size to search over
lims = opt.Limits(0.1, 3.0, 0.025)

# Number of PLDs to optimise for
nPLD = 6

# Type of optimisation
#opttype = opt.LOptimal([[1, 0],  [0, 0]])
opttype = opt.DOptimal()

# Run the optimisation
# Note: the output bestminVariance is not comparable between D-optimal and L-optimal
bestPLD, numAv, bestminVariance = opt.optimize(opttype, params, att_dist, scan, nPLD, lims)

## Example plotting of the CRLB SDs

params.pld = bestPLD
params.bat = att_dist.exclude_taper

covOptProtocol = np.zeros((scan.slices, len(params.bat), 2, 2))
for sliceN in range(scan.slices):
    params.t = params.tau + np.array(params.pld) + scan.slicedt*sliceN
    #numAv, TotalTR = TRWeightingOrNAveFloor(params, scanTime, 1, sliceN, scan.slicedt, tReadout)
    covOptProtocol[sliceN, ...], _ = opt.covOEDNAveFloor(params, scan, sliceN)

# Plot optimized protocol
plt.subplot(1, 2, 1)
plt.title("CBFSD")
plt.ylabel('SD (ml/100g/min)')
plt.xlabel("ATT (s)")
plt.ylim(0, 9)
plt.plot(params.bat, np.squeeze(np.sqrt(np.abs(covOptProtocol[..., 0, 0]))), 
         label="Optimised Protocol")
plt.legend()

plt.subplot(1, 2, 2)
plt.title("BATSD")
plt.ylabel("SD (s)")
plt.xlabel("ATT (s)")
plt.ylim(0, 0.25)
plt.plot(params.bat, np.squeeze(np.sqrt(np.abs(covOptProtocol[..., 1, 1]))), 
         label="Optimised Protocol")
plt.legend()

plt.show()
