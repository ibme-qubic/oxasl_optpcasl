"""
Test cases for PLD optimization

Ground truth is generally taken from the MATLAB code
"""
import numpy as np

import oxasl_optpcasl as opt

def test_doptimal():
    """ Test D-optimal optimization method """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.DOptimal()

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.2, 0.7, 0.725, 1.55, 1.875, 2.075])

def test_loptimal():
    """ Test L-optimal optimization method """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.LOptimal(A=[[1, 0], [0, 0]])

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.2, 1.175, 1.8, 2.025, 2.1, 2.1])
