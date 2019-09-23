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

def test_loptimal_cbf():
    """ Test L-optimal optimization method for CBF """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.LOptimal(A=[[1, 0], [0, 0]])

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.2, 1.175, 1.8, 2.025, 2.1, 2.1])

def test_loptimal_att():
    """ Test L-optimal optimization method for ATT """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.LOptimal(A=[[0, 0], [0, 1]])

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.1, 0.475, 0.7, 1.025, 1.725, 2.1])

def test_doptimal_multislice():
    """ Test D-optimal optimization method with 2D readout """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5, slices=10, slicedt=0.0452)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.DOptimal()

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.1, 0.575, 0.725, 1.4, 1.75, 2.025])

def test_loptimal_cbf_multislice():
    """ Test L-optimal optimization method for CBF with 2D readout """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5, slices=10, slicedt=0.0452)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.LOptimal(A=[[1, 0], [0, 0]])

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.1, 1.025, 1.625, 1.8, 1.95, 2.1])

def test_loptimal_att_multislice():
    """ Test L-optimal optimization method for ATT with 2D readout """
    params = opt.ASLParams('var_multi_pCASL', f=50.0/6000)
    att_dist = opt.BATDist(0.2, 2.1, 0.001, 0.3)
    scan = opt.Scan(duration=300, npld=6, readout=0.5, slices=10, slicedt=0.0452)
    lims = opt.Limits(0.1, 3.0, 0.025)
    opttype = opt.LOptimal(A=[[0, 0], [0, 1]])

    best_plds, num_av, best_min_var = opt.optimize(opttype, params, att_dist, scan, lims)
    assert np.allclose(best_plds, [0.1, 0.375, 0.7, 1.075, 1.65, 2.000])
