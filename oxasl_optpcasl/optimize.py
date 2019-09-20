"""
OXASL_OPTPCASL - Functions for optimizing the PLDs of a multi-pld pCASL protocol

Copyright 2019 University of Oxford
"""
import sys
import copy

import numpy as np

from ._version import __version__
from .structures import VAR_MULTI_PCASL, VAR_TE_PCASL, LOOK_LOCKER, VAR_TE_PCASL_NPLD

def optimize(opttype, params, att_dist, scan, lims, log=sys.stdout):
    """
    Optimise the PLDs of a multi-pld PCASL protocol

    :param opttype: Optimization type object (L-optimal or D-optimal)
    :param params: ASLParams instance
    :param att_dist: BATDist instance giving ATT time distribution and weighting
    :param scan: Scan instance giving details of scan to optimize for
    :param lims: Limits instance, giving pld limiting values (max/min)
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
    best_min_variance = 1e99
    count = 0
    cost_history = []
    pld_history = []

    # Initialise the PLDs and check that the number of avearges > 0. Reduce
    # the pld range if not
    factor = 1.0/lims.step
    pld_subtract = 0
    tr_weight = 0
    while tr_weight < 1:

        # Initial sequence of PLDs spaced evenly between upper and lower bound
        pld = np.linspace(lims.lb, lims.ub - pld_subtract, scan.npld)
        pld = np.round(pld*factor) / factor
        params.pld = pld

        # Sequence of TIs corresponding to the PLDs. FIXME tau could be sequence?
        params.t = params.tau + params.pld
        tr_weight, _ = TRWeightingOrNAveFloor(params, scan, 0, 0)

        if tr_weight < 1:
            pld_subtract += 0.1

    # Optimisation loop
    while 1:
        # Each time through all of the samples, save them
        old_pld = np.copy(pld)

        for pld_idx in range(scan.npld):
            # Work out the set of available PLDs to search over. This is a series of
            # values between the pld before the current pld and the pld after, with a
            # step size of lims.step however special case for first and last PLDs
            # where we use the upper/lower bounds instead of the previous/next pld
            if pld_idx == 0:
                start, stop = lims.lb, pld[pld_idx+1]
            elif pld_idx == scan.npld-1:
                start, stop = pld[pld_idx-1], lims.ub
            else:
                start, stop = pld[pld_idx-1], pld[pld_idx+1]

            pld_trial = np.round(np.arange(start, stop+0.001, lims.step), 5)

            if len(pld_trial) > 0:
                # Current PLDs list is tiled with 1 column per trial pld
                params.pld = np.tile(pld[:, np.newaxis], (1, len(pld_trial)))

                # Distweight sequence is tiled with 1 row per trial pld
                dist_weight_curr = np.tile(att_dist.weight, (len(pld_trial), 1))

                # Flag matrix shape with 1 row per dist weight and 1 column per slice
                dist_ind = np.zeros((att_dist.length, scan.slices), dtype=np.bool)

                # Variance matrix shape (#trial PLDs, #distweights, #slices)
                variance = np.zeros((len(pld_trial), att_dist.length, scan.slices), dtype=np.float)

                for slice in range(scan.slices): # Loop over the slices

                    # Can only calculate the cost function for bat >
                    # shortest pld, otherwise we don't see the inflow and
                    # determinant blows up.
                    #
                    # NB Matlab uses min() here on lims.lb suggesting it could be a sequence?
                    min_pld_possible = lims.lb + slice*scan.slicedt

                    # For this slice dist_ind = True if corresponding bat is above the minimum pld
                    dist_ind[:, slice] = att_dist.dist > min_pld_possible
                    params.bat = att_dist.dist[dist_ind[:, slice]]

                    # TRWeightingOrNAveFloor only adjusts the readout time.
                    # Adjust the PLDs (and t) here. Add slice*scan.slicedt to
                    # the PLDs. If first slice, add 0s. If 12th slice, add
                    # 0.4972s.

                    # By specifying the slice, this calculation will ensure
                    # that the shortest pld allows for the preceding
                    # slices. The PLDs given out will be for the first
                    # slice.
                    other_ind = np.concatenate((np.arange(pld_idx), np.arange(pld_idx+1, scan.npld)))
                    params.pld[other_ind, :] = np.tile(pld[other_ind, np.newaxis], (1, len(pld_trial))) + (slice*scan.slicedt)
                    params.pld[pld_idx, :] = pld_trial + (slice*scan.slicedt)
                    params.t = params.tau + params.pld

                    variance[:, dist_ind[:, slice], slice] = opttype.hessian_var(params, scan, slice)

                variance = variance * 6000 * 6000  # Change into (ml/100g/min)^2

                # Take mean of generalised variance across slices
                variance_mean = np.zeros((len(pld_trial), att_dist.length))
                for idx in range(att_dist.length):
                    variance_mean[:, idx] = np.mean(variance[:, idx, dist_ind[idx, :]], 1)

                # Correct for artificial zeros (probably only happened when
                # I was using the pseudo-inverse, but leave here for
                # robustness)
                if np.any(variance_mean == 0):
                    log.write('WARNING: Some variances are zero. Setting to inf...\n')
                    variance_mean[variance_mean == 0] = np.inf

                # Take mean of generalised variance across the bat distribution
                cost = dist_weight_curr * variance_mean                       # Weight by the ATT distribution
                cost[dist_weight_curr == 0] = 0                              # To correct for 0*nan in att_dist.weight * CovCurr
                cost[np.isnan(cost)] = np.inf                                 # To correct for 0*inf in att_dist.weight * CovCurr

                #slice_weight = sum(dist_ind, 1)/allSliceL;               % Weighting based on how many slices contribute
                slice_weight = np.sum(dist_ind, 1)/np.arange(1, scan.slices+1)  # Weighting based on how many slices contribute

                cost_weighted = cost * np.tile(slice_weight[np.newaxis, :], (len(pld_trial), 1))   # Apply the slice weighting
                cost_mean = np.mean(cost_weighted, 1)                        # Weighted mean

                # Find the pld that leads to the minimum generalised variance
                min_variance = np.min(cost_mean)
                min_variance_idx = np.argmin(cost_mean)

                pld[pld_idx] = pld_trial[min_variance_idx]

                # Save the results from each step for visualising
                cost_history.append(min_variance)
                pld_history.append(pld)

                count += 1

        # If the PLDs stay constant, then exit
        if not np.any(pld != old_pld):
            log.write("\nFinished optimization after %i iters - PLDs unchanged\n" % count)
            break

    if (best_min_variance-min_variance) / best_min_variance > 1e-12:

        best_pld = sorted(pld)
        log.write('Optimal PLDs: %s\n' % best_pld)

        best_min_variance = min_variance

        params.t = params.tau + pld
        num_av, total_tr = TRWeightingOrNAveFloor(params, scan, 0, 0)
        log.write('num_av = %i\n' % num_av)
        log.write('Scan time = %f\n' % (num_av*total_tr))

    return best_pld, num_av, best_min_variance

def TRWeightingOrNAveFloor(params, scan, time_dim=0, slice=0):
    # I round the tr since there are occaisionally computational rounding
    # errors. I have used 5 decimal places to all dense bat sampling, but I am
    # very unlikely to go finer than 0.001 density.
    if params.asltype == VAR_TE_PCASL:
        te_ind = max(params.t[:params.num_enc, 0, 0, 0, 0], [], 0)
        tr = np.squeeze(params.t[te_ind, :, :, :, :]) + scan.readout - (slice*scan.slicedt)
        # Multiply by the number of images that have to be acquired
        total_tr = round(tr*(params.num_enc+1), 5)

    elif params.asltype == VAR_TE_PCASL_NPLD:
        te_ind = max(params.t[:params.num_enc, 0, 0, 0, 0], [], 0)
        tr = 0
        for ii in range(params.multiPLD):
            ind = te_ind + ii*params.num_enc
            tr = tr + np.squeeze(params.t[ind, :, :, :, :]) + scan.readout - (slice*scan.slicedt)

        # Multiply by the number of images that have to be acquired
        total_tr = round(tr*(params.num_enc+1), 5)
    elif params.asltype == VAR_MULTI_PCASL:
        tr = params.t + scan.readout - (slice*scan.slicedt)
        total_tr = np.round((2*(np.sum(tr, time_dim))), 5)

    elif params.asltype == LOOK_LOCKER:
        tr = params.t(end, 1) + (scan.readout/2) - (slice*scan.slicedt)
        total_tr = np.round(tr*2, 5)
    else:
        raise ValueError("Unrecognized ASL type: %s" % params.asltype)

    return np.floor(np.round(np.squeeze(scan.duration/total_tr), 5)), total_tr

class OptimizationMethod(object):
    """
    Optimization method base class
    
    Defines common methods used by L-optimal and D-optimal methods
    """

    def resize_inputs(self, params):
        """
        Returns a modified version of params with t, tau, f and bat
        all set to a common shape
        """
        params = copy.copy(params)
        params.tau = np.full(params.t.shape, params.tau)

        # Tile t and tau so third dimension has size of bat
        params.t = np.squeeze(np.tile(params.t[..., np.newaxis], [1, 1, params.bat.shape[0]]))
        params.tau = np.squeeze(np.tile(params.tau[..., np.newaxis], [1, 1, params.bat.shape[0]]))

        # Make f and bat the same dimensions as t
        # f is constant, bat is 1D (which should be preserved in the 3rd dimension)
        params.f = np.full(params.t.shape, params.f)
        if len(params.t.shape) == 2:
            params.bat = np.tile(params.bat[np.newaxis, :], list(params.t.shape)[:-1] + [1])
        else:
            params.bat = np.tile(params.bat[np.newaxis, np.newaxis, :], list(params.t.shape)[:-1] + [1])

        # FIXME deleted MATLAB error handling here for 1 pld - test?
        return params

    def sensitivity(self, params):
        """
        This function calculates the sensitivty functions of the Buxtion CASL
        model (Buxton et al. MRM 1998) and are given in Woods et al. MRM 2019.
        The bat sensitivty function assumes a fixed outflow in t1_prime in order
        to simplify the equations.
        """
        # This sensitivty function includes an assumed fixed t1_prime, so we fix the
        # f in t1_prime to a sensible value.
        f_fixed = 50.0/6000

        t1_prime = 1.0/((1.0/params.t1t) + (f_fixed/params.lam))
        M = 2*params.m0b * params.alpha * t1_prime * np.exp(-params.bat/params.t1b)

        # Initialise
        df = np.zeros(params.t.shape, dtype=np.float32)
        dbat = np.zeros(params.t.shape, dtype=np.float32)

        # for t between deltaT and tau plus deltaT
        t_region = np.logical_and(params.t > params.bat, params.t <= (params.tau + params.bat))
        t = params.t[t_region]
        # if sum(size(params.bat)>1) > 1 # Check whether multiple values are being used
        if np.any(np.array(params.bat.shape) > 1):
            bat = params.bat[t_region]
            f = params.f[t_region]
            m_temp = M[t_region]
        else:
            bat = params.bat
            f = params.f
            m_temp = M

        t1_prime_temp = t1_prime

        df[t_region] = m_temp * (1 - np.exp((bat-t) / t1_prime_temp))
        dbat[t_region] = m_temp * f * ((-1.0/params.t1b) - np.exp((bat-t) / t1_prime_temp) * ((1.0/t1_prime_temp) - (1.0/params.t1b)))

        # for t greater than tau plus deltaT
        t_region = params.t > (params.tau+params.bat)
        t = params.t[t_region]
        if np.any(np.array(params.tau.shape) > 1):
            tau = params.tau[t_region]
        else:
            tau = params.tau

        # if sum(size(params.bat)>1) > 1
        if np.any(np.array(params.bat.shape) > 1):
            bat = params.bat[t_region]
            f = params.f[t_region]
            m_temp = M[t_region]
        else:
            bat = params.bat
            f = params.f
            m_temp = M
        t1_prime_temp = t1_prime

        df[t_region] = m_temp * np.exp((-t+tau+bat) / t1_prime_temp) * (1-np.exp(-tau/t1_prime_temp))
        dbat[t_region] = m_temp * f * (1-np.exp(-tau/t1_prime_temp)) * np.exp((bat+tau-t)/t1_prime_temp) * ((1.0/t1_prime_temp)-(1.0/params.t1b))

        return df, dbat

    def hessian(self, params, scan, slice, time_dim=0, **kwargs):
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
        tr_weight = kwargs.get("nAverage", None)
        if tr_weight is None:
            tr_weight, _ = TRWeightingOrNAveFloor(params, scan, time_dim, slice)

        tr_weight = tr_weight/(params.noise**2)
        df, dbat = self.sensitivity(params)

        # Form the Hessian
        # Hessian dimensions can be: 2 x 2 x pld x Tau x bat x f
        hess = np.zeros(list(df.shape)[1:]+ [2, 2])
        hess[..., 0, 0] = tr_weight * np.squeeze(np.sum(df*df, time_dim))
        hess[..., 0, 1] = tr_weight * np.squeeze(np.sum(df*dbat, time_dim))
        hess[..., 1, 0] = tr_weight * np.squeeze(np.sum(dbat*df, time_dim))
        hess[..., 1, 1] = tr_weight * np.squeeze(np.sum(dbat*dbat, time_dim))
        return hess

    def cov_oedn_ave_floor(self, params, scan, slice, **kwargs):
        hess = self.hessian(params, scan, slice, **kwargs)

        # Numpy determinant function batches over leading dimensions
        cov = np.linalg.inv(hess)

        # Change into (ml/100g/min)
        cov[..., 0, 0] = cov[..., 0, 0] * 6000 * 6000
        return cov, hess

class LOptimal(OptimizationMethod):
    """
    Choose CBF variance to minimise
    """
    def __init__(self, A):
        self.A = A
        self.name = 'L-optimal'

    def hessian_var(self, params, scan, slice):
        hess = self.hessian(params, scan, slice)

        # Numpy matrix ops batch over leading dimensions
        cov = np.abs(np.matmul(self.A, np.linalg.inv(hess)))

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

    def hessian_var(self, params, scan, slice):
        hess = self.hessian(params, scan, slice)

        # Numpy determinant function batches over leading dimensions
        det_h = np.linalg.det(hess)

        return 1.0/np.abs(det_h)
    