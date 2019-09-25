"""
OXASL_OPTPCASL - Functions for optimizing the PLDs of a multi-pld pCASL protocol

Copyright 2019 University of Oxford
"""
import sys
import copy

import numpy as np

from ._version import __version__
from .structures import VAR_MULTI_PCASL, VAR_TE_PCASL, LOOK_LOCKER, VAR_TE_PCASL_NPLD

class OptimizationOutput(object):
    """
    Stores the output of the optimization process as attributes
    """

    def __init__(self):
        self.cost_history = []
        self.pld_history = []
        self.best_min_variance = 1e99

class Optimizer(object):
    """
    Optimization base class
    
    Defines the basic algorithm for optimization. The L-optimal and D-optimal subclasses
    implement the particular cost functions for these methods
    """

    def __init__(self, params, scan, att_dist, lims, log=sys.stdout):
        """
        :param params: ASLParams instance
        :param scan: ASLScan instance giving details of scan to optimize for
        :param att_dist: ATTDist instance giving ATT time distribution and weighting
        :param lims: Limits instance, giving pld limiting values (max/min)
        :param log: Stream-like object for logging output
        """
        self.params = params
        self.scan = scan
        self.att_dist = att_dist
        self.lims = lims
        self.log = log

    def optimize(self):
        """
        Optimise the PLDs of a multi-pld PCASL protocol
        """
        self.log.write("Optimizing PLDs for %s\n" % self.scan)
        self.log.write("PLD search limits: %s\n" % self.lims)
        self.log.write("Optimizing for %i PLDs\n" % self.scan.npld)
        self.log.write("%s\n" % str(self.att_dist))
        self.log.write("Optimization method: %s\n" % self.name)

        # Initialise
        output = OptimizationOutput()
        count = 0

        # Initialise the PLDs and check that the number of avearges > 0. Reduce
        # the pld range if not
        factor = 1.0/self.lims.step
        pld_subtract = 0
        tr_weight = 0
        while tr_weight < 1:

            # Initial sequence of PLDs spaced evenly between upper and lower bound
            output.plds = np.linspace(self.lims.lb, self.lims.ub - pld_subtract, self.scan.npld)
            output.plds = np.round(output.plds*factor) / factor

            # Sequence of TIs corresponding to the PLDs. FIXME tau could be sequence?
            tr_weight, _ = self.TRWeightingOrNAveFloor(self.params.tau + output.plds, 0)
            if tr_weight < 1:
                pld_subtract += 0.1

        # Optimisation loop
        while 1:
            # Each time through all of the samples, save them
            old_pld = np.copy(output.plds)

            for pld_idx in range(self.scan.npld):
                # Work out the set of available PLDs to search over. This is a series of
                # values between the pld before the current pld and the pld after, with a
                # step size of self.lims.step however special case for first and last PLDs
                # where we use the upper/lower bounds instead of the previous/next pld
                if pld_idx == 0:
                    start, stop = self.lims.lb, output.plds[pld_idx+1]
                elif pld_idx == self.scan.npld-1:
                    start, stop = output.plds[pld_idx-1], self.lims.ub
                else:
                    start, stop = output.plds[pld_idx-1], output.plds[pld_idx+1]

                trial_values = np.round(np.arange(start, stop+0.001, self.lims.step), 5)

                if len(trial_values) > 0:
                    # Current PLDs list is tiled with 1 column per trial plds
                    trial_plds = np.tile(output.plds[:, np.newaxis], (1, len(trial_values)))

                    # Distweight sequence is tiled with 1 row per trial pld
                    dist_weight_curr = np.tile(self.att_dist.weight, (len(trial_values), 1))

                    # Flag matrix shape with 1 row per dist weight and 1 column per slice
                    dist_ind = np.zeros((self.att_dist.length, self.scan.slices), dtype=np.bool)

                    # Variance matrix shape (#trial PLDs, #distweights, #slices)
                    variance = np.zeros((len(trial_values), self.att_dist.length, self.scan.slices), dtype=np.float)

                    for slice in range(self.scan.slices): # Loop over the slices

                        # Can only calculate the cost function for att >
                        # shortest pld, otherwise we don't see the inflow and
                        # determinant blows up.
                        #
                        # NB Matlab uses min() here on self.lims.lb suggesting it could be a sequence?
                        min_pld_possible = self.lims.lb + slice*self.scan.slicedt

                        # For this slice dist_ind = True if corresponding att is above the minimum pld
                        dist_ind[:, slice] = self.att_dist.dist > min_pld_possible
                        att = self.att_dist.dist[dist_ind[:, slice]]

                        # TRWeightingOrNAveFloor only adjusts the readout time.
                        # Adjust the PLDs (and t) here. Add slice*self.scan.slicedt to the PLDs. 

                        # By specifying the slice, this calculation will ensure
                        # that the shortest plds allows for the preceding
                        # slices. The PLDs given out will be for the first
                        # slice.
                        other_ind = np.concatenate((np.arange(pld_idx), np.arange(pld_idx+1, self.scan.npld)))
                        trial_plds[other_ind, :] = np.tile(output.plds[other_ind, np.newaxis], (1, len(trial_values))) + (slice*self.scan.slicedt)
                        trial_plds[pld_idx, :] = trial_values + (slice*self.scan.slicedt)

                        variance[:, dist_ind[:, slice], slice] = self.hessian_var(self.params.tau + trial_plds, att, slice)

                    # Take mean of generalised variance across slices
                    variance_mean = np.zeros((len(trial_values), self.att_dist.length))
                    for idx in range(self.att_dist.length):
                        variance_mean[:, idx] = np.mean(variance[:, idx, dist_ind[idx, :]], 1)

                    # Correct for artificial zeros (probably only happened when
                    # I was using the pseudo-inverse, but leave here for
                    # robustness)
                    if np.any(variance_mean == 0):
                        self.log.write('WARNING: Some variances are zero. Setting to inf...\n')
                        variance_mean[variance_mean == 0] = np.inf

                    # Take mean of generalised variance across the att distribution
                    cost = dist_weight_curr * variance_mean                       # Weight by the ATT distribution
                    cost[dist_weight_curr == 0] = 0                              # To correct for 0*nan in self.att_dist.weight * CovCurr
                    cost[np.isnan(cost)] = np.inf                                 # To correct for 0*inf in self.att_dist.weight * CovCurr

                    #slice_weight = sum(dist_ind, 1)/allSliceL;               % Weighting based on how many slices contribute
                    slice_weight = np.sum(dist_ind, 1)/self.scan.slices            # Weighting based on how many slices contribute

                    cost_weighted = cost * np.tile(slice_weight[np.newaxis, :], (len(trial_values), 1))   # Apply the slice weighting
                    cost_mean = np.mean(cost_weighted, 1)                        # Weighted mean

                    # Find the pld that leads to the minimum generalised variance
                    min_variance = np.min(cost_mean)
                    min_variance_idx = np.argmin(cost_mean)

                    output.plds[pld_idx] = trial_values[min_variance_idx]

                    # Save the results from each step for visualising
                    output.cost_history.append(min_variance)
                    output.pld_history.append(output.plds)

                    count += 1

            # If the PLDs stay constant, then exit
            if not np.any(output.plds != old_pld):
                self.log.write("\nFinished optimization after %i iters - PLDs unchanged\n" % count)
                break

        if (output.best_min_variance - min_variance) / output.best_min_variance > 1e-12:
            output.plds = np.array(sorted(output.plds))
            output.times = self.params.tau + output.plds
            output.best_min_variance = min_variance
            output.num_av, output.total_tr = self.TRWeightingOrNAveFloor(output.times, 0)
            output.scan_time = output.total_tr * output.num_av
            output.att = self.att_dist.exclude_taper
            output.cov_optimized = np.zeros((self.scan.slices, len(output.att), 2, 2))
            for slice_idx in range(self.scan.slices):
                times = output.times + self.scan.slicedt*slice_idx
                output.cov_optimized[slice_idx, ...] = self.cov_lower_bound(times, output.att, slice_idx)

        return output

    def TRWeightingOrNAveFloor(self, times, slice=0):
        # I round the tr since there are occaisionally computational rounding
        # errors. I have used 5 decimal places to all dense att sampling, but I am
        # very unlikely to go finer than 0.001 density.
        #
        # FIXME the num_enc parameter is not defined so the first two scan types will
        # not work presently. This issue is also in the source MATLAB code. So we will
        # disable these scan types for now until we know what this should be
        #if self.scan.asltype == VAR_TE_PCASL:
        #    te_ind = max(times[:self.params.num_enc, 0, 0, 0, 0], [], 0)
        #    tr = np.squeeze(times[te_ind, :, :, :, :]) + self.scan.readout - (slice*self.scan.slicedt)
        #    # Multiply by the number of images that have to be acquired
        #    total_tr = round(tr*(self.params.num_enc+1), 5)
        #
        #elif self.scan.asltype == VAR_TE_PCASL_NPLD:
        #    te_ind = max(times[:self.params.num_enc, 0, 0, 0, 0], [], 0)
        #    tr = 0
        #    for ii in range(self.params.multiPLD):
        #        ind = te_ind + ii*self.params.num_enc
        #        tr = tr + np.squeeze(times[ind, :, :, :, :]) + self.scan.readout - (slice*self.scan.slicedt)
        #
        #    # Multiply by the number of images that have to be acquired
        #    total_tr = round(tr*(self.params.num_enc+1), 5)
        if self.scan.asltype == VAR_MULTI_PCASL:
            tr = times + self.scan.readout - (slice*self.scan.slicedt)
            total_tr = np.squeeze(np.round((2*(np.sum(tr, 0))), 5))
        elif self.scan.asltype == LOOK_LOCKER:
            tr = times[-1] + (self.scan.readout/2) - (slice*self.scan.slicedt)
            total_tr = np.round(tr*2, 5)
        else:
            raise ValueError("Unrecognized ASL type: %s" % self.scan.asltype)

        return np.atleast_1d(np.floor(np.round(np.squeeze(self.scan.duration/total_tr), 5))), total_tr

    def sensitivity(self, times, att):
        """
        This function calculates the sensitivty functions of the Buxtion CASL
        model (Buxton et al. MRM 1998) and are given in Woods et al. MRM 2019.
        The att sensitivty function assumes a fixed outflow in t1_prime in order
        to simplify the equations.
        """
        # This sensitivty function includes an assumed fixed t1_prime, so we fix the
        # f in t1_prime to a sensible value.
        f_fixed = 50.0/6000

        t1_prime = 1.0/((1.0/self.params.t1t) + (f_fixed/self.params.lam))
        M = 2*self.params.m0b * self.params.alpha * t1_prime * np.exp(-att / self.params.t1b)

        # Initialise
        times = np.repeat(times[..., np.newaxis], len(att), -1)
        att = att[np.newaxis, :]
        if times.ndim > 2:
            att = att[np.newaxis, :]
        df = np.zeros(times.shape, dtype=np.float32)
        datt = np.zeros(times.shape, dtype=np.float32)

        # for t between deltaT and tau plus deltaT
        t_during = np.logical_and(times > att, times <= (self.params.tau + att))
        df_during = M * (1 - np.exp((att - times) / t1_prime))
        datt_during = M * self.params.f * ((-1.0/self.params.t1b) - np.exp((att - times) / t1_prime) * ((1.0/t1_prime) - (1.0/self.params.t1b)))
        df[t_during] = df_during[t_during]
        datt[t_during] = datt_during[t_during]

        # for t greater than tau plus deltaT
        t_after = times > (self.params.tau + att)
        df_after = M * np.exp((self.params.tau + att - times) / t1_prime) * (1 - np.exp(-self.params.tau/t1_prime))
        datt_after = M * self.params.f * (1 - np.exp(-self.params.tau/t1_prime)) * np.exp((self.params.tau + att - times)/t1_prime) * (1.0/t1_prime - 1.0/self.params.t1b)
        df[t_after] = df_after[t_after]
        datt[t_after] = datt_after[t_after]

        return df, datt

    def hessian(self, times, att, slice, **kwargs):
        """
        This function calculates an approximation of the Hessian for optimising
        an experimental design. param is the struct containing the variables,

        :param times: 1D or 2D array of PLDs, if 2D calculate for each set of PLDs
        :param att: 1D array of ATT times to calculate Hessian at
        :param slice: Slice index
        """
        # Work out how many averages we can fit in and divide by the noise SD
        tr_weight, _ = self.TRWeightingOrNAveFloor(times, slice)
        tr_weight = tr_weight/(self.params.noise**2)

        # Calculate derivatives of ASL kinetic model
        df, datt = self.sensitivity(times, att)

        # Form the Hessian - dimensions are: num PLD sets x num ATTs x 2 x 2
        # Note Numpy matrix functions batch over leading dimensions
        hess = np.zeros(list(df.shape)[1:]+ [2, 2])
        tr_weight = tr_weight[:, np.newaxis]
        hess[..., 0, 0] = tr_weight * np.squeeze(np.sum(df*df, 0))
        hess[..., 0, 1] = tr_weight * np.squeeze(np.sum(df*datt, 0))
        hess[..., 1, 0] = tr_weight * np.squeeze(np.sum(datt*df, 0))
        hess[..., 1, 1] = tr_weight * np.squeeze(np.sum(datt*datt, 0))
        return hess

    def cov_lower_bound(self, times, att, slice, **kwargs):
        """
        :return: Lower bound on covariance matrix
        """
        hess = self.hessian(times, att, slice, **kwargs)
        cov = np.linalg.inv(hess)

        # Correct for inf*0 errors in A*inverse
        cov[np.isnan(cov)] = np.inf

        # Change into (ml/100g/min)
        cov[..., 0, 0] = cov[..., 0, 0] * 6000 * 6000
        cov[..., 0, 1] = cov[..., 0, 1] * 6000
        cov[..., 1, 0] = cov[..., 1, 0] * 6000
        return cov

class LOptimal(Optimizer):
    """
    Choose CBF or ATT variance to minimise
    """
    def __init__(self, A, *args, **kwargs):
        Optimizer.__init__(self, *args, **kwargs)
        self.A = A
        self.name = 'L-optimal'

    def hessian_var(self, times, att, slice):
        cov = np.abs(np.matmul(self.A, self.cov_lower_bound(times, att, slice)))

        # Force trace function to batch across leading dimensions
        return np.trace(cov, axis1=-1, axis2=-2)

class DOptimal(Optimizer):
    """
    Optimize for both CBF and ATT variance
    """
    def __init__(self, *args, **kwargs):
        Optimizer.__init__(self, *args, **kwargs)
        self.name = 'D-optimal'

    def hessian_var(self, times, att, slice):
        hess = self.hessian(times, att, slice)
        det_h = np.linalg.det(hess)
        return 1.0/np.abs(det_h)
