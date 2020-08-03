"""
OXASL_OPTPCASL - ASL Scan protocols

Copyright 2019 University of Oxford
"""
import numpy as np

class Protocol(object):
    """
    A scan protocol that can be optimized

    A protocol has a set of variable parameters, can generate 'trial' values for
    each parameter, and can calculate the cost of a set of parameters
    """

    def __init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims, ld_lims=None):
        self.kinetic_model = kinetic_model
        self.cost_model = cost
        self.scan_params = scan_params
        self.att_dist = att_dist
        self.pld_lims = pld_lims
        self.ld_lims = ld_lims

    def initial_params(self):
        """
        Get the initial parameter set

        Parameters may be a set of PLDs, or possibly with additional labelling durations.
        It's up to the scan type to use the same parameter format in other methods
        e.g. ``cost``

        :return: Initial parameters [NParams]
        """
        raise NotImplementedError()

    def trial_params(self, params, idx):
        """
        Get a set of trial parameter values by varying one of the parameters

        :param params: Current initial param values [NParams]
        :param idx: Index of parameter to vary

        :return: Trial parameter values [Trials, NParams]
        """
        trial_values = self.trial_param_values(params, idx)
        trial_params = np.tile(params[np.newaxis, :], (len(trial_values), 1))
        trial_params[:, idx] = trial_values
        return trial_params

    def trial_param_values(self, params, idx):
        """
        Get a set of trial values for a single parameter

        :param params: Current initial param values [NParams]
        :param idx: Index of parameter to vary

        :return: Trial values of the parameter to vary [Trials]
        """
        raise NotImplementedError()

    def name_params(self, params):
        """
        Unpack a set of parameter values into named groups (e.g. PLDs, LDs)

        :param params: Parameter values [NParams]
        :return: Mapping from parameter name to value or array of values
        """
        raise NotImplementedError()

    def cost(self, params):
        """
        Get the cost for a set of parameters

        :param params: Parameter set [NTrials, NParams] or [NParams]
        :return Cost [NTrials] or scalar
        """
        raise NotImplementedError()

    def repeats_total_tr(self, params):
        """
        Get the number of averages possible and the total TR.

        :param params: Set of trial parameters [Trials, NParams] or [NParams]
        :return Tuple of number of averages, Total TR [Trials] or scalars
        """
        raise NotImplementedError()

    def param_bounds(self):
        """
        Get the max and min values for all parameters

        :return: Sequence of tuples (lower bound, upper bound), one for each parameter
        """
        raise NotImplementedError()

class PCASLProtocol(Protocol):
    """
    A PCASL protocol which may have multiple PLDs and independent
    labelling durations associated with them.
    """

    def __init__(self, *args, **kwargs):
        Protocol.__init__(self, *args, **kwargs)

        # Slice time offset [NSlices]
        self.slicedt = np.arange(self.scan_params.nslices, dtype=np.float) * self.scan_params.slicedt

        # We can only calculate the cost function for ATT > shortest PLD, otherwise we don't
        # see the inflow and the determinant blows up. So we modify the ATT distribution
        # weight to be zero at ATTs which are not relevant to a slice.
        # Note that we do not renormalize the ATT distribution weights because we want slices
        # to contribute more when they have more relevant ATTs
        min_pld_possible = self.pld_lims.lb + self.slicedt
        relevant_atts = self.att_dist.atts[np.newaxis, ...] > min_pld_possible[..., np.newaxis]
        self.att_weights = self.att_dist.weight * relevant_atts

    def cost(self, params):
        # Hessian matrix for sensitivity of kinetic model [NTrials, NSlices, NATTs, 2, 2]
        hessian = self._hessian(params)

        # Cost at each time point and for each ATT in distribution [NTrials, NSlices, NATTs]
        cost = self.cost_model.cost(hessian)
        #print("cost", cost.shape)
        #print(cost)

        # Weight each component of the cost according to ATT distribution weight
        # This also takes into account which ATTs are relevant to which slices
        att_weights = self.att_weights
        if params.ndim == 2:
            att_weights = np.tile(self.att_weights[np.newaxis, ...], (params.shape[0], 1, 1))

        cost *= att_weights
        cost[att_weights == 0] = 0      # To correct for 0*nan
        cost[np.isnan(cost)] = np.inf   # To correct for 0*inf
        #print("wcost", cost.shape)
        #print(cost)

        # Take mean of cost across slices and ATTs
        return np.mean(cost, axis=(-1, -2))

    def repeats_total_tr(self, params):
        # Allow for label/control image at each time point
        lds, plds = self._timings(params)
        tr = lds + plds + self.scan_params.readout
        total_tr = 2*np.sum(tr, axis=-1)

        # I round the tr since there are occasionally computational rounding
        # errors. I have used 5 decimal places to allow dense att sampling, but I am
        # very unlikely to go finer than 0.001 density.
        return np.floor(self.scan_params.duration/total_tr), np.round(total_tr, 5)

    def _timings(self, params):
        """
        Get the effective labelling duration and PLDs for a set of trial params

        :param params: Parameter values [NTrials, NParams] or [NParams]
        :return Tuple of labelling duration and PLDs, each [NTrials, NPLDs] or [NPLDs]
        """
        raise NotImplementedError()

    def _hessian(self, params):
        """
        FIXME noise scaling for TE data

        :param params: Parameters [NTrials, NParams] or [NParams]

        :return: Hessian matrices of second derivatives wrt cbf
                 and att [NTrials, NSlices, NATTs, 2, 2] or [NSlices, NATTs, 2, 2]
        """
        # Time points to evaluate the sensitivity of the kinetic model at
        # [NTrials, NPLDs, NSlices]
        lds, plds = self._timings(params)
        lds = np.repeat(lds[..., np.newaxis], len(self.slicedt), axis=-1)
        plds = plds[..., np.newaxis]
        slicedt = self.slicedt[np.newaxis, ...]
        while slicedt.ndim != plds.ndim:
            slicedt = slicedt[np.newaxis, ...]
        times = lds + plds + slicedt

        # Work out how many averages we can fit in and divide by the noise SD
        num_repeats, _tr = self.repeats_total_tr(params)
        num_repeats = num_repeats/(self.scan_params.noise**2)

        # Calculate derivatives of ASL kinetic model [NTrials, NPLDs, NSlices, NATTs]
        att = self.att_dist.atts
        df, datt = self.kinetic_model.sensitivity(lds.flatten(), times.flatten(), att)
        df = df.reshape(list(times.shape) + [len(att)])
        datt = datt.reshape(list(times.shape) + [len(att)])
        #print("df, datt", df.shape)
        #print(df)
        #print(datt)
        #print(times)

        # Form the Hessian [Trials, Slices, ATTs, 2, 2], summed over PLDs
        # Note Numpy matrix functions batch over leading dimensions
        if params.ndim == 1:
            hess = np.zeros([df.shape[1], df.shape[2], 2, 2])
            pld_idx = 0
        else:
            hess = np.zeros([df.shape[0], df.shape[2], df.shape[3], 2, 2])
            pld_idx = 1
        num_repeats = np.atleast_1d(num_repeats)
        while num_repeats.ndim < hess.ndim - 2:
            num_repeats = num_repeats[:, np.newaxis]
        hess[..., 0, 0] = num_repeats * np.sum(df*df, axis=pld_idx)
        hess[..., 0, 1] = num_repeats * np.sum(df*datt, axis=pld_idx)
        hess[..., 1, 0] = num_repeats * np.sum(datt*df, axis=pld_idx)
        hess[..., 1, 1] = num_repeats * np.sum(datt*datt, axis=pld_idx)
        return hess

class FixedLDPCASLProtocol(PCASLProtocol):
    """
    PCASL protocol with a single fixed labelling duration
    """
    def __str__(self):
        return "Multi-PLD PCASL protocol with fixed label duration"

    def name_params(self, params):
        return {
            "plds" : params
        }

    def initial_params(self):
        # Parameters in this case are a set of PLDs
        if self.scan_params.plds is not None:
            return self.scan_params.plds
        else:
            factor = 1.0/self.pld_lims.step
            max_pld = self.pld_lims.ub
            while 1:

                # Initial sequence of PLDs spaced evenly between upper and lower bound
                plds = np.linspace(self.pld_lims.lb, max_pld, self.scan_params.npld)
                plds = np.round(plds*factor) / factor

                # See how long the TR is - if it is larger than the maximum, reduce the maximum PLD
                total_tr = 2*np.sum(self.scan_params.ld + plds + self.scan_params.readout)
                if total_tr <= self.scan_params.duration:
                    break
                max_pld -= 0.1

            return plds

    def trial_param_values(self, params, idx):
        # For the first and last PLDs we use the upper/lower bounds instead of the
        # previous/next pld.
        start, stop = self.pld_lims.lb, self.pld_lims.ub
        if idx > 0:
            start = params[idx-1]
        if idx < self.scan_params.npld-1:
            stop = params[idx+1]

        return np.round(np.arange(start, stop+0.001, self.pld_lims.step), 5)

    def param_bounds(self):
        return [
            (self.pld_lims.lb, self.pld_lims.ub)
            for idx in range(self.scan_params.npld)
        ]

    def _timings(self, params):
        return np.full(params.shape, self.scan_params.ld), params

class MultiPLDPcaslVarLD(FixedLDPCASLProtocol):
    """
    PCASL protocol with multiple PLDs and single variable LD
    """

    def __str__(self):
        return "Multi-PLD PCASL protocol with single variable label duration"

    def name_params(self, params):
        return {
            "plds" : params[:self.scan_params.npld],
            "lds" : params[self.scan_params.npld],
        }

    def initial_params(self):
        plds = FixedLDPCASLProtocol.initial_params(self)

        # Add variable label duration to parameters
        return np.array(list(plds) + [self.scan_params.ld])

    def trial_params_values(self, params, idx):
        if idx < self.scan_params.npld:
            # We are varying a PLD. Get the trial values from the base class
            # and just tack on the current label duration
            return FixedLDPCASLProtocol.trial_param_values(self, params, idx)
        else:
            # We are varying the labelling duration.
            return np.round(np.arange(self.ld_lims.lb, self.ld_lims.ub+0.001, self.ld_lims.step), 5)

    def param_bounds(self):
        return [
            (self.pld_lims.lb, self.pld_lims.ub)
            for idx in range(self.scan_params.npld)
        ] + [
            (self.ld_lims.lb, self.ld_lims.ub)
        ]

    def _timings(self, params):
        plds = params[..., :-1]
        lds = np.zeros(plds.shape, dtype=np.float)
        lds[:] = params[..., -1][..., np.newaxis]
        return lds, plds

class MultiPLDPcaslMultiLD(MultiPLDPcaslVarLD):
    """
    PCASL protocol with multiple PLDs and multiple variable LDs
    """

    def __str__(self):
        return "Multi-PLD PCASL protocol with variable label durations (one per PLD)"

    def name_params(self, params):
        return {
            "plds" : params[:self.scan_params.npld],
            "lds" : params[self.scan_params.npld:],
        }

    def initial_params(self):
        plds = FixedLDPCASLProtocol.initial_params(self)

        # Add variable label durations to parameters
        lds = [self.scan_params.ld] * self.scan_params.npld
        return np.array(list(plds) + lds)

    def param_bounds(self):
        return [
            (self.pld_lims.lb, self.pld_lims.ub)
            for idx in range(self.scan_params.npld)
        ] + [
            (self.ld_lims.lb, self.ld_lims.ub)
            for idx in range(self.scan_params.npld)
        ]

    def _timings(self, params):
        return params[..., self.scan_params.npld:], params[..., :self.scan_params.npld]

class HadamardFixedLD(MultiPLDPcaslMultiLD):
    """
    Hadamard time-encoded protocol with single (variable) LD and single (variable) PLD

    The LD in this case is the sub-boli labelling duration
    """
    def __init__(self, *args, **kwargs):
        Protocol.__init__(self, *args, **kwargs)
        if self.scan_params.npld != 1:
            raise ValueError("Hadamard protocol must have a single PLD")
        self.had_size = kwargs.get("had_size", 8)

        # Slice time offset [NSlices]
        self.slicedt = np.arange(self.scan_params.nslices, dtype=np.float) * self.scan_params.slicedt

        # We can only calculate the cost function for ATT > shortest PLD, otherwise we don't
        # see the inflow and the determinant blows up. So we modify the ATT distribution
        # weight to be zero at ATTs which are not relevant to a slice.
        # Note that we do not renormalize the ATT distribution weights because we want slices
        # to contribute more when they have more relevant ATTs
        min_pld_possible = self.pld_lims.lb + self.slicedt
        relevant_atts = self.att_dist.atts[np.newaxis, ...] > min_pld_possible[..., np.newaxis]
        self.att_weights = self.att_dist.weight * relevant_atts

    def __str__(self):
        return "Hadamard time-encoded protocol with constant sub-boli label durations and PLD"

    def name_params(self, params):
        return {
            "plds" : params[0],
            "lds" : params[1],
        }

    def trial_param_values(self, params, idx):
        if idx == 0:
            lims = self.pld_lims
        else:
            lims = self.ld_lims

        return np.round(np.arange(lims.lb, lims.ub+0.001, lims.step), 5)

    def param_bounds(self):
        return [
            (self.pld_lims.lb, self.pld_lims.ub),
            (self.ld_lims.lb, self.ld_lims.ub)
        ]

    def repeats_total_tr(self, params):
        pld = params[..., 0]
        ld = params[..., 1]

        total_tr = self.had_size * ((self.had_size-1) * ld + pld + self.scan_params.readout)
        return np.floor(self.scan_params.duration/total_tr), np.round(total_tr, 5)

    def cost(self, params):
        # For comparison with other protocols need to scale
        # cost since each repetition gives 8 volumes of information
        # rather than 2 for a label/control acquisition
        return FixedLDPCASLProtocol.cost(self, params) / (self.had_size/2)

    def _timings(self, params):
        pld = params[..., 0][..., np.newaxis]
        ld = params[..., 1][..., np.newaxis]

        # Calculate the effective PLDs and LDs of the data after they have undergone
        # the addition/subtraction process to reduce the TE images into single-delay
        # images. In this case each of the effective LDs will be constant but the
        # effective PLDs will vary (e.g. a long PLD resulting from the isolation of
        # the first sub-bolus)
        effective_lds = np.repeat(ld, self.had_size-1, axis=-1)
        effective_plds = np.repeat(pld, self.had_size-1, axis=-1)
        effective_plds += np.arange(0, self.had_size-1, dtype=np.float32) * ld

        return effective_lds, effective_plds

class HadamardT1Decay(HadamardFixedLD):
    """
    Hadamard time-encoded protocol with variable sub-boli durations
    to account for T1 decay

    The LD in this case is the *first* sub-boli labelling duration
    """

    def _timings(self, params):
        pld = params[..., 0]
        first_ld = params[..., 1]

        eff_lds = self._effective_lds(first_ld)
        eff_plds = self._effective_plds(eff_lds, pld)

        return eff_lds, eff_plds

    def _effective_plds(self, lds, pld=0):
        """
        Generate the effective PLDs for each sub-bolus

        This is based on the fact that each sub-bolus has a delay
        from all the sub-boli that come after it plus the global PLD

        :param lds: Sub-boli LDs [NTrials, NBlocks] or [NBlocks]
        :param pld: Global PLD [NTrials] or scalar
        """
        eff_plds = np.zeros(lds.shape, dtype=np.float32)
        eff_plds += pld[..., np.newaxis]
        eff_plds[..., :-1] += np.cumsum(lds[..., :0:-1], -1)[..., ::-1]
        return eff_plds

    def _effective_lds(self, first_ld, step_size=0.025):
        """
        Generate optimal sub-bolus durations given the length of the first
        sub-bolus.

        We use a different method to that used in Joe's Matlab code - instead
        of optimizing the sub-boli durations to equalize the net Mz we calculate them
        by construction.

        :param first_ld: First sub-bolus duration [NTrials] or scalar
        :param step_size: Step size for durations - output is rounded to this precision
        """
        if np.array(first_ld).ndim == 0:
            first_ld = np.array(first_ld)[..., np.newaxis]

        # Derivation of the following is left as an exercise to the reader...
        t1b = self.kinetic_model._phys_params.t1b
        t = np.zeros(first_ld.shape)
        tau = np.zeros((first_ld.shape[0], self.had_size-1), dtype=np.float32)
        tau[:, 0] = first_ld
        v0 = np.exp(-t/t1b) * (1-np.exp(-first_ld/t1b))
        for idx in range(1, self.had_size-1):
            factor = np.exp(-t/t1b)
            x = v0 / factor + 1
            tauv = t1b*np.log(x)
            t -= tauv
            tau[:, idx] = tauv

        tau = np.round(tau / step_size) * step_size
        return np.squeeze(tau)
