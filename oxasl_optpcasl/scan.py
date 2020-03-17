import numpy as np

class Protocol(object):
    """
    A scan protocol that can be optimized

    A protocol has a set of variable parameters, can generate 'trial' values for 
    each parameter, and can calculate the cost of a set of parameters
    """

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
        Get a set of trial values for a single parameter

        :param params: Current initial param values [NParams]
        :param idx: Index of parameter to vary

        :return: Trial parameter values [Trials, NParams]
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

class MultiPLDPcasl(Protocol):
    """
    A multi-PLD PCASL protocol with fixed labelling duration
    """

    def __init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims):
        self.kinetic_model = kinetic_model
        self.cost_model = cost
        self.scan_params = scan_params
        self.att_dist = att_dist
        self.pld_lims = pld_lims

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
        return "Multi-PLD PCASL protocol with fixed label duration"

    def name_params(self, params):
        return {
            "plds" : params
        }

    def initial_params(self):
        # Parameters in this case are a set of PLDs
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
    
    def cost(self, params):
        #print("cparams", params.shape)
        # Hessian matrix for sensitivity of kinetic model [NTrials, NSlices, NATTs, 2, 2]
        hessian = self._hessian(params)
        #print("hessian", hessian.shape)
        #print(hessian)

        # Cost at each time point and for each ATT in distribution [NTrials, NSlices, NATTs]
        cost = self.cost_model.cost(hessian)
        #print("cost", cost.shape)
        #print(cost)

        # Weight each component of the cost according to ATT distribution weight
        # This also takes into account which ATTs are relevant to which slices
        att_weights = self.att_weights
        #print("att_weights")
        #print(att_weights)
        if params.ndim == 2:
            att_weights = np.tile(self.att_weights[np.newaxis, ...], (params.shape[0], 1, 1))
        #print(att_weights)
        cost *= att_weights
        cost[att_weights == 0] = 0      # To correct for 0*nan
        cost[np.isnan(cost)] = np.inf   # To correct for 0*inf
        #print("wcost", cost.shape)

        # Take mean of cost across slices and ATTs
        return np.mean(cost, axis=(-1, -2))
        
    def trial_params(self, params, idx):
        trial_values = self.trial_param_values(params, idx)
        
        # Trial PLDs list [Trials, NParams]
        # The PLD that we are testing in this loop gets a different trial value in each
        # column, the other PLDs are constant in each column
        # Adjust the PLDs to the slice delay time and calculate the TI (including labelling duration)
        trial_plds = np.tile(params[np.newaxis, :], (len(trial_values), 1))
        trial_plds[:, idx] = trial_values
        return trial_plds

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
        Get the labelling duration and PLDs for a set of trial params
        
        :param params: Parameter values [NTrials, NParams] or [NParams]
        :return Tuple of labelling duration and PLDs, each [NTrials, NPLDs] or [NPLDs]
        """
        return np.full(params.shape, self.scan_params.ld), params

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
        num_repeats, _ = self.repeats_total_tr(params)
        num_repeats = num_repeats/(self.scan_params.noise**2)#
        
        # Calculate derivatives of ASL kinetic model [NTrials, NPLDs, NSlices, NATTs]
        att = self.att_dist.atts
        df, datt = self.kinetic_model.sensitivity(lds.flatten(), times.flatten(), att)
        df = df.reshape(list(times.shape) + [len(att)])
        datt = datt.reshape(list(times.shape) + [len(att)])
        #print("df, datt\n", df, datt)
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

class MultiPLDPcaslVarLD(MultiPLDPcasl):
    """
    PCASL protocol with multiple PLDs and single variable LD
    """

    def __init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims, ld_lims):
        MultiPLDPcasl.__init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims)
        self.ld_lims = ld_lims

    def __str__(self):
        return "Multi-PLD PCASL protocol with single variable label duration"

    def name_params(self, params):
        return {
            "plds" : params[:self.scan_params.npld],
            "lds" : params[self.scan_params.npld],
        }

    def initial_params(self):
        plds = MultiPLDPcasl.initial_params(self)

        # Add variable label duration to parameters
        return np.array(list(plds) + [self.scan_params.ld])

    def trial_params(self, params, idx):
        if idx < self.scan_params.npld:
            # We are varying a PLD. Get the trial values from the base class
            # and just tack on the current label duration
            trial_plds = MultiPLDPcasl.trial_params(self, params[:-1], idx)
            trial_params = np.zeros((trial_plds.shape[0], self.scan_params.npld+1), dtype=np.float)
            trial_params[:, :self.scan_params.npld] = trial_plds
            trial_params[:, self.scan_params.npld] = params[self.scan_params.npld]
        else:
            # We are varying the labelling duration.
            trial_lds = np.round(np.arange(self.ld_lims.lb, self.ld_lims.ub+0.001, self.ld_lims.step), 5)
            trial_params = np.zeros((len(trial_lds), self.scan_params.npld+1), dtype=np.float)
            trial_params[:, :] = params
            trial_params[:, idx] = trial_lds

        return trial_params

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

class MultiPLDPcaslMultiLD(MultiPLDPcasl):
    """
    PCASL protocol with multiple PLDs and multiple variable LDs
    """

    def __init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims, ld_lims):
        MultiPLDPcasl.__init__(self, kinetic_model, cost, scan_params, att_dist, pld_lims)
        self.ld_lims = ld_lims

    def __str__(self):
        return "Multi-PLD PCASL protocol with variable label durations (one per PLD)"

    def name_params(self, params):
        return {
            "plds" : params[:self.scan_params.npld],
            "lds" : params[self.scan_params.npld:],
        }

    def initial_params(self):
        plds = MultiPLDPcasl.initial_params(self)

        # Add variable label durations to parameters
        lds = [self.scan_params.ld] * self.scan_params.npld

        return np.array(list(plds) + lds)

    def trial_params(self, params, idx):
        if idx < self.scan_params.npld:
            # We are varying a PLD. Get the trial values from the base class
            # and just tack on the current label durations
            trial_plds = MultiPLDPcasl.trial_params(self, params, idx)
            trial_params = np.zeros((trial_plds.shape[0], self.scan_params.npld*2), dtype=np.float)
            trial_params[:, :trial_plds.shape[1]] = trial_plds
            trial_params[:, trial_plds.shape[1]:] = params[trial_plds.shape[1]:]
        else:
            # We are varying a labelling duration.
            trial_lds = np.round(np.arange(self.ld_lims.lb, self.ld_lims.ub+0.001, self.ld_lims.step), 5)
            trial_params = np.zeros((len(trial_lds), self.scan_params.npld*2), dtype=np.float)
            trial_params[:, :] = params
            trial_params[:, idx] = trial_lds

        return trial_params

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
