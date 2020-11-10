"""
OXASL_OPTPCASL - ASL Scan protocols

Copyright 2019 University of Oxford
"""
import numpy as np
import scipy.linalg

def get_protocol(options):
    """
    Get the protocol class from command line options
    """
    if options.protocol == "pcasl":
        if options.optimize_ld:
            if options.multi_ld:
                return MultiPLDPcaslMultiLD
            else:
                return MultiPLDPcaslVarLD
        else:
            return FixedLDPcaslProtocol
    elif options.protocol == "hadamard":
        return HadamardSingleLd
    elif options.protocol == "hadamardt1":
        return HadamardT1Decay
    elif options.protocol == "hadamardvar":
        return HadamardMultiLd
    elif options.protocol == "hadamardfl":
        return HadamardFreeLunch
    else:
        raise ValueError("Unrecognized protocol: %s" % options.protocol)

class Protocol(object):
    """
    A scan protocol that can be optimized

    A protocol has a set of variable parameters, can generate 'trial' values for
    each parameter, and can calculate the cost of a set of parameters
    """

    def __init__(self, kinetic_model, scan_params, att_dist, pld_lims, ld_lims=None):
        self.kinetic_model = kinetic_model
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

    def cost(self, params, cost_model):
        """
        Get the cost for a set of parameters

        :param params: Parameter set [NTrials, NParams] or [NParams]
        :param cost_model: Cost model
        :return Cost [NTrials] or scalar
        """
        raise NotImplementedError()

    def hessian(self, params):
        """
        Calculate Hessian matrix for ATT and CBF

        :param params: Parameters [NTrials, NParams] or [NParams]

        :return: Hessian matrices of second derivatives wrt cbf
                 and att [NTrials, NSlices, NATTs, 2, 2] or [NSlices, NATTs, 2, 2]
        """
        raise NotImplementedError()

    def cov(self, params):
        """
        Calculate covariance matrix for ATT and CBF

        This is the inverse of the Hessian with unit conversion and inf/NaN 
        protection

        :param params: Parameters [NTrials, NParams] or [NParams]

        :return: Covariance matrix for CBF and ATT in standard units
                 and att [NTrials, NSlices, NATTs, 2, 2] or [NSlices, NATTs, 2, 2]
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

    def protocol_summary(self, params):
        """
        Get a summary of the protocol for a set of parameters

        :param params: Parameter values [NParams]
        :return: Sequence of tuples. Each tuple describes an acquisition and consists of 
                 the following elements: name, sequence of labelling durations (typically
                 only one, but may be multiple in case of time encoded protocols), matching 
                 sequence of 1 or 0 indicating tag/control, PLD and
                 readout time (all in seconds).
        """
        raise NotImplementedError()

class PcaslProtocol(Protocol):
    """
    A PCASL protocol which may have multiple PLDs and independent
    labelling durations associated with them.

    The first ``scan_params.npld`` parameters are the PLDs and any remaining
    parameters are labelling durations (may be none, one or multiple)
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

        self.nld = 0

    def name_params(self, params):
        ret = {"plds" : params[:self.scan_params.npld]}
        if len(params) > self.scan_params.npld:
            ret["lds"] = params[self.scan_params.npld:]
        return ret

    def initial_params(self):
        param_values = self._initial_plds()
        
        if self.nld > 0:
            # Add variable label durations to parameters
            param_values = list(param_values) + list(self._initial_lds())
        return np.array(param_values)

    def param_bounds(self):
        return [
            (self.pld_lims.lb, self.pld_lims.ub)
            for idx in range(self.scan_params.npld)
        ] + [
            (self.ld_lims.lb, self.ld_lims.ub)
            for idx in range(self.nld)
        ]

    def trial_param_values(self, params, idx):
        if idx < self.scan_params.npld:
            # We are varying a PLD. Generate values between the previous and following PLDs
            # to keep them in increasing order.For the first and last PLDs we use the upper/lower 
            # bounds instead of the previous/next pld.
            start, stop = self.pld_lims.lb, self.pld_lims.ub
            if idx > 0:
                start = params[idx-1]
            if idx < self.scan_params.npld-1:
                stop = params[idx+1]

            return np.round(np.arange(start, stop+0.001, self.pld_lims.step), 5)
        else:
            # We are varying a labelling duration. Generate values between upper and lower bounds
            # since LDs shouldn't be kept in any particular order
            start, stop = self.ld_lims.ub, self.ld_lims.lb
            ret = np.round(np.arange(start, stop+0.001, -self.ld_lims.step), 5)
            return ret

    def cov(self, params):
        hessian = self.hessian(params)
        det = np.linalg.det(hessian)
        #print("det", det)
        cov = np.zeros(hessian.shape)
        cov[det != 0] = np.linalg.inv(hessian[det != 0])
        cov[det == 0] = np.inf
        #print("cov\n", cov)
        # Correct for inf*0 errors in A*inverse
        #cov[np.isnan(cov)] = np.inf

        # Change into (ml/100g/min)
        cov[..., 0, 0] = cov[..., 0, 0] * 6000 * 6000
        cov[..., 0, 1] = cov[..., 0, 1] * 6000
        cov[..., 1, 0] = cov[..., 1, 0] * 6000

        return cov

    def cost(self, params, cost_model):
        # Covariance matrix for sensitivity of kinetic model [NTrials, NSlices, NATTs, 2, 2]
        cov = self.cov(params)

        # Cost at each time point and for each ATT in distribution [NTrials, NSlices, NATTs]
        cost = cost_model.cost(cov)
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
        lds, plds = self.timings(params)
        tr = lds + plds + self.scan_params.readout
        total_tr = 2*np.sum(tr, axis=-1)

        # I round the tr since there are occasionally computational rounding
        # errors. I have used 5 decimal places to allow dense att sampling, but I am
        # very unlikely to go finer than 0.001 density.
        return np.floor(self.scan_params.duration/total_tr), np.round(total_tr, 5)

    def _initial_plds(self):
        if self.scan_params.plds is not None:
            return self.scan_params.plds
        else:
            factor = 1.0/self.pld_lims.step
            max_pld = self.pld_lims.ub
            while 1:

                # Initial sequence of PLDs spaced evenly between upper and lower bound
                plds = np.linspace(self.pld_lims.lb, max_pld, self.scan_params.npld)
                plds = np.round(plds*factor) / factor

                # See how long the TR is - if it is larger than the maximum, reduce the maximum PLD FIXME lds
                total_tr = 2*np.sum(max(self.scan_params.ld) + plds + self.scan_params.readout)
                if total_tr <= self.scan_params.duration:
                    break
                max_pld -= 0.1

            return plds

    def _initial_lds(self):
        nld = len(self.scan_params.ld)
        if nld == self.nld:
            return self.scan_params.ld
        elif nld == 1 and self.nld > 1:
            # Initial LD is 1.8s unless upper bound is below this
            return [min(self.ld_lims.ub, 1.8)] * self.nld
        else:
            raise ValueError("Invalid number of initial labelling durations passed (%i provided, expected %i" % (nld, self.nld))

    def timings(self, params):
        """
        Get the effective labelling duration and PLDs for a set of trial params

        :param params: Parameter values [NTrials, NParams] or [NParams]
        :return Tuple of labelling duration and PLDs, each [NTrials, NPLDs] or [NPLDs]
        """
        raise NotImplementedError()

    def hessian(self, params):
        # Time points to evaluate the sensitivity of the kinetic model at
        # [NTrials, NPLDs, NSlices]
        lds, plds = self.timings(params)
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

    def all_lds(self, lds, step_size=0.025):
        """
        Return the full set of LDs given the LD parameter(s)
        """
        return lds

class FixedLDPcaslProtocol(PcaslProtocol):
    """
    PCASL protocol with a single fixed labelling duration
    """
    
    def __init__(self, *args, **kwargs):
        PcaslProtocol.__init__(self, *args, **kwargs)
        self.nld = 0

    def __str__(self):
        return "Multi-PLD PCASL protocol with fixed label duration"

    def timings(self, params):
        return np.full(params.shape, self.scan_params.ld[0]), params

    def protocol_summary(self, params):
        ret = []
        for pld in params:
            ret.append(("", self.scan_params.ld, [1], pld, self.scan_params.readout))
            ret.append(("", self.scan_params.ld, [0], pld, self.scan_params.readout))
        return ret

class MultiPLDPcaslVarLD(PcaslProtocol):
    """
    Pcasl protocol with multiple PLDs and single variable LD
    """

    def __init__(self, *args, **kwargs):
        PcaslProtocol.__init__(self, *args, **kwargs)
        self.nld = 1

    def __str__(self):
        return "Multi-PLD PCASL protocol with single variable label duration"

    def timings(self, params):
        plds = params[..., :-1]
        lds = np.zeros(plds.shape, dtype=np.float)
        lds[:] = params[..., -1][..., np.newaxis]
        return lds, plds

    def protocol_summary(self, params):
        ret = []
        ld = params[self.scan_params.npld]
        for pld in params[:self.scan_params.npld]:
            ret.append(("", [ld], [1], pld, self.scan_params.readout))
            ret.append(("", [ld], [0], pld, self.scan_params.readout))
        return ret

class MultiPLDPcaslMultiLD(PcaslProtocol):
    """
    PCASL protocol with multiple PLDs and multiple variable LDs
    """

    def __init__(self, *args, **kwargs):
        PcaslProtocol.__init__(self, *args, **kwargs)
        self.nld = self.scan_params.npld

    def __str__(self):
        return "Multi-PLD PCASL protocol with variable label durations (one per PLD)"

    def timings(self, params):
        return params[..., self.scan_params.npld:], params[..., :self.scan_params.npld]

    def protocol_summary(self, params):
        ret = []
        plds = params[:self.scan_params.npld]
        lds = params[self.scan_params.npld:]
        for pld, ld in zip(plds, lds):
            ret.append(("", [ld], [1], pld, self.scan_params.readout))
            ret.append(("", [ld], [0], pld, self.scan_params.readout))
        return ret

class Hadamard(PcaslProtocol):
    """
    Hadamard time-encoded protocol
    """
    def __init__(self, *args, **kwargs):
        PcaslProtocol.__init__(self, *args, **kwargs)
        self.had_size = self.scan_params.had_size
        self.nld = 1

    def repeats_total_tr(self, params):
        plds = params[..., :self.scan_params.npld]
        sub_boli = self.all_lds(params[..., self.scan_params.npld:], self.ld_lims.step)
        total_tr = 0

        for pld_idx in range(self.scan_params.npld):
            pld = plds[..., pld_idx]
            total_tr += self.had_size * (np.sum(sub_boli, axis=-1) + pld + self.scan_params.readout)
        return np.floor(self.scan_params.duration/total_tr), np.round(total_tr, 5)

    def cost(self, params, cost_model):
        # For comparison with other protocols need to scale
        # cost since each repetition gives self.had_size volumes of information
        # rather than 2 for a label/control acquisition
        return PcaslProtocol.cost(self, params, cost_model) / (self.had_size/2)

    def protocol_summary(self, params):
        plds = params[..., :self.scan_params.npld]
        lds = self.all_lds(params[self.scan_params.npld:], self.ld_lims.step)
        had = scipy.linalg.hadamard(self.had_size)
        ret = []
        for pld in plds:
            for row in had:
                ret.append(("", lds, row[1:], pld, self.scan_params.readout))
        return ret

    def timings(self, params):
        plds = params[..., :self.scan_params.npld]
        eff_lds = self.all_lds(params[..., self.scan_params.npld:], self.ld_lims.step)
        eff_lds_full = np.tile(eff_lds, self.scan_params.npld)
        eff_plds_full = np.zeros(eff_lds_full.shape)
        for pld in range(self.scan_params.npld):
            eff_plds_full[..., pld*(self.had_size-1):(pld+1)*(self.had_size-1)] = self._effective_plds(eff_lds, plds[..., pld])

        return eff_lds_full, eff_plds_full

    def _effective_plds(self, lds, pld=0):
        """
        Generate the effective PLDs for each sub-bolus

        This is based on the fact that each sub-bolus has a delay
        from all the sub-boli that come after it plus the global PLD

        :param lds: Sub-boli LDs [NTrials, NBlocks] or [NBlocks]
        :param pld: Global PLD [NTrials] or scalar
        :return effective PLD for each sub-bolus [NTrials, NBlocks]
        """
        pld = np.array(pld)[..., np.newaxis]
        eff_plds = np.repeat(pld, self.had_size-1, axis=-1)
        eff_plds[..., :-1] += np.cumsum(lds[..., :0:-1], -1)[..., ::-1]
        return eff_plds

class HadamardSingleLd(Hadamard):
    """
    Hadamard time-encoded protocol with single (variable) LD and single (variable) PLD

    The LD in this case is the sub-boli labelling duration
    """

    def __str__(self):
        return "Hadamard time-encoded protocol with equal sub-boli label durations and single PLD"

    def all_lds(self, ld, step_size=0.025):
        return np.repeat(ld, self.had_size-1, axis=-1)

class HadamardT1Decay(Hadamard):
    """
    Hadamard time-encoded protocol with variable sub-boli durations
    to account for T1 decay

    The LD in this case is the *first* sub-boli labelling duration
    """

    def __str__(self):
        return "Hadamard time-encoded protocol with sub-boli label durations chosen to equalize overall effect of each sub-bolus given T1 decay"

    def all_lds(self, lds, step_size=0.025):
        """
        Generate optimal sub-bolus durations given the length of the first
        sub-bolus.

        We use a different method to that used in Joe's Matlab code - instead
        of optimizing the sub-boli durations to equalize the net Mz we calculate them
        by construction.

        :param lds: Sub-bolus durations [NTrials, NLDs] or [NLDs] - note NLDs should be 1
        :param step_size: Step size for durations - output is rounded to this precision
        :return Sub-bolus durations [NTrials, NBlocks]
        """
        first_ld = lds[..., 0]
        if first_ld.ndim == 0:
            first_ld = first_ld[..., np.newaxis]

        # Derivation of the following is left as an exercise to the reader...
        t1b = self.kinetic_model.phys_params.t1b
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
        if lds.ndim == 1:
            tau = np.squeeze(tau)
        return tau

class HadamardFreeLunch(HadamardT1Decay):
    """
    Hadamard time-encoded protocol that uses the 'spare' post labelling delay from
    a single-delay acquisition to cram in some more label/control TE sequences
    
    The PLD in this case is the notional single-delay PLD while the LD is the 
    *first* sub-boli labelling duration after the initial labelling. The remainder
    are chosen using the T1-decay approach.

    The initial labelling duration is fixed by the scan parameters  
    """
    def __str__(self):
        return "Hadamard free-lunch protocol with post-initial sub-boli label durations chosen to equalize overall effect of each sub-bolus given T1 decay"

    def _initial_lds(self):
        # For free-lunch we deal with special case where two initial LDs are given. The
        # first is the fixed LD and the second the variable first sub-bolus LD. If only
        # one is given we use it for both these cases
        if len(self.scan_params.ld) == 2:
            return [self.scan_params.ld[1]]
        else:
            return PcaslProtocol._initial_lds(self)

    def all_lds(self, lds, step_size=0.025):
        """
        Generate optimal sub-bolus durations given the length of the first
        sub-bolus.

        We use a different method to that used in Joe's Matlab code - instead
        of optimizing the sub-boli durations to equalize the net Mz we calculate them
        by construction.

        :param lds: Sub-bolus durations [NTrials, NLDs] or [NLDs] - note NLDs should be 1
        :param step_size: Step size for durations - output is rounded to this precision
        :return Sub-bolus durations [NTrials, NBlocks]
        """
        first_ld = lds[..., 0]
        if first_ld.ndim == 0:
            first_ld = first_ld[..., np.newaxis]

        # Derivation of the following is left as an exercise to the reader...
        t1b = self.kinetic_model.phys_params.t1b
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

        # Discard the last sub-bolus becuase we're going to substitute the fixed 
        # LD at the start
        tau = np.squeeze(np.round(tau / step_size) * step_size)[..., :-1]

        initial_ld = np.full(tau[..., 0].shape, self.scan_params.ld[0])[..., np.newaxis]
        return np.concatenate((initial_ld, tau), axis=-1)
        
class HadamardMultiLd(Hadamard):
    """
    Hadamard time-encoded protocol where all sub-bolus LDs are allowed to vary independently and there
    is also a single (variable) PLD

    There number of LDs is the Hadamard size - 1 (e.g. 7 for a 8x7 Hadamard matrix)
    """
    def __init__(self, *args, **kwargs):
        Hadamard.__init__(self, *args, **kwargs)
        self.nld = self.had_size - 1

    def __str__(self):
        return "Hadamard time-encoded protocol with multiple variable sub-boli label durations and single PLD"
