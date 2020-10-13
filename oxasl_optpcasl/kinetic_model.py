"""
OXASL_OPTPCASL - Kinetic models for PCASL optimizer
"""
import numpy as np
import scipy.special

class KineticModel(object):
    """
    Base class for a kinetic model

    The kinetic model must be able to calculate the sensitivity of CBF and ATT to specified
    labelling durations and time points at a range of ATTs
    """

    def __init__(self, phys_params):
        self.phys_params = phys_params

    def signal(self, ld, times, att, m0=1.0, pldmax=5.0):
        """
        This function calculates the output of the kinetic model i.e. the
        predicted signal
        
        :param ld: Labelling duration(s). This must have a shape compatible with
                   ``times``, i.e one LD per time point or a single scalar LD.
        :param times: Time points [NT]
        :param att: ATT distribution [ATTs]

        :return: Signal [NT, ATTs]
        """
        raise NotImplementedError()
        
    def sensitivity(self, ld, times, att):
        """
        This function calculates the sensitivities (partial derivative) of the kinetic
        model

        :param ld: Labelling duration(s). This must have a shape compatible with
                   ``times``, i.e one LD per time point or a single scalar LD.
        :param times: Time points [NT]
        :param att: ATT distribution [ATTs]

        :return: Tuple of df [NT, ATTs], datt [NT, ATTs]
        """
        raise NotImplementedError()

class BuxtonPcasl(KineticModel):
    """
    Basic one-compartment Buxton model
    """
    def __init__(self, phys_params):
        """
        This class implements the Buxton CASL model (Buxton et al. MRM 1998) 
        given in Woods et al. MRM 2019.

        Note that we assume a fixed flow value in t1_prime in order
        to simplify the equations. It is shown in Chappell et al (FIXME ref) that
        this causes negligible error in the kinetic model and avoids a circular
        dependency of the model on its own output.
        """
        KineticModel.__init__(self, phys_params)
        self.t1_prime = 1.0/((1.0/self.phys_params.t1t) + (self.phys_params.f/self.phys_params.lam))

    def _preproc(self, ld, times, att):
        # Set up output arrays - take advantage of Numpy broadcasting to combine times with ATTs
        # All these arrays have shape [NT, ATTs]
        att = np.atleast_1d(att)
        ld = np.repeat(np.atleast_1d(ld)[..., np.newaxis], len(att), -1)
        times = np.repeat(times[..., np.newaxis], len(att), -1)

        # Weights to show which times were during or after bolus arrival. If using ERF pseudo-dispersion
        # these are floating point weights. If not, they are effectively boolean weights (1 or 0)
        USE_ERF = False
        if not USE_ERF:
            # When the TI is effectively equal to LD + ATT we get instabilities due to
            # rounding which can cause Matlab and Python to categorize time points
            # differently. This array is designed to define TIs which are 'effectively'
            # equal so the rounding differences do not matter so much
            times_equal = np.isclose(times, ld + att)
            weight_during = np.logical_and(times > att, np.logical_or(times < (ld + att), times_equal)).astype(np.int)
            weight_after = np.logical_and(times > ld + att, np.logical_not(times_equal)).astype(np.int)
        else:
            # The alternative is to smoothly blend the before/after bolus arrival solutions
            # using an error function with a short time scale (0.01s in this case)
            weight_during = ((1+scipy.special.erf((times-att)/0.01))/2) * (scipy.special.erfc((times-ld-att)/0.01)/2)
            weight_after = (1+scipy.special.erf((times-att-ld)/0.01))/2

        return ld, times, weight_during, weight_after

    def signal(self, ld, times, att):
        ld, times, weight_during, weight_after = self._preproc(ld, times, att)

        M = np.zeros(times.shape, dtype=np.float32)

        M += weight_during * np.exp(-att/self.phys_params.t1b) * (1-np.exp(-(times-att)/self.t1_prime))
        M += weight_after * np.exp(-att/self.phys_params.t1b) * np.exp(-(times-ld-att)/self.t1_prime) * (1-np.exp(-ld/self.t1_prime))

        M *= 2 * self.phys_params.f * self.phys_params.m0b * self.phys_params.alpha * self.t1_prime 
        return M

    def sensitivity(self, ld, times, att):
        ld, times, weight_during, weight_after = self._preproc(ld, times, att)
        M = 2 * self.phys_params.m0b * self.phys_params.alpha * self.t1_prime * np.exp(-att / self.phys_params.t1b) # [ATTs]

        # For t between deltaT and label duration plus deltaT
        df_during = M * (1 - np.exp((att - times) / self.t1_prime))
        datt_during = M * self.phys_params.f * ((-1.0/self.phys_params.t1b) - np.exp((att - times) / self.t1_prime) * ((1.0/self.t1_prime) - (1.0/self.phys_params.t1b)))

        # for t greater than ld plus deltaT
        df_after = M * np.exp((ld + att - times) / self.t1_prime) * (1 - np.exp(-ld/self.t1_prime))
        datt_after = M * self.phys_params.f * (1 - np.exp(-ld/self.t1_prime)) * np.exp((ld + att - times)/self.t1_prime) * (1.0/self.t1_prime - 1.0/self.phys_params.t1b)

        df = df_during * weight_during + df_after * weight_after
        datt = datt_during * weight_during + datt_after * weight_after

        return df, datt
