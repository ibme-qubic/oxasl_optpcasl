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
        self._phys_params = phys_params

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

    def sensitivity(self, ld, times, att):
        """
        This function calculates the sensitivties (partial derivative) of the Buxton CASL
        model (Buxton et al. MRM 1998) - given in Woods et al. MRM 2019.

        Note that we assume a fixed flow value in t1_prime in order
        to simplify the equations. It is shown in Chappell et al (FIXME ref) that
        this causes negligible error in the kinetic model and avoids a circular
        dependency of the model on its own output.
        """
        t1_prime = 1.0/((1.0/self._phys_params.t1t) + (self._phys_params.f/self._phys_params.lam))
        M = 2*self._phys_params.m0b * self._phys_params.alpha * t1_prime * np.exp(-att / self._phys_params.t1b) # [ATTs]

        # Set up output arrays - take advantage of Numpy broadcasting to combine times with ATTs
        # All these arrays have shape [NT, ATTs]
        times = np.repeat(times[..., np.newaxis], len(att), -1)
        ld = np.repeat(ld[..., np.newaxis], len(att), -1)

        # For t between deltaT and label duration plus deltaT
        df_during = M * (1 - np.exp((att - times) / t1_prime))
        datt_during = M * self._phys_params.f * ((-1.0/self._phys_params.t1b) - np.exp((att - times) / t1_prime) * ((1.0/t1_prime) - (1.0/self._phys_params.t1b)))

        # for t greater than ld plus deltaT
        df_after = M * np.exp((ld + att - times) / t1_prime) * (1 - np.exp(-ld/t1_prime))
        datt_after = M * self._phys_params.f * (1 - np.exp(-ld/t1_prime)) * np.exp((ld + att - times)/t1_prime) * (1.0/t1_prime - 1.0/self._phys_params.t1b)

        USE_ERF = False
        if not USE_ERF:
            # When the TI is effectively equal to LD + ATT we get instabilities due to
            # rounding which can cause Matlab and Python to categorize time points
            # differently. This array is designed to define TIs which are 'effectively'
            # equal so the rounding differences do not matter so much
            times_equal = np.isclose(times, ld + att)

            df = np.zeros(times.shape, dtype=np.float32)
            datt = np.zeros(times.shape, dtype=np.float32)
            t_during = np.logical_and(times > att, np.logical_or(times < (ld + att), times_equal))

            df[t_during] = df_during[t_during]
            datt[t_during] = datt_during[t_during]

            t_after = np.logical_and(times > ld + att, np.logical_not(times_equal))
            df[t_after] = df_after[t_after]
            datt[t_after] = datt_after[t_after]
        else:
            # The alternative is to smoothly blend the before/after bolus arrival solutions
            # using an error function with a short time scale (0.01s in this case)
            weight_during = ((1+scipy.special.erf((times-att)/0.01))/2) * (scipy.special.erfc((times-ld-att)/0.01)/2)
            weight_after = (1+scipy.special.erf((times-att-ld)/0.01))/2

            df = df_during * weight_during + df_after * weight_after
            datt = datt_during * weight_during + datt_after * weight_after

        return df, datt
