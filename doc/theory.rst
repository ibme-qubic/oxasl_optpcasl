Theory
======

This is a **brief** overvies of the theory behind multi-PLD PCASL optimization. 
For a more detailed account see Woods et al, 2018 [1]_, 2020 [2]_

``OXASL_OPTPCASL`` is based on minimising the Cramer-Rao lower bound on the variance 
of the parameters estimated from the data which for a multi-PLD PCASL experiment is
typically the CBF and the ATT.

In this case the lower bound is obtained from the inverse of the Fisher information 
matrix which takes the following form:

.. math::
    F(t; \theta)_{jk} = \frac{A}{\sigma^2} \sum_{i=1}^N \frac{\partial \Delta M(t_i; \theta)}{\partial \theta_j} \frac{\partial \Delta M(t_i; \theta)}{\partial \theta_k} 

Here :math:`N` is the number of acquisitions (PLDs), :math:`A` is the number of averages
(repeats) for each PLD, :math:`\sigma^2` is the normally distributed noise variance
and :math:`t` are the timings (i.e. PLDs).

:math:`\Delta M(t; \theta)` is a model for the ASL signal (magnitude difference between 
labelled and control images) based on model parameters :math:`\theta` which in this
case we will take as being CBF and ATT only.

``OPTPCASL`` supports two approaches to minimising the Cramer-Rao lower bound. In
each case the elements of the Fisher information matrix are evaluated by differentiating
the basic PCASL kinetic model:

.. math::
    \begin{array}{rlr}
    \Delta M(t) & = 0                                                                                                                                           & 0 < t < \Delta t \\
                & = 2M_{0B} f T_1^\prime \alpha \exp{(\frac{-\Delta t}{T_{1b}})} (1 - \exp{\frac{-(t-\Delta t)}{T_1^\prime}})                                   & \Delta t < t < \tau + \Delta t \\
                & = 2M_{0B} f T_1^\prime \alpha \exp{(\frac{-\Delta t}{T_{1b}})} \exp{\frac{-t-\tau-\Delta t}{T_1^\prime}} (1 - \exp{\frac{-\tau}{T_1^\prime}}) & \tau + \Delta t < t
    \end{array}

The apparent :math:`T_1` relaxation time :math:`T_1^\prime` is calculated assuming 
a fixed CBF of 50 ml/100g/min - this limitation in practice leads to an insignificant
error in the calculation of the kinetic model [1]_.

Optimising for both CBF and ATT variance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case we seek to minimise the cost function:

.. math::
    C = \frac{1}{det(F(t; \theta)})

i.e. minimising the overall magnitude of the covariance matrix

Optimising for either CBF or ATT variance only
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case we simply minimise a particular component of :math:`F^{-1}` which 
corresponds to the variance of the parameter of interest - e.g. if CBF is the 
first parameter we would minimise :math:`F^{-1}_{00}` to optimize for CBF and minimise
:math:`F^{-1}_{11}` to optimize for ATT.

Priors
~~~~~~

Since the elements of F depend on the parameters being estimated we adopt a 
pseudo-Bayesian approach where the cost function being minimised is averaged
across a prior distribution for the parameters. 

In the case of the ATT, a essentially uniform distribution is used with minimum
and maximum ATT values. In addition, a *tapering* parameter is supported which reduces
the probability density linearly within a fixed distance from the minimum and maximum
times (default: 0.3s).

A prior distribution for CBF is not required since (with the use of a fixed CBF in the
estimation of :math:`T_1^\prime`), the cost function is essentially independent of CBF.

Implementation of the cost minimisation algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Minimisiation is accomplished by an iterative search method looping over the PLDs in turn. For 
each PLD the optimal value is obtained for that PLD allowing it to vary between the
current values of the preceding and successive PLD. This optimization loop is repeated
until the set of PLDs is unchanged.

References
----------

.. [1] Woods JG, Chappell MA, Okell TW. A general framework for optimizing
       arterial spin labelling MRI experiments. Magn Reson Med. 2019;81:2474-2488.
       https://doi.org/10.1002/mrm.27580
       
.. [2] Woods JG, Chappell MA, Okell TW. Designing and comparing optimized 
       pseudo-continuous Arterial Spin Labeling protocols for measurement of 
       cerebral blood flow. NeuroImage 2020 223:1053-8119
       doi:10.1016/j.neuroimage.2020.117246
