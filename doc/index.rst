
.. image:: images/oxasl.png
   :scale: 30 %
   :alt: OXASL logo
   :align: right

OXASL_OPTPCASL - Tool for optimizing PCASL acquisitions
=======================================================

OXASL_OPTPCASL is a package for generating optimal PLDs for 
PCASL experiments in order to maximise sensitivity to CBF, ATT
or both.

Installation
------------

If using FSL, it is simplest to install using ``fslpython`` (however
FSL is not required)::

    fslpython -m pip install oxasl_optpcasl --user

Otherwise can be installed into any Python 3 environment. Use of Conda 
or virtualenv is recommended::

    python -m pip install oxasl_optpcasl

You may need to add ``--user`` if using system python.

References
----------

 - Woods JG, Chappell MA, Okell TW. *A general framework for optimizing arterial
   spin labeling MRI experiments.* Magn Reson Med. 2019;81(4):2474-2488. 
   doi:10.1002/mrm.27580

 - Woods JG, Chappell MA, Okell TW. *Designing and comparing optimized 
   pseudo-continuous Arterial Spin Labeling protocols for measurement of 
   cerebral blood flow* NeuroImage 2020 223:1053-8119
   doi:10.1016/j.neuroimage.2020.117246

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   theory
   example_cli
   example_gui
