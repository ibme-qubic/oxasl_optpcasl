Example use of the GUI
======================

The ``OXASL_OPTPCASL`` GUI is started using the command ``oxasl_optpcasl_gui``

Setting the scan parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: /images/gui_scan.png
    :alt: Setting scan parameters

On the first tab we set the parameters for the scan we want to optimize.
If we are using a 2D multi-slice readout, additional options for the 
time per slice and number of slices are enabled:

.. image:: /images/gui_scan_2d.png
    :alt: Setting scan parameters

Optimization configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

On the second tab we can control the optimization parameters. This includes
the prior distribution of ATT and the search limits for PLDs. However the
most likely setting to change here is whether we want to optimize for
CBF, ATT or both.

.. image:: /images/gui_opt.png
    :alt: Setting scan parameters

Running the optimization
~~~~~~~~~~~~~~~~~~~~~~~~

Click the ``Run Optimization`` button to start - normally it takes a few
seconds - possibly as much as a minute with multi-slice or with very 
fine PLD search steps.

Viewing output analysis
~~~~~~~~~~~~~~~~~~~~~~~

The plot to the right of the configuration can display the following
analysis (selected from the menu above the plot):

 - CBF error displayed as a plot of the estimated standard deviation on
   the CBF measurement at a range of (true) ATT values
 - ATT error displayed as a plot of standard deviation on ATT measurement
   at a range of (true) ATT values
 - A plot of the ASL kinetic curve at a range of ATT values showing 
   the positions of the PLDs selected

For example, these plots were obtained by optimizing
for ATT and CBF on a 3D acquisition:

.. image:: /images/gui_cbferr_both.png
    :alt: CBF error

.. image:: /images/gui_atterr_both.png
    :alt: ATT error

.. image:: /images/gui_curve_both.png
    :alt: Kinetic curve

The following equivalent plots were obtained by optimizing for CBF
alone:

.. image:: /images/gui_cbferr_cbf.png
    :alt: CBF error

.. image:: /images/gui_atterr_cbf.png
    :alt: ATT error

.. image:: /images/gui_curve_cbf.png
    :alt: Kinetic curve

Note that in this case the estimated error on CBF is reduced, while the
estimated error on ATT has increased. The optimizer has placed more
PLDs on or after the signal curve peak when estimating CBF rather than
before the peak (which is a better location for estimating ATT).
