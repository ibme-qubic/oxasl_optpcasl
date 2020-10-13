Command line tool user guide
============================

The command line tool is ``oxasl_optpcasl``. The default arguments optimize for a 6 PLD PCASL
acquisition::

    $ oxasl_optpcasl 
    OXASL - PCASL Optimizer 0.0.1.post12
    ====================================
    Optimizing PLDs for 300s 3D scan with readout time 0.500000s
    PLD search limits: PLDs between 0.10s and 3.00s in steps of 0.02500s
    Optimizing for 6 PLDs
    BAT distribution: 1901 values between 0.20s and 2.100000s (weight taper=0.30s)
    Optimization method: D-optimal

    Finished optimization after 30 iters - PLDs unchanged
    Optimal PLDs: [0.2, 0.7, 0.725, 1.55, 1.875, 2.075]
    num_av = 8
    Scan time = 296.400000
    DONE

Varying the number of PLDs
~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the command line option ``--scan-npld=<n>``::

    $ oxasl_optpcasl --scan-npld=8
    OXASL - PCASL Optimizer 0.0.1.post12
    ====================================
    Optimizing PLDs for 300s 3D scan with readout time 0.500000s
    PLD search limits: PLDs between 0.10s and 3.00s in steps of 0.02500s
    Optimizing for 8 PLDs
    BAT distribution: 1901 values between 0.20s and 2.100000s (weight taper=0.30s)
    Optimization method: D-optimal

    Finished optimization after 48 iters - PLDs unchanged
    Optimal PLDs: [0.2, 0.7, 0.7, 0.725, 1.5, 1.775, 1.95, 2.1]
    num_av = 6
    Scan time = 298.200000
    DONE

Optimizing for CBF or ATT individually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the command line options ``--optimize=CBF`` or ``--optimize=ATT``::

    $ oxasl_optpcasl --optimize=CBF
    OXASL - PCASL Optimizer 0.0.1.post12
    ====================================
    Optimizing PLDs for 300s 3D scan with readout time 0.500000s
    PLD search limits: PLDs between 0.10s and 3.00s in steps of 0.02500s
    Optimizing for 6 PLDs
    BAT distribution: 1901 values between 0.20s and 2.100000s (weight taper=0.30s)
    Optimization method: L-optimal

    Finished optimization after 24 iters - PLDs unchanged
    Optimal PLDs: [0.2, 1.175, 1.8, 2.025, 2.1, 2.1]
    num_av = 7
    Scan time = 291.200000
    DONE

Varying the ATT prior distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ATT prior distribution determines the range of arterial transit times
that are considered relevant when optimizing the performance of the protocol.
The default range is from 0.2s to 2.3s, with the weighting of the upper and 
lower 0.3s being linearly tapered to zero.

This is quite a wide range, and in some cases it may be more sensible to 
restric the ATT prior to a more limited range. For example the following 
will perform optimization based on ATTs between 1 and 2s, with the range
1s to 1.1s and the range 1.9s to 2s tapered in weight to zero::

    $ oxasl_optpcasl --att-start=1.0 --att-end=2.0 --att-step=0.001 --att-taper=0.1
    OXASL - PCASL Optimizer 0.0.1.post21
    ====================================
    Optimizing PLDs for 300s 3D scan with readout time 0.500000s
    PLD search limits: PLDs between 0.10s and 3.00s in steps of 0.02500s
    Optimizing for 6 PLDs
    BAT distribution: 1000 values between 1.00s and 2.000000s (weight taper=0.10s)
    Optimization method: D-optimal

    Finished optimization after 24 iters - PLDs unchanged
    Optimal PLDs: [0.6, 0.6, 0.6, 1.575, 1.8, 2.025]
    num_av = 8
    Scan time = 297.600000
    DONE

Varying the PLD search limits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example uses finer search parameters to determine optimal PLDs. Using
very fine search limits gives a more accurate result but the optimization
will take longer. In this example we also increase the search range of PLDs
up to 5s (although it turns out that such long PLDs are not optimal in this
case and the largest value returned is 2.08s)::

    $ oxasl_optpcasl --pld-min=0.1 --pld-max=5 --pld-step=0.01
    OXASL - PCASL Optimizer 0.0.1.post12
    ====================================
    Optimizing PLDs for 300s 3D scan with readout time 0.500000s
    PLD search limits: PLDs between 0.10s and 5.00s in steps of 0.01000s
    Optimizing for 6 PLDs
    BAT distribution: 1901 values between 0.20s and 2.100000s (weight taper=0.30s)
    Optimization method: D-optimal


    Finished optimization after 30 iters - PLDs unchanged
    Optimal PLDs: [0.2, 0.7, 0.71, 1.54, 1.87, 2.08]
    num_av = 8
    Scan time = 296.000000
    DONE

Varying other scan properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we optimize for a longer scan and modify the readout time (in seconds). 
Note that the number of averages (repeats) has increased::

    $ oxasl_optpcasl --scan-readout=0.75 --scan-duration=500
    OXASL - PCASL Optimizer 0.0.1.post12
    ====================================
    Optimizing PLDs for 500s 3D scan with readout time 0.750000s
    PLD search limits: PLDs between 0.10s and 3.00s in steps of 0.02500s
    Optimizing for 6 PLDs
    BAT distribution: 1901 values between 0.20s and 2.100000s (weight taper=0.30s)
    Optimization method: D-optimal

    Finished optimization after 36 iters - PLDs unchanged
    Optimal PLDs: [0.2, 0.7, 0.7, 0.925, 1.8, 2.0]
    num_av = 13
    Scan time = 499.850000
    DONE
