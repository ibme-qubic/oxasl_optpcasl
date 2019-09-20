"""
oxasl.gui.options.py

TabPage widgets which allow the user to change options

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from ..structures import VAR_MULTI_PCASL, VAR_TE_PCASL, LOOK_LOCKER, VAR_TE_PCASL_NPLD, ASLParams, Scan, BATDist, Limits
from ..optimize import LOptimal, DOptimal
from .widgets import TabPage

class ScanOptions(TabPage):
    """
    Tab page for specifying the ASL scan to optimize
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")

        self.section("ASL scan to optimize for")
        self._asltype = self.choice("ASL scan type", choices=[VAR_MULTI_PCASL, VAR_TE_PCASL, LOOK_LOCKER, VAR_TE_PCASL_NPLD])
        self._fvalue = self.number("Estimated perfusion (ml/100g/min)", minval=0, maxval=100.0, initial=50.0)
        self._duration = self.number("Approximate scan duration (s)", minval=0, maxval=1000, initial=300)
        self._readout = self.number("Readout time (s)", minval=0, maxval=2.0, initial=0.5)
        self._nplds = self.integer("Number of PLDs", minval=1, maxval=20, initial=6)

        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)

    def aslparams(self):
        return ASLParams(self._asltype.GetString(self._asltype.GetSelection()), 
                         self._fvalue.GetValue()/6000.0)

    def scan(self):
        return Scan(duration=self._duration.GetValue(), 
                    npld=self._nplds.GetValue(), 
                    readout=self._readout.GetValue())

class OptimizerOptions(TabPage):
    """
    Tab page for specifying options for the optimization
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")

        self.section("Optimization type")
        self._opttype = self.choice("Method", choices=["Optimize CBF and ATT", "Optimize CBF only", "Optimize ATT only"])

        self.section("ATT prior distribution")
        self._att_start = self.number("Starting value (s)", minval=0, maxval=1.0, initial=0.2)
        self._att_end = self.number("Starting value (s)", minval=0, maxval=5.0, initial=2.1)
        self._att_step = self.number("Step (s)", minval=0, maxval=0.01, initial=0.001, digits=4)
        self._att_taper = self.number("Taper value (s)", minval=0, maxval=1.0, initial=0.3)
        
        self.section("PLD search limits")
        self._pld_min = self.number("Min PLD (s)", minval=0, maxval=1.0, initial=0.1)
        self._pld_max = self.number("Max PLD (s)", minval=1.0, maxval=5.0, initial=3.0)
        self._pld_step = self.number("Search step (s)", minval=0, maxval=0.1, initial=0.025, digits=4)
        
        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)

    def attdist(self):
        return BATDist(self._att_start.GetValue(), self._att_end.GetValue(), 
                       self._att_step.GetValue(), self._att_taper.GetValue())

    def pldlimits(self):
        return Limits(self._pld_min.GetValue(), self._pld_max.GetValue(), self._pld_step.GetValue())

    def optimizer(self):
        if self._opttype.GetSelection() == 1:
            return LOptimal([[1, 0],  [0, 0]])
        if self._opttype.GetSelection() == 2:
            return LOptimal([[0, 0],  [0, 1]])
        else:
            return DOptimal()
