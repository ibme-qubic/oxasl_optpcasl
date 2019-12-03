"""
oxasl.gui.options.py

TabPage widgets which allow the user to change options

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from ..structures import VAR_MULTI_PCASL, VAR_TE_PCASL, LOOK_LOCKER, VAR_TE_PCASL_NPLD, ASLParams, ASLScan, ATTDist, Limits
from ..optimize import LOptimal, DOptimal
from .widgets import TabPage, NumberChooser

class ScanOptions(TabPage):
    """
    Tab page for specifying the ASL scan to optimize
    """

    SCAN_TYPES = {
        "Multi-PLD PCASL" : VAR_MULTI_PCASL,
        #"Multi-TE PCASL" : VAR_TE_PCASL,
        "Look-Locker" : LOOK_LOCKER,
        #"VAR_TE_PCASL_NPLD" : VAR_TE_PCASL_NPLD,
    }

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")

        self.section("ASL scan to optimize for")
        self._asltype = self.choice("ASL scan type", choices=list(self.SCAN_TYPES.keys()))
        self._fvalue = self.number("Estimated perfusion (ml/100g/min)", minval=0, maxval=100.0, initial=50.0)
        self._duration = self.number("Approximate scan duration (s)", minval=0, maxval=1000, initial=300)
        self._readout_time = self.number("Readout time (s)", minval=0, maxval=2.0, initial=0.5)
        self._tau = self.number("Bolus duration (s)", minval=0, maxval=5, initial=1.4)
        self._nplds = self.integer("Number of PLDs", minval=1, maxval=20, initial=6)

        self.readout_ch = self.choice("Readout", choices=["3D (eg GRASE)", "2D multi-slice (eg EPI)"])
        self.readout_ch.SetSelection(0)
        self.readout_ch.Bind(wx.EVT_CHOICE, self._readout_changed)

        self._nslices = self.integer("Number of slices", minval=1, maxval=100, initial=10)
        self._slicedt = self.number("Time per slice (ms)", minval=0, maxval=50, step=1, initial=10)
        self._readout_changed()

        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)

    def _readout_changed(self, _event=None):
        self._nslices.Enable(self.readout_ch.GetSelection() == 1)
        self._slicedt.Enable(self.readout_ch.GetSelection() == 1)

    def slicedt(self):
        if self.readout_ch.GetSelection() == 1:
            return self._slicedt.GetValue() / 1000
        else:
            return 0.0

    def nslices(self):
        if self.readout_ch.GetSelection() == 1:
            return self._nslices.GetValue()
        else:
            return 1

    def aslparams(self):
        return ASLParams(f=self._fvalue.GetValue()/6000.0)

    def scan(self):
        return ASLScan(self.SCAN_TYPES[self._asltype.GetString(self._asltype.GetSelection())],
                       duration=self._duration.GetValue(),
                       tau=self._tau.GetValue(),
                       npld=self._nplds.GetValue(), 
                       readout=self._readout_time.GetValue(),
                       slices=self.nslices(),
                       slicedt=self.slicedt())

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
        return ATTDist(self._att_start.GetValue(), self._att_end.GetValue(), 
                       self._att_step.GetValue(), self._att_taper.GetValue())

    def pldlimits(self):
        return Limits(self._pld_min.GetValue(), self._pld_max.GetValue(), self._pld_step.GetValue())

    def optimizer(self, *args, **kwargs):
        if self._opttype.GetSelection() == 1:
            return LOptimal([[1, 0],  [0, 0]], *args, **kwargs)
        if self._opttype.GetSelection() == 2:
            return LOptimal([[0, 0],  [0, 1]], *args, **kwargs)
        else:
            return DOptimal(*args, **kwargs)
