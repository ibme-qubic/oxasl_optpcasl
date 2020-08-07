"""
OXASL_OPTPCASL: Widget to control the parameters of the PCASL scan to optimize

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from ..structures import ScanParams, PhysParams, ATTDist, Limits
from ..scan import PROTOCOLS
from .widgets import TabPage, NumberChooser

class ScanOptions(TabPage):
    """
    Tab page for specifying the ASL scan to optimize
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Input Data", idx, n, name="input")

        self.section("ASL scan to optimize for")
        self._asltype = self.choice("Scan protocol", choices=list(PROTOCOLS.keys()))
        self._fvalue = self.number("Estimated perfusion (ml/100g/min)", minval=0, maxval=100.0, initial=50.0)
        self._duration = self.number("Approximate scan duration (s)", minval=0, maxval=1000, initial=300)
        self._readout_time = self.number("Readout time (s)", minval=0, maxval=2.0, initial=0.5)
        self._ld = self.number("Label duration (s)", minval=0, maxval=5, initial=1.4)
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
                       ld=self._ld.GetValue(),
                       npld=self._nplds.GetValue(), 
                       readout=self._readout_time.GetValue(),
                       slices=self.nslices(),
                       slicedt=self.slicedt())
