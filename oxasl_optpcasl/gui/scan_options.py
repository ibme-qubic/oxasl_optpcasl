"""
OXASL_OPTPCASL: Widget to control the parameters of the PCASL scan to optimize

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from ..structures import ScanParams, PhysParams, ATTDist, Limits
from ..scan import *
from .widgets import TabPage, NumberChooser

PROTOCOLS = [
    {
        "name" : "PCASL",
        "desc" : "Standard label/control PCASL acquisition (single or multi PLD)"
    },
    {
        "name" : "Hadamard",
        "desc" : "Time-encoded pCASL acquisition using a Hadamard encoding pattern"
    },
    {
        "name" : "Free lunch",
        "desc" : "Acquisition using a fixed long labelling bolus with the long PLD filled with Hadamard encoded sub-boli"
    },
]

class ScanOptions(TabPage):
    """
    Tab page for specifying the ASL scan to optimize
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Protocol", idx, n, name="protocol")
        self._errors = {}

        self.section("Scan parameters")
        self._protocol = self.choice("Scan protocol", choices=[p["name"] for p in PROTOCOLS], handler=self._protocol_changed)
        self._duration = self.number("Maximum scan duration (s)", minval=0, maxval=1000, initial=300)
        self._readout_time = self.number("Readout time (s)", minval=0, maxval=2.0, initial=0.5)
        self._readout = self.choice("Readout", choices=["3D (eg GRASE)", "2D multi-slice (eg EPI)"], handler=self._readout_changed)
        self._readout.SetSelection(0)

        self._nslices = self.integer("Number of slices", minval=1, maxval=100, initial=10)
        self._slicedt = self.number("Time per slice (ms)", minval=0, maxval=50, step=1, initial=10)

        self._nplds = self.integer("Number of PLDs", minval=1, maxval=20, initial=1, handler=self._nplds_changed)
        self._plds = self.number_list("Initial PLDs (s)", minval=0.1, maxval=3.0, initial=[2.0], handler=self._plds_changed)

        self._ld = self.choice("Label duration", choices=["Fixed", "Single variable", "Multiple variable (one per PLD)"], handler=self._ld_changed)
        self._lds = self.number_list("Initial label durations (s)", minval=0.1, maxval=3.0, initial=[1.8], handler=self._lds_changed)
        
        self._had_label = self.section("Time encoding")
        self._had_size = self.choice("Hadamard matrix size", choices=["4", "8", "12"], initial=1)
        self._had_ld = self.choice("Sub-boli", choices=["All equal", "T1-adjusted", "Unconstrained"], handler=self._had_ld_changed)
        self._had_lds = self.number_list("Sub-bolus duration (s)", minval=0, maxval=3, initial=[0.5], handler=self._had_lds_changed)
        
        self._set_btn = self.button("Set protocol", handler=self._set)
        self._error = self.text(" ", textcol=wx.TheColourDatabase.Find("RED"), bold=True, span=2)

        self.next_prev()
        self.sizer.AddGrowableCol(1)
        self._readout_changed()
        self._protocol_changed()
        self._ld_changed()

    def _readout_changed(self, _event=None):
        self._nslices.Enable(self._readout.GetSelection() == 1)
        self._slicedt.Enable(self._readout.GetSelection() == 1)

    def _protocol_changed(self, _event=None):
        had = self._protocol.GetSelection() in (1, 2)
        had_method = "Show" if had else "Hide"
        nonhad_method = "Hide" if had else "Show"
        
        for w in (self._had_size, self._had_ld, 
                  self._had_lds, self._had_label):
            getattr(w, had_method)()
            if hasattr(w, "label"):
                getattr(w.label, had_method)()
        
        for w in (self._lds, self._ld):
            getattr(w, nonhad_method)()
            getattr(w.label, nonhad_method)()
        
        self.sizer.Layout()

    def _ld_changed(self, _event=None):
        ld_type = self._ld.GetSelection()
        if ld_type == 0:
            self._lds.label.SetLabel("Label duration (s)")
        elif ld_type == 1:
            self._lds.label.SetLabel("Initial label duration (s)")
        else:
            self._lds.label.SetLabel("Initial label durations (s)")
        self._validate_lds()
        self.sizer.Layout()

    def _lds_changed(self, _event=None):
        self._validate_lds()
    
    def _plds_changed(self, _event=None):
        self._validate_plds()
    
    def _nplds_changed(self, _event=None):
        self._validate_plds()

    def _validate_lds(self):
        self._errors.pop("LD", None)
        
        if self._ld.GetSelection() in (0, 1):
            expected_nlds = 1
        else:
            expected_nlds = self._nplds.GetValue()
        
        try:
            lds = self._lds.GetValue()
            if len(lds) != expected_nlds:
                self._errors["LD"] = "Wrong number of LDs - expected %i" % expected_nlds
        except ValueError:
            self._errors["LD"] = "LDs must be space or comma separated numbers"
        self._update_errors()
        
    def _validate_plds(self):
        self._errors.pop("PLD", None)
        
        expected_nplds = self._nplds.GetValue()
        try:
            plds = self._plds.GetValue()
            if len(plds) != expected_nplds:
                self._errors["PLD"] = "Wrong number of PLDs - expected %i" % expected_nplds
        except ValueError:
            self._errors["PLD"] = "PLDs must be space or comma separated numbers"
        self._update_errors()

    def _had_ld_changed(self, _event=None):
        ld_type = self._had_ld.GetSelection()
        if ld_type == 0:
            self._had_lds.label.SetLabel("Initial sub-bolus duration (s)")
        elif ld_type == 1:
            self._had_lds.label.SetLabel("Initial first sub-bolus duration (s)")
        else:
            self._had_lds.label.SetLabel("Initial sub-bolus durations (s)")
        self._validate_had_lds()
        self.sizer.Layout()

    def _had_lds_changed(self, _event=None):
        self._validate_had_lds()
           
    def _validate_had_lds(self):
        self._errors.pop("HADLD", None)
        
        if self._had_ld.GetSelection() in (0, 1):
            expected_nlds = 1
        else:
            expected_nlds = self.had_size-1
        
        try:
            lds = self._had_lds.GetValue()
            if len(lds) != expected_nlds:
                self._errors["HADLD"] = "Wrong number of sub-bolus durations - expected %i" % expected_nlds
        except ValueError:
            self._errors["HADLD"] = "Sub-bolus durations must be space or comma separated numbers"
        self._update_errors()
        
    def _update_errors(self):
        errors = "\n".join(["ERROR: %s" % v for v in self._errors.values()])
        self._error.SetLabel(errors)
        if self._errors:
            self._set_btn.Disable()
        else:
            self._set_btn.Enable()

    @property
    def slicedt(self):
        if self._readout.GetSelection() == 1:
            return self._slicedt.GetValue() / 1000
        else:
            return 0.0

    @property
    def nslices(self):
        if self._readout.GetSelection() == 1:
            return self._nslices.GetValue()
        else:
            return 1

    @property
    def had_size(self):
        return int(self._had_size.GetString(self._had_size.GetSelection()))

    @property
    def scan_params(self):
        if self._protocol.GetSelection() == 0:
            ld = self._lds.GetValue()
        else:
            ld = self._had_lds.GetValue()

        return ScanParams(
            duration=self._duration.GetValue(),
            npld=self._nplds.GetValue(),
            nslices=self.nslices,
            slicedt=self.slicedt,
            readout=self._readout_time.GetValue(),
            ld=ld,
            #noise=,
            plds=self._plds.GetValue(),
            had_size=self.had_size,
        )

    def _set(self, _event=None):
        self.notebook.win.set_scan()

    def get(self, kinetic_model, opt):
        protocol_type = self._protocol.GetSelection()
        ld_type = self._ld.GetSelection()
        had_ld_type = self._had_ld.GetSelection()
        if protocol_type == 0:
            # PCASL
            scan_class = [
                FixedLDPcaslProtocol,
                MultiPLDPcaslVarLD,
                MultiPLDPcaslMultiLD,
            ][ld_type]
        elif protocol_type == 1:
            # Hadamard
            scan_class = [
                HadamardSingleLd,
                HadamardT1Decay,
                HadamardMultiLd,
            ][had_ld_type]
        elif protocol_type == 2:
            # Free lunch
            scan_class = HadamardFreeLunch
        else:
            raise ValueError("Protocol type: %i" % protocol_type)

        return scan_class(kinetic_model,  opt.cost_model, self.scan_params,  opt.att_dist,  opt.pld_lims,  opt.ld_lims)
