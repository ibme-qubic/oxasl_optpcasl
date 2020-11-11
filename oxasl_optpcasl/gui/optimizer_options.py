"""
OXASL_OPTPCASL: Widget to control the optimization process

Copyright (c) 2019 University of Oxford
"""
import os

import wx
import wx.grid

from ..structures import ScanParams, PhysParams, ATTDist, Limits
from ..cost import CBFCost, ATTCost, DOptimalCost
from .widgets import TabPage, NumberChooser

class OptimizerOptions(TabPage):
    """
    Tab page for specifying options for the optimization
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Optimization", idx, n, name="optimization")

        self.section("Optimization type")
        self._opttype = self.choice("Method", choices=["Optimize CBF and ATT", "Optimize CBF only", "Optimize ATT only"])
        
        self.section("ATT prior distribution")
        self._att_start = self.number("Starting value (s)", minval=0, maxval=1.0, initial=0.2, digits=4)
        self._att_end = self.number("Starting value (s)", minval=0, maxval=5.0, initial=2.3, digits=4)
        self._att_step = self.number("Step (s)", minval=0, maxval=0.01, initial=0.001, digits=4)
        self._att_taper = self.number("Taper value (s)", minval=0, maxval=1.0, initial=0.3, digits=4)
        
        self.section("PLD search limits")
        self._pld_min = self.number("Min PLD (s)", minval=0, maxval=1.0, initial=0.075, digits=4)
        self._pld_max = self.number("Max PLD (s)", minval=1.0, maxval=5.0, initial=2.3, digits=4)
        self._pld_step = self.number("Search step (s)", minval=0, maxval=0.1, initial=0.025, digits=4)
        
        self.section("LD search limits")
        self._ld_min = self.number("Min LD (s)", minval=0, maxval=1.0, initial=0.1, digits=4)
        self._ld_max = self.number("Max LD (s)", minval=1.0, maxval=5.0, initial=1.8, digits=4)
        self._ld_step = self.number("Search step (s)", minval=0, maxval=0.1, initial=0.025, digits=4)
        
        self.section("Optimization loop")
        self._niters = self.integer("Number of times to repeat optimization", minval=1, maxval=100, initial=10)
        self._ngridpoints = self.integer("Initial grid search - number of points", minval=100, maxval=1000000, initial=10000, optional=True, handler=self._gridpts_changed)

        self._set_btn = self.button("Optimize", handler=self._optimize)
        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)

    @property
    def att_dist(self):
        return ATTDist(self._att_start.GetValue(), self._att_end.GetValue(), 
                       self._att_step.GetValue(), self._att_taper.GetValue())

    @property
    def pld_lims(self):
        return Limits(self._pld_min.GetValue(), self._pld_max.GetValue(), self._pld_step.GetValue())

    @property
    def ld_lims(self):
        return Limits(self._ld_min.GetValue(), self._ld_max.GetValue(), self._ld_step.GetValue())

    @property
    def gridpts(self):
        if self._ngridpoints.checkbox.GetValue():
            return self._ngridpoints.GetValue()
        else:
            return None

    @property
    def cost_model(self):
        if self._opttype.GetSelection() == 1:
            return CBFCost()
        if self._opttype.GetSelection() == 2:
            return ATTCost()
        else:
            return DOptimalCost()
    
    def _gridpts_changed(self, _event=None):
        if self._ngridpoints.checkbox.GetValue():
            self._ngridpoints.Enable()
        else:
            self._ngridpoints.Disable()

    def _optimize(self, _event=None):
        self._set_btn.Disable()
        self.notebook.win.optimize(self._niters.GetValue())
