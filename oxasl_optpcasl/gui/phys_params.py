"""
OXASL_OPTPCASL: Widget to control physiological parameters

Copyright (c) 2019 University of Oxford
"""
from ..structures import PhysParams
from .widgets import TabPage

class PhysParamOptions(TabPage):
    """
    Tab page for specifying the ASL scan to optimize
    """

    def __init__(self, parent, idx, n):
        TabPage.__init__(self, parent, "Params", idx, n, name="params")

        self.section("Physiological parameters")
        self._f = self.number("Estimated perfusion (ml/100g/min)", digits=3, minval=0, maxval=100.0, initial=50.0)
        self._t1t = self.number("T1 (tissue)", digits=3, minval=0, maxval=5, initial=1.445)
        self._t1b = self.number("T1 (blood)", digits=3, minval=0, maxval=5, initial=1.65)
        self._alpha = self.number("Inversion efficiency", digits=3, minval=0, maxval=1, initial=0.85)
        self._lam = self.number("Partition coefficient", digits=3, minval=0, maxval=1.0, initial=0.9)

        self.sizer.AddGrowableCol(1, 1)
        self.SetSizer(self.sizer)

    def get(self):
        return PhysParams(f=self._f.GetValue() / 6000, t1t=self._t1t.GetValue(), t1b=self._t1b.GetValue(), 
                          alpha=self._alpha.GetValue(), lam=self._lam.GetValue())
