"""
OXASL_OPTPCASL: Widget that displays the sensitivity of the protocol to CBF and ATT
Copyright (c) 2019 University of Oxford
"""
import numpy as np

import wx

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

from .widgets import NumberChooser
from ..kinetic_model import BuxtonPcasl

class SensitivityPlot(wx.Panel):
    """
    Displays plots illustrating the optimized protocol
    """
    def __init__(self, parent):
        self._params = None
        self._phys_params = None

        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        self._figure = Figure(figsize=(3.5, 3.5), dpi=100, facecolor='white')
        self._canvas = FigureCanvas(self, -1, self._figure)
        sizer.Add(self._canvas, 2, border=5, flag=wx.EXPAND | wx.ALL)

        self.Layout()

    def set(self, phys_params, scan, params):
        self._phys_params = phys_params
        self._scan = scan
        self._params = params
        self._paramdict = scan.name_params(params)
        self._refresh_plot()

    def _refresh_plot(self):
        self._figure.clf()

        if self._params is not None:
            self._update_plot()

        self._canvas.draw()
        self._canvas.Refresh()

class CBFSensitivityPlot(SensitivityPlot):
    
    def _update_plot(self, _evt=None):
        cov = self._scan.cov(self._params)
        cbf_var = np.squeeze(np.mean(np.sqrt(np.abs(cov[..., 0, 0])), axis=0))
        plot_axes = self._figure.add_subplot(111)
        plot_axes.set_ylim(0, min(20, cbf_var.max()))
        plot_axes.set_title("Estimated CBF std.dev.")
        plot_axes.set_ylabel("CBF std.dev. (ml/100g/min)")
        plot_axes.set_xlabel("ATT (s)")
        plot_axes.plot(self._scan.att_dist.atts, cbf_var, label="")

class ATTSensitivityPlot(SensitivityPlot):
    
    def _update_plot(self, _evt=None):
        cov = self._scan.cov(self._params)
        att_var = np.squeeze(np.mean(np.sqrt(np.abs(cov[..., 1, 1])), axis=0))
        self._plot_axes = self._figure.add_subplot(111)
        self._plot_axes.set_ylim(0, min(att_var.max(), 1.0))
        self._plot_axes.set_title("Estimated ATT std.dev")
        self._plot_axes.set_ylabel("ATT std.dev. (s)")
        self._plot_axes.set_xlabel("ATT (s)")
        self._plot_axes.plot(self._scan.att_dist.atts, att_var, label="")

class KineticCurve(SensitivityPlot):
    
    def __init__(self, parent, model=BuxtonPcasl):
        SensitivityPlot.__init__(self, parent)
        
        self._att_num = NumberChooser(self, label="ATT", minval=0.3, maxval=2.5, initial=1.3, changed_handler=self._att_changed)
        self.GetSizer().Add(self._att_num, border=5, flag=wx.EXPAND | wx.ALL)
        self.Layout()

        self._model = model
        self._att = 1.3
    
    def _update_plot(self, _evt=None):
        times = np.linspace(0, 5.0, 200)
        model = self._model(self._phys_params)
        ld = self._paramdict.get("lds", [self._scan.scan_params.ld])[0]

        lds, plds = self._scan.timings(self._params)
        multi_ld = np.min(lds) != np.max(lds)
        rows = np.ceil(np.sqrt(len(lds)))
        cols = np.ceil(len(lds) / rows)
        idx = 1
        plot_axes = None

        for ld, pld in zip(lds, plds):
            if multi_ld:
                plot_axes = self._figure.add_subplot(rows, cols, idx)
                plot_axes.set_title("LD=%.2f" % ld)
                idx += 1
            
            if multi_ld or plot_axes is None:
                if plot_axes is None:
                    plot_axes = self._figure.add_subplot(111)

                ydata = model.signal(ld, times, self._att)
                plot_axes.plot(times, ydata, linestyle='-', color='blue')
                plot_axes.set_yticklabels([])
                plot_axes.set_ylabel("Relative signal")
                plot_axes.set_xlabel("Time (s)")

            plot_axes.axvline(pld+ld, linestyle='--', color='green')
        self._figure.tight_layout()

    def _att_changed(self, _event=None):
        self._att = self._att_num.GetValue()
        self._refresh_plot()
