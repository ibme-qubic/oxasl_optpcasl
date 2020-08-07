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

class PlotOutputPanel(wx.Panel):
    """
    Displays plots illustrating the optimized protocol
    """
    def __init__(self, parent):
        self._opt_output = None
        self._params = None

        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        plds_panel = wx.Panel(self)
        plds_sizer = wx.BoxSizer(wx.HORIZONTAL)
        plds_panel.SetSizer(plds_sizer)

        font = wx.Font(8, wx.FONTFAMILY_TELETYPE, wx.NORMAL, wx.NORMAL, False)
        label = wx.StaticText(plds_panel, label="Optimized PLDs (s)")
        plds_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._plds_text = wx.TextCtrl(plds_panel, style=wx.TE_READONLY)
        self._plds_text.SetFont(font)
        plds_sizer.Add(self._plds_text, 5, wx.ALL, 5)

        label = wx.StaticText(plds_panel, label="Scan time (s)")
        plds_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._scantime_text = wx.TextCtrl(plds_panel, style=wx.TE_READONLY)
        self._scantime_text.SetFont(font)
        plds_sizer.Add(self._scantime_text, 1, wx.ALL, 5)
        plds_panel.Layout()

        sizer.Add(plds_panel, 0, wx.EXPAND)

        hpanel = wx.Panel(self)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hpanel.SetSizer(hsizer)
        plot_choice_label = wx.StaticText(hpanel, label="Plot selection")
        hsizer.Add(plot_choice_label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self._plot_choice = wx.Choice(hpanel, choices=["CBF error", "ATT error", "Kinetic curve"])
        self._plot_choice.SetSelection(0)
        self._plot_choice.Bind(wx.EVT_CHOICE, self._update_plot)
        hsizer.Add(self._plot_choice, 1, wx.ALL|wx.ALIGN_RIGHT|wx.EXPAND, 5)
        #hpanel.Layout()

        sizer.Add(hpanel, 0)

        figure = Figure(figsize=(3.5, 3.5), dpi=100, facecolor='white')
        self._plot_axes = figure.add_subplot(111)
        figure.subplots_adjust(bottom=0.2)
        self._canvas = FigureCanvas(self, -1, figure)
        sizer.Add(self._canvas, 2, border=5, flag=wx.EXPAND | wx.ALL)

        self.Layout()

    def set_optimized_scan(self, params, scan, opt_output):
        self._opt_output = opt_output
        self._params = params
        self._scan = scan

        self._plds_text.Clear()
        self._plds_text.AppendText(" ".join([str(pld) for pld in opt_output.plds]))

        self._scantime_text.Clear()
        self._scantime_text.AppendText(str(opt_output.scan_time))

        self._update_plot()

    def _update_plot(self, _evt=None):
        if self._opt_output is None:
            return

        cbf_var = np.squeeze(np.mean(np.sqrt(np.abs(self._opt_output.cov_optimized[..., 0, 0])), axis=0))
        att_var = np.squeeze(np.mean(np.sqrt(np.abs(self._opt_output.cov_optimized[..., 1, 1])), axis=0))
        self._plot_axes.clear()
        if self._plot_choice.GetSelection() == 0:
            self._plot_axes.set_title("Estimated CBF error")
            self._plot_axes.set_ylabel('SD (ml/100g/min)')
            self._plot_axes.set_xlabel("ATT (s)")
            self._plot_axes.set_ylim(0, 10)
            self._plot_axes.plot(self._opt_output.att, cbf_var, label="Optimized protocol")
            self._plot_axes.legend()
        elif self._plot_choice.GetSelection() == 1:
            self._plot_axes.set_title("Estimated ATT error")
            self._plot_axes.set_ylabel("SD (s)")
            self._plot_axes.set_xlabel("ATT (s)")
            self._plot_axes.set_ylim(0, 0.25)
            self._plot_axes.plot(self._opt_output.att, att_var, label="Optimized protocol")
            self._plot_axes.legend()
        else:
            self._plot_axes.set_title("ASL kinetic curve")
            self._plot_axes.set_ylabel("Relative signal")
            self._plot_axes.set_xlabel("Time (s)")
            atts = np.linspace(1.0, 1.6, 3)
            for att in atts:
                xdata, ydata = self._kinetic_model(att, self._scan.ld, pldmax=max(self._opt_output.plds))
                self._plot_axes.plot(xdata, ydata, label="ATT=%.2fs" % att)
            for pld in self._opt_output.plds:
                self._plot_axes.axvline(pld+self._scan.ld, linestyle='--', color='green')
            self._plot_axes.legend()

        self._canvas.draw()
        self._canvas.Refresh()

    def _kinetic_model(self, att, ld, f=50.0, lam=0.9, m0=1.0, alpha=0.85, t1b=1.65, t1t=1.445, pldmax=5.0):
        t_all = np.linspace(0, ld+pldmax+1, 50)
        M = np.zeros(len(t_all))
        f = f / 6000 # Fix units
        t1prime = 1/(1.0/t1t + f/lam)

        # During bolus
        relevant_ts = np.logical_and(t_all>att, t_all<=ld+att)
        t = t_all[relevant_ts]
        M[relevant_ts] = np.exp(-att/t1b) * (1-np.exp(-(t-att)/t1prime))

        # After bolus
        relevant_ts = t_all > att+ld
        t = t_all[relevant_ts]
        M[relevant_ts] = np.exp(-att/t1b) * np.exp(-(t-ld-att)/t1prime) * (1-np.exp(-ld/t1prime))

        M *= 2 * m0 * f * t1prime * alpha
        return t_all, M
