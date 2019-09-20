#!/usr/bin/env python
"""
Simple wxpython based GUI front-end to OXASL_OPTPCASL

Requirements (beyond OXASL requirements):
    - wxpython
    - matplotlib

Copyright (c) 2019 University of Oxford
"""

import os
import sys

import wx

from .widgets import PlotOutputPanel
from .options import ScanOptions, OptimizerOptions
from .. import optimize, TRWeightingOrNAveFloor

class OptPCASLGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="OXASL PCASL Optimizer", size=(1100, 600), style=wx.DEFAULT_FRAME_STYLE)
        main_panel = wx.Panel(self)
        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        hpanel = wx.Panel(main_panel)
        hpanel.SetBackgroundColour((180, 189, 220))
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((180, 189, 220))
        banner_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        hsizer.Add(banner)

        # Dumb hack to make banner images align to left and right
        spacer = wx.StaticText(hpanel, label="")
        spacer.SetBackgroundColour((180, 189, 220))
        hsizer.Add(spacer, 10)

        banner = wx.Panel(main_panel, size=(-1, 80))
        banner.SetBackgroundColour((180, 189, 220))
        banner_fname = os.path.join(os.path.abspath(os.path.dirname(__file__)), "oxasl.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        hsizer.Add(banner)

        hpanel = wx.Panel(main_panel)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        notebook = wx.Notebook(hpanel, id=wx.ID_ANY, style=wx.BK_DEFAULT)

        self._scan_options = ScanOptions(notebook, 0, 1)
        notebook.AddPage(self._scan_options, "Scan Options")

        self._opt_options = OptimizerOptions(notebook, 0, 1)
        notebook.AddPage(self._opt_options, "Optimizer Options")

        hsizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)

        self._plot = PlotOutputPanel(hpanel)
        hsizer.Add(self._plot, 1, wx.EXPAND)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        run_panel = wx.Panel(main_panel)
        runsizer = wx.BoxSizer(wx.HORIZONTAL)
        run_label = wx.StaticText(run_panel, label="")
        runsizer.Add(run_label, 1, wx.EXPAND)
        run_btn = wx.Button(run_panel, label="Run Optimization")
        run_btn.Bind(wx.EVT_BUTTON, self._dorun)
        runsizer.Add(run_btn, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        run_panel.SetSizer(runsizer)
        main_vsizer.Add(run_panel, 0, wx.EXPAND)
        run_panel.Layout()

        main_panel.SetSizer(main_vsizer)
        self.Layout()

    def _dorun(self, _):
        try:
            
            params = self._scan_options.aslparams()
            scan = self._scan_options.scan()

            att_dist = self._opt_options.attdist()
            lims = self._opt_options.pldlimits()
            optimizer = self._opt_options.optimizer()
            
            # Run the optimisation
            best_plds, num_av, best_min_variance = optimize(optimizer, params, att_dist, scan, lims)
            num_av, total_tr = TRWeightingOrNAveFloor(params, scan, 0, 0)
            scantime = num_av * total_tr
            self._plot.set_optimized_scan(best_plds, scantime, params, att_dist, scan, optimizer)
            
        except (RuntimeError, ValueError) as exc:
            sys.stderr.write("ERROR: %s\n" % str(exc))
            import traceback
            traceback.print_exc()

def main():
    """
    Program entry point
    """
    app = wx.App(redirect=False)
    top = OptPCASLGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
