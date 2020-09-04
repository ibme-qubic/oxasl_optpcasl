#!/usr/bin/env python
"""
OXASL_OPTPCASL: Main GUI window

Copyright (c) 2020 University of Oxford
"""

import os
import sys

import wx

from .scan_options import ScanOptions
from .phys_params import PhysParamOptions
from .optimizer_options import OptimizerOptions
from .scan_summary import ScanSummary
from .sensitivity_plot import CBFSensitivityPlot, ATTSensitivityPlot, KineticCurve

from ..kinetic_model import BuxtonPcasl
from ..optimize import Optimizer

class OptPCASLGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="OXASL PCASL Optimizer", size=(1100, 700), style=wx.DEFAULT_FRAME_STYLE)
        self._panel = wx.Panel(self)

        local_dir = os.path.abspath(os.path.dirname(__file__))
        icon = wx.Icon()
        icon.CopyFromBitmap(wx.Bitmap(os.path.join(local_dir, "icon.png"), wx.BITMAP_TYPE_ANY))
        self.SetIcon(icon)

        main_vsizer = wx.BoxSizer(wx.VERTICAL)

        hpanel = wx.Panel(self._panel)
        hpanel.SetBackgroundColour((180, 189, 220))
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        banner = wx.Panel(self._panel, size=(-1, 80))
        banner.SetBackgroundColour((180, 189, 220))
        banner_fname = os.path.join(local_dir, "banner.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        hsizer.Add(banner)

        # Dumb hack to make banner images align to left and right
        spacer = wx.StaticText(hpanel, label="")
        spacer.SetBackgroundColour((180, 189, 220))
        hsizer.Add(spacer, 10)

        banner = wx.Panel(self._panel, size=(-1, 80))
        banner.SetBackgroundColour((180, 189, 220))
        banner_fname = os.path.join(local_dir, "oxasl.png")
        wx.StaticBitmap(banner, -1, wx.Bitmap(banner_fname, wx.BITMAP_TYPE_ANY))
        hsizer.Add(banner)

        hpanel = wx.Panel(self._panel)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hpanel.SetSizer(hsizer)
        main_vsizer.Add(hpanel, 2, wx.EXPAND)

        notebook = wx.Notebook(hpanel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        hsizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        notebook.win = self

        self._protocol = ScanOptions(notebook, 0, 1)
        notebook.AddPage(self._protocol, "Scan protocol")

        self._phys_params = PhysParamOptions(notebook, 0, 1)
        notebook.AddPage(self._phys_params, "Physiological parameters")

        self._opt = OptimizerOptions(notebook, 0, 1)
        notebook.AddPage(self._opt, "Optimization")

        notebook = wx.Notebook(hpanel, id=wx.ID_ANY, style=wx.BK_DEFAULT)
        hsizer.Add(notebook, 1, wx.ALL|wx.EXPAND, 5)
        notebook.win = self
        
        self._ss = ScanSummary(notebook)
        notebook.AddPage(self._ss, "Scan summary")
        
        self._cbf = CBFSensitivityPlot(notebook)
        notebook.AddPage(self._cbf, "CBF sensitivity")

        self._att = ATTSensitivityPlot(notebook)
        notebook.AddPage(self._att, "ATT sensitivity")

        self._curve = KineticCurve(notebook)
        notebook.AddPage(self._curve, "Kinetic curve")

        self._panel.SetSizer(main_vsizer)
        self.Layout()

    def set_scan(self):
        """
        """
        phys_params = self._phys_params.get()
        kinetic_model = BuxtonPcasl(phys_params)
        protocol = self._protocol.get(kinetic_model,  self._opt)
        params = protocol.initial_params()
        for plot in (self._att, self._curve, self._cbf, self._ss):
            plot.set(phys_params, protocol, params)
    
    def optimize(self, niters=1):
        """
        """
        phys_params = self._phys_params.get()
        kinetic_model = BuxtonPcasl(phys_params)
        protocol = self._protocol.get(kinetic_model,  self._opt)
        params = protocol.initial_params()
        opt = Optimizer(protocol)
        output = opt.optimize(params, niters)
        for plot in (self._att, self._curve, self._cbf, self._ss):
            plot.set(phys_params, protocol, output["params"])
    
    def _dorun(self, _event):
        try:
            self._run_btn.Enable(False)
            self._panel.SetCursor(wx.Cursor(wx.CURSOR_WAIT))
            wx.Yield()
            params = self._scan_options.aslparams()
            scan = self._scan_options.scan()
            att_dist = self._opt_options.attdist()
            lims = self._opt_options.pldlimits()
            optimizer = self._opt_options.optimizer(params, scan, att_dist, lims)
            
            # Run the optimisation and display outcome
            output = optimizer.optimize()
            self._plot.set_optimized_scan(params, scan, output)
            
        except (RuntimeError, ValueError) as exc:
            sys.stderr.write("ERROR: %s\n" % str(exc))
            import traceback
            traceback.print_exc()
        finally:
            self._run_btn.Enable(True)
            self._panel.SetCursor(wx.Cursor(wx.CURSOR_ARROW))

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
