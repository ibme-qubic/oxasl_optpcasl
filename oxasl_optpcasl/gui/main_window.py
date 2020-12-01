#!/usr/bin/env python
"""
OXASL_OPTPCASL: Main GUI window

Copyright (c) 2020 University of Oxford
"""

import os
import sys
import webbrowser

# Check we have WX python because this is an optional dependency
try:
    import wx
except ImportError:
    sys.stderr.write("You need to install wxpython to use the GUI\n")
    sys.stderr.write("Try: pip install wxpython\n")
    sys.stderr.write(" or: conda install wxpython if you are using Conda\n")
    sys.exit(1)

from wx.lib.pubsub import pub

from oxasl.reporting import Report

from .scan_options import ScanOptions
from .phys_params import PhysParamOptions
from .optimizer_options import OptimizerOptions
from .scan_summary import ScanSummary
from .sensitivity_plot import CBFSensitivityPlot, ATTSensitivityPlot, KineticCurve
from .runner import OptimizationRunner

from ..kinetic_model import BuxtonPcasl

class OptPCASLGui(wx.Frame):
    """
    Main GUI window
    """

    def __init__(self):
        wx.Frame.__init__(self, None, title="OXASL PCASL Optimizer", size=(1280, 700), style=wx.DEFAULT_FRAME_STYLE)
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

        self.opt = OptimizerOptions(notebook, 0, 1)
        notebook.AddPage(self.opt, "Optimization")

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

        self._runner = OptimizationRunner(self)
        pub.subscribe(self._opt_finished, "opt_finished")

        self._panel.SetSizer(main_vsizer)
        self.Layout()

    def set_scan(self):
        """
        """
        phys_params = self._phys_params.get()
        kinetic_model = BuxtonPcasl(phys_params)
        protocol = self._protocol.get(kinetic_model,  self.opt)
        params = protocol.initial_params()
        for plot in (self._att, self._curve, self._cbf, self._ss):
            plot.set(phys_params, protocol, params, None)
    
    def optimize(self, niters=1):
        """
        """
        phys_params = self._phys_params.get()
        kinetic_model = BuxtonPcasl(phys_params)
        protocol = self._protocol.get(kinetic_model,  self.opt)
        params = protocol.initial_params()
        self._runner.run(protocol, self.opt.cost_model, initial_params=params, reps=niters, gridpts=self.opt.gridpts)

    def generate_report(self, report_path):
        """
        """
        report = Report("OXASL PCASL optimizer - Protocol report", include_timings=False)
        for plot in (self._ss, self._cbf, self._att, self._curve):
            plot.add_to_report(report)
        report.generate(report_path)
        from webbrowser import open
        webbrowser.open("file://" + os.path.abspath(report_path) + "/index.html")

    def _opt_finished(self, output):
        self.opt._set_btn.Enable()
        if output is not None:
            phys_params = self._phys_params.get()
            kinetic_model = BuxtonPcasl(phys_params)
            protocol = self._protocol.get(kinetic_model,  self.opt)
            for plot in (self._att, self._curve, self._cbf, self._ss):
                plot.set(phys_params, protocol, output["params"], self.opt.cost_model)

def create_app():
    """
    Create the wx App checking for dumb 'framework build' error on Mac

    If we find this error, try to fix it if possible
    """
    try:
        return wx.App(redirect=False)
    except SystemExit as exc:
        if sys.platform == "darwin" and "Framework" in str(exc):
            # Ok this is a mega hack. If we are on Mac and it looks like it is the classic
            # 'framework build on Conda' problem then we will first check the Conda framework
            # build of python is installed, and if it is we will try to patch the launcher
            # script to use it. There is a lot that can go wrong here so try to check for
            # problems and bail out if we find any preferably without doing too much damage
            # first
            path_parts = os.path.normpath(sys.executable).split(os.sep)
            conda_python_app = ""
            for part in path_parts:
                if part == "bin": break
                conda_python_app += part + os.sep
            conda_python_app += "python.app/Contents/MacOS/python"
            if not os.path.exists(conda_python_app):
                sys.stderr.write("You seem to have a Conda installation on Mac. We need\n")
                sys.stderr.write("an extra package to run GUIs on Mac Conda. Please try:\n\n")
                sys.stderr.write("  conda install python.app\n\n")
                sys.stderr.write("then try running this program again\n")
                sys.exit(1)
            else:
                launcher = sys.argv[0]
                sys.stderr.write("WARNING: We have found a known problem with running GUI\n")
                sys.stderr.write("applications on Mac. I will try to fix it now by modifying\n")
                sys.stderr.write("the file:\n\n")
                sys.stderr.write("  %s\n\n" % launcher)
                head, tail = os.path.split(launcher)
                with open(launcher) as launcher_file:
                    launcher_code = launcher_file.readlines()
                try:
                    patched = False
                    with open(launcher, "w") as new_launcher_file:
                        for idx, line in enumerate(launcher_code):
                            if idx == 0 and line.startswith("#!"):
                                new_launcher_file.write("#!%s\n" % conda_python_app)
                                patched = True
                            else:
                                new_launcher_file.write(line)

                    if patched:
                        sys.stderr.write("SUCCESS - trying to execute new launcher\n")
                        sys.exit(os.execv(sys.argv[0], sys.argv))
                    else:
                        sys.stderr.write("FAILED - file listed above does not look like a normal GUI launcher\n")
                        sys.stderr.write("Re-throwing original error - you may need to manually run this\n")
                        sys.stderr.write("progam using pythonw or %s instead of just python\n" % conda_python_app)
                        raise
                except PermissionError:
                    sys.stderr.write("FAILED to modify the launcher file because of insufficient permissions\n")
                    sys.stderr.write("You can try either of the following\n")
                    sys.stderr.write(" - Running this program ONCE using sudo (i.e. sudo %s)\n" % launcher)
                    sys.stderr.write(" - OR editing the file above as an administrator and replacing\n")
                    sys.stderr.write("   the first line with:\n\n")
                    sys.stderr.write("#!%s\n\n" % conda_python_app)
                    sys.exit(1)

def main():
    """
    Program entry point
    """
    app = create_app()
    top = OptPCASLGui()
    top.Show()
    app.MainLoop()

if __name__ == '__main__':
    main()
