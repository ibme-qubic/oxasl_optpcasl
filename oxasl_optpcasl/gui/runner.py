"""
OXASL_OPTPCASL - Code to run optimization, and display output

Copyright (C) 2020 University of Oxford
"""

from threading import Thread

import wx
from wx.lib.pubsub import pub

from ..optimize import Optimizer

class OptimizationRunner(wx.Frame):
    """
    Runs oxford_asl (and any other required commands) based on configured options
    """

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Run", size=(600, 400), style=wx.DEFAULT_FRAME_STYLE)
        self._optimizer_thread = None

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.output_text = wx.TextCtrl(self, style=wx.TE_READONLY | wx.TE_MULTILINE)
        font = wx.Font(10, wx.TELETYPE, wx.NORMAL, wx.NORMAL)
        self.output_text.SetFont(font)
        self.sizer.Add(self.output_text, 2, flag=wx.EXPAND)

        self.cancel = wx.Button(self, label="Cancel")
        self.cancel.Bind(wx.EVT_BUTTON, self._cancel)
        self.sizer.Add(self.cancel, 0)
        self.cancel.Disable()

        self.SetSizer(self.sizer)
        self.Bind(wx.EVT_CLOSE, self._close)
        pub.subscribe(self._write_output, "opt_log")
        pub.subscribe(self._finished, "opt_finished")

    def _write_output(self, line):
        self.output_text.AppendText(line)

    def _close(self, _event):
        self.Hide()

    def _cancel(self, _event):
        if self._optimizer_thread is not None:
            self._optimizer_thread.optimizer.cancel = True
            self.cancel.Disable()

    def _finished(self, output):
        if output is None:
            self._write_output("\nWARNING: optimization failed\n")
        self.output = output
        self._optimizer_thread = None
        self.cancel.Disable()

    def run(self, protocol, cost_model, **kwargs):
        self.Show()
        self.Raise()
        self.cancel.Enable()
        self.output_text.Clear()
        self._optimizer_thread = OptimizationThread(protocol, cost_model, **kwargs)
        self._optimizer_thread.start()

class OptimizationThread(Thread):
    """
    Run the optimization in the background
    """
    def __init__(self, protocol, cost_model, **kwargs):
        """
        :param cmds: Sequence of ``Cmd`` instances
        """
        Thread.__init__(self)
        self.optimizer = Optimizer(protocol, cost_model, log=self)
        self.kwargs = kwargs

    def run(self):
        output = None
        try:
            gridpts = self.kwargs.pop("gridpts", None)
            if gridpts is not None:
                self.kwargs["initial_params"] = self.optimizer.gridsearch(gridpts)
            output = self.optimizer.optimize(**self.kwargs)
        except Exception as exc:
            print(exc)
            wx.CallAfter(pub.sendMessage, "opt_log", line=str(exc))
        finally:
            wx.CallAfter(pub.sendMessage, "opt_finished", output=output)

    def write(self, line):
        wx.CallAfter(pub.sendMessage, "opt_log", line=line)