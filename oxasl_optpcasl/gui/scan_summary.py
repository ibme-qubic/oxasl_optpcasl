"""
OXASL_OPTPCASL: Widget that displays summary of the scan protocol

Copyright (c) 2019 University of Oxford
"""
import numpy as np

import wx

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

class ScanSummary(wx.Panel):
    """
    Displays plots illustrating the optimized protocol
    """
    def __init__(self, parent):
        self._scan = None
        self._params = None

        wx.Panel.__init__(self, parent, size=wx.Size(300, 600))
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        self._vis = ScanVisualisation(self)
        sizer.Add(self._vis, 2, border=5, flag=wx.EXPAND | wx.ALL)

        plds_panel = wx.Panel(self)
        plds_sizer = wx.BoxSizer(wx.HORIZONTAL)
        plds_panel.SetSizer(plds_sizer)
        sizer.Add(plds_panel, 0, wx.EXPAND)

        monospace_font = wx.Font(8, wx.FONTFAMILY_TELETYPE, wx.NORMAL, wx.NORMAL, False)
        label = wx.StaticText(plds_panel, label="PLDs (s)")
        plds_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._plds_text = wx.TextCtrl(plds_panel, style=wx.TE_READONLY)
        self._plds_text.SetFont(monospace_font)
        plds_sizer.Add(self._plds_text, 5, wx.ALL, 5)

        lds_panel = wx.Panel(self)
        lds_sizer = wx.BoxSizer(wx.HORIZONTAL)
        lds_panel.SetSizer(lds_sizer)
        sizer.Add(lds_panel, 0, wx.EXPAND)

        label = wx.StaticText(lds_panel, label="LDs (s)")
        lds_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._lds_text = wx.TextCtrl(lds_panel, style=wx.TE_READONLY)
        self._lds_text.SetFont(monospace_font)
        lds_sizer.Add(self._lds_text, 5, wx.ALL, 5)

        t_panel = wx.Panel(self)
        t_sizer = wx.BoxSizer(wx.HORIZONTAL)
        t_panel.SetSizer(t_sizer)
        sizer.Add(t_panel, 0, wx.EXPAND)

        label = wx.StaticText(t_panel, label="TR (s)")
        t_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._tr_text = wx.TextCtrl(t_panel, style=wx.TE_READONLY)
        self._tr_text.SetFont(monospace_font)
        t_sizer.Add(self._tr_text, 1, wx.ALL, 5)

        label = wx.StaticText(t_panel, label="Number of repeats")
        t_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._rpts_text = wx.TextCtrl(t_panel, style=wx.TE_READONLY)
        self._rpts_text.SetFont(monospace_font)
        t_sizer.Add(self._rpts_text, 1, wx.ALL, 5)

        label = wx.StaticText(t_panel, label="Total scan time (s)")
        t_sizer.Add(label, 0, wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self._scantime_text = wx.TextCtrl(t_panel, style=wx.TE_READONLY)
        self._scantime_text.SetFont(monospace_font)
        t_sizer.Add(self._scantime_text, 1, wx.ALL, 5)
        
        self.Layout()

    def set(self, phys_params, scan, params):
        self._phys_params = phys_params
        self._scan = scan
        self._params = params

        self._plds_text.Clear()
        self._lds_text.Clear()
        self._tr_text.Clear()
        self._rpts_text.Clear()
        self._scantime_text.Clear()
        
        if self._params is not None:
            paramdict = self._scan.name_params(self._params)
            rpts, tr = self._scan.repeats_total_tr(params)
            self._plds_text.AppendText(" ".join([str(pld) for pld in paramdict.get("plds", [])]))
            self._lds_text.AppendText(" ".join([str(ld) for ld in paramdict.get("lds", [self._scan.scan_params.ld])]))
            self._tr_text.AppendText(str(tr))
            self._rpts_text.AppendText(str(int(rpts)))
            self._scantime_text.AppendText(str(tr * rpts))
            
            desc = self._scan.protocol_summary(params)
            self._vis._summary = desc
            self._vis.Refresh()

class ScanVisualisation(wx.Panel):
    """
    Visual preview of the scan protocol
    """

    def __init__(self, parent):
        wx.Panel.__init__(self, parent, size=wx.Size(300, 150))
        self.SetBackgroundStyle(wx.BG_STYLE_CUSTOM)
        self.Bind(wx.EVT_SIZE, self._on_size)
        self.Bind(wx.EVT_PAINT, self._on_paint)
        self._summary = None
        
        self.hfactor = 0.95
        self.vfactor = 0.95

    def _on_size(self, event):
        event.Skip()
        self.Refresh()

    def _on_paint(self, _event):
        dc = wx.AutoBufferedPaintDC(self)
        dc.Clear()

        if self._summary is None:
            return

        width, height = self.GetClientSize()
        row_height = min(20, 0.8*self.vfactor*height / (len(self._summary) + 4))
        row_width = self.hfactor*width

        ox = width*(1-self.hfactor)/2
        oy = height*(1-self.vfactor)/2

        # Calculate time scale of x axis (pixels per second)
        total_time = 0
        for label, lds, tcs, pld, readout in self._summary:
            total_time = max(total_time, sum(lds) + pld + readout)
        total_time = float(round(total_time*2 + 0.5)) / 2
        scale = row_width / total_time

        y = oy
        self._centered_text(dc, "Protocol schematic", ox + row_width / 2, oy + row_height / 2)

        y += row_height * 2
        for label, lds, tcs, pld, readout in self._summary:
            x = ox
            # Label/control timings
            for ld, tc in zip(lds, tcs):
                rect = wx.Rect(x, y, int(ld * scale)-1, row_height-1)
                col = wx.Colour(128, 128, 255)
                if tc == 1:
                    dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("BLUE"), wx.SOLID))
                else:
                    dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("GREY"), wx.SOLID))
                dc.DrawRectangle(*rect.Get())
                    
                x += ld*scale
            x += pld*scale

            # Readout
            rect = wx.Rect(x, y, readout*scale, row_height-1)
            dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("RED"), wx.SOLID))
            dc.DrawRectangle(*rect.Get())

            y += row_height

        # Scale
        y += 5
        dc.DrawLine(ox, y, ox+row_width, y)
        y += 10
        t = 0.0
        while t < total_time + 0.1:
            x = ox + t * scale
            dc.DrawLine(x, y-5, x, y-10)
            self._centered_text(dc, "%.1f" % t, x, y)
            t += 0.5
        y += row_height - 15

        # Legend
        key_width = row_width / 5
        x = ox + key_width
        y += row_height
        rect = wx.Rect(x, y, 40, row_height-1)
        dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("BLUE"), wx.SOLID))
        dc.DrawRectangle(*rect.Get())
        dc.DrawText("Label", x + 45, y)
        
        x += key_width
        rect = wx.Rect(x, y, 40, row_height-1)
        dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("GREY"), wx.SOLID))
        dc.DrawRectangle(*rect.Get())
        dc.DrawText("Control", x + 45, y)

        x += key_width
        rect = wx.Rect(x, y, 40, row_height-1)
        dc.SetBrush(wx.Brush(wx.TheColourDatabase.Find("RED"), wx.SOLID))
        dc.DrawRectangle(*rect.Get())
        dc.DrawText("Readout", x + 45, y)
        
    def _centered_text(self, dc, text, x, y):
        text_size = dc.GetTextExtent(text)
        dc.DrawText(text, x-text_size.x/2, y-text_size.y/2)
