"""
oxasl.gui.widgets.py

Useful wx widgets for building the OXASL_OPTPCASL GUI

Copyright (c) 2019 University of Oxford
"""
import os
import colorsys

import numpy as np

import wx
import wx.grid

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

class PlotOutputPanel(wx.Panel):
    """
    Displays plots illustrating the optimized protocol
    """
    def __init__(self, parent):
        self._params = None
        self._cov_optimized = None

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

    def set_optimized_scan(self, best_plds, scantime, params, att_dist, scan, optimizer):
        self._plds_text.Clear()
        self._plds_text.AppendText(" ".join([str(pld) for pld in best_plds]))

        self._scantime_text.Clear()
        self._scantime_text.AppendText(str(scantime))

        params.pld = best_plds
        params.bat = att_dist.exclude_taper

        self._cov_optimized = np.zeros((scan.slices, len(params.bat), 2, 2))
        for slice_idx in range(scan.slices):
            params.t = params.tau + np.array(params.pld) + scan.slicedt*slice_idx
            self._cov_optimized[slice_idx, ...], _ = optimizer.cov_oedn_ave_floor(params, scan, slice_idx)
        
        self._params = params
        self._update_plot()

    def _update_plot(self, _evt=None):
        if self._params is None:
            return

        self._plot_axes.clear()
        if self._plot_choice.GetSelection() == 0:
            self._plot_axes.set_title("Estimated CBF error")
            self._plot_axes.set_ylabel('SD (ml/100g/min)')
            self._plot_axes.set_xlabel("ATT (s)")
            self._plot_axes.set_ylim(0, 9)
            self._plot_axes.plot(self._params.bat, np.squeeze(np.sqrt(np.abs(self._cov_optimized[..., 0, 0]))),
                                 label="Optimized protocol")
            self._plot_axes.legend()
        else:
            self._plot_axes.set_title("Estimated ATT error")
            self._plot_axes.set_ylabel("SD (s)")
            self._plot_axes.set_xlabel("ATT (s)")
            self._plot_axes.set_ylim(0, 0.25)
            self._plot_axes.plot(self._params.bat, np.squeeze(np.sqrt(np.abs(self._cov_optimized[..., 1, 1]))),
                                 label="Optimized protocol")
            self._plot_axes.legend()

        self._canvas.draw()
        self._canvas.Refresh()

class TabPage(wx.Panel):
    """
    Shared methods used by the various tab pages in the GUI
    """
    def __init__(self, parent, title, idx, n, name=None):
        wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY)
        self.notebook = parent
        self.idx = idx
        self.n = n
        self.sizer = wx.GridBagSizer(vgap=5, hgap=5)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name

    def options(self):
        """
        :return: Dictionary of OXASL options from this tab page
        """
        return {}

    def next_prev(self):
        """
        Add next/previous buttons
        """
        if self.idx < self.n-1:
            self.next_btn = wx.Button(self, label="Next", id=wx.ID_FORWARD)
            self.next_btn.Bind(wx.EVT_BUTTON, self._next)
        else:
            self.next_btn = wx.StaticText(self, label="")

        if self.idx > 0:
            self.prev_btn = wx.Button(self, label="Previous", id=wx.ID_BACKWARD)
            self.prev_btn.Bind(wx.EVT_BUTTON, self._prev)
        else:
            self.prev_btn = wx.StaticText(self, label="")

        self.pack(" ")
        self.sizer.AddGrowableRow(self.row-1, 1)
        self.sizer.Add(self.prev_btn, pos=(self.row, 0), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 1), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(wx.StaticText(self, label=""), pos=(self.row, 2), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_LEFT)
        self.sizer.Add(self.next_btn, pos=(self.row, 3), border=10, flag=wx.ALIGN_CENTRE_VERTICAL | wx.ALIGN_RIGHT)

    def _next(self, _):
        """
        Move to the next tab in the notebook
        """
        self.notebook.SetSelection(self.idx+1)

    def _prev(self, _):
        """
        Move to the previous tab in the notebook
        """
        self.notebook.SetSelection(self.idx-1)

    def pack(self, label, *widgets, **kwargs):
        """
        Add a horizontal line to the tab with a label and series of widgets

        If label is empty, first widget is used instead (usually to provide a checkbox)
        """
        col = 0
        border = kwargs.get("border", 10)
        font = self.GetFont()
        if "size" in kwargs:
            font.SetPointSize(kwargs["size"])
        if kwargs.get("bold", False):
            font.SetWeight(wx.BOLD)

        if label != "":
            text = wx.StaticText(self, label=label)
            text.SetFont(font)
            self.sizer.Add(text, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.LEFT)
            col += 1
        else:
            text = None

        for w in widgets:
            span = (1, 1)
            w.label = text
            if hasattr(w, "span"):
                span = (1, w.span)
            w.SetFont(font)
            w.Enable(col == 0 or kwargs.get("enable", True))
            self.sizer.Add(w, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.EXPAND | wx.LEFT, span=span)
            col += span[1]
        self.row += 1

    def file_picker(self, label, pick_dir=False, handler=None, optional=False, initial_on=False, pack=True, **kwargs):
        """
        Add a file picker to the tab
        """
        if not handler:
            handler = self._changed
        if pick_dir:
            picker = wx.DirPickerCtrl(self, style=wx.DIRP_USE_TEXTCTRL)
            picker.Bind(wx.EVT_DIRPICKER_CHANGED, handler)
        else:
            picker = wx.FilePickerCtrl(self)
            picker.Bind(wx.EVT_FILEPICKER_CHANGED, handler)
        picker.span = 2
        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, handler)
            picker.checkbox = cb
            if pack:
                self.pack("", cb, picker, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, picker, **kwargs)

        return picker

    def choice(self, label, choices, initial=0, optional=False, initial_on=False, handler=None, pack=True, **kwargs):
        """
        Add a widget to choose from a fixed set of options
        """
        if not handler:
            handler = self._changed
        ch = wx.Choice(self, choices=choices)
        ch.SetSelection(initial)
        ch.Bind(wx.EVT_CHOICE, handler)
        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, self._changed)
            ch.checkbox = cb
            if pack:
                self.pack("", cb, ch, enable=initial_on, **kwargs)
        elif pack:
            self.pack(label, ch, **kwargs)
        return ch

    def number(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        if not handler:
            handler = self._changed
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, minval=1, maxval=100, **kwargs):
        """
        Add a widget to choose an integer
        """
        if not handler:
            handler = self._changed
        spin = wx.SpinCtrl(self, min=minval, max=maxval, **kwargs)
        spin.SetValue(kwargs.get("initial", 0))
        spin.Bind(wx.EVT_SPINCTRL, handler)
        if pack:
            self.pack(label, spin)
        return spin

    def checkbox(self, label, initial=False, handler=None, **kwargs):
        """
        Add a simple on/off option
        """
        cb = wx.CheckBox(self, label=label)
        cb.span = 2
        cb.SetValue(initial)
        if handler:
            cb.Bind(wx.EVT_CHECKBOX, handler)
        else: cb.Bind(wx.EVT_CHECKBOX, self._changed)
        self.pack("", cb, **kwargs)
        return cb

    def section(self, label):
        """
        Add a section heading
        """
        self.pack(label, bold=True)

    def _changed(self, _):
        self.update()

    def update(self):
        """
        Update the run module, i.e. when options have changed
        """
        if hasattr(self, "run"):
            self.run.update()
            if hasattr(self, "preview"):
                self.preview.run = self.run

    def image(self, label, fname):
        """
        Check if a specified filename is a valid Nifti image
        """
        if not os.path.exists(fname):
            raise OptionError("%s - no such file or directory" % label)
        try:
            return Image(fname)
        except:
            raise OptionError("%s - invalid image file" % label)

class NumberChooser(wx.Panel):
    """
    Widget for choosing a floating point number
    """

    def __init__(self, parent, label=None, minval=0, maxval=1, initial=0.5, step=0.1, digits=2, changed_handler=None):
        super(NumberChooser, self).__init__(parent)
        self.minval, self.orig_min, self.maxval, self.orig_max = minval, minval, maxval, maxval
        self.handler = changed_handler
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        if label is not None:
            self.label = wx.StaticText(self, label=label)
            self.hbox.Add(self.label, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
        # Set a very large maximum as we want to let the user override the default range
        #self.spin = wx.SpinCtrl(self, min=0, max=100000, initial=initial)
        #self.spin.Bind(wx.EVT_SPINCTRL, self._spin_changed)
        self.spin = wx.SpinCtrlDouble(self, min=0, max=100000, inc=step, initial=initial)
        self.spin.SetDigits(digits)
        self.spin.Bind(wx.EVT_SPINCTRLDOUBLE, self._spin_changed)
        self.slider = wx.Slider(self, value=initial, minValue=0, maxValue=100)
        self.slider.SetValue(100*(initial-self.minval)/(self.maxval-self.minval))
        self.slider.Bind(wx.EVT_SLIDER, self._slider_changed)
        self.hbox.Add(self.slider, proportion=1, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.hbox.Add(self.spin, proportion=0, flag=wx.EXPAND | wx.ALIGN_CENTRE_VERTICAL)
        self.SetSizer(self.hbox)

    def GetValue(self):
        """
        :return: numeric value selected
        """
        return self.spin.GetValue()

    def SetValue(self, val):
        """
        Set the numeric value displayed
        """
        self.spin.SetValue(val)
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))

    def _slider_changed(self, event):
        v = event.GetInt()
        val = self.minval + (self.maxval-self.minval)*float(v)/100
        self.spin.SetValue(val)
        if self.handler:
            self.handler(event)
        event.Skip()

    def _spin_changed(self, event):
        """ If user sets the spin outside the current range, update the slider range
        to match. However if they go back inside the current range, revert to this for
        the slider"""
        val = event.GetValue()
        if val < self.minval:
            self.minval = val
        elif val > self.orig_min:
            self.minval = self.orig_min
        if val > self.maxval:
            self.maxval = val
        elif val < self.orig_max:
            self.maxval = self.orig_max
        self.slider.SetValue(100*(val-self.minval)/(self.maxval-self.minval))
        if self.handler:
            self.handler(event)
        event.Skip()
