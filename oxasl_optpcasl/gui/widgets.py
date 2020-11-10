"""
OXASL_OPTPCASL: Useful wx widgets for building the GUI

Copyright (c) 2019 University of Oxford
"""
import os

import wx

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
        self.SetSizer(self.sizer)
        self.row = 0
        self.title = title
        if name is None:
            self.name = title.lower()
        else:
            self.name = name

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
            if "textcol" in kwargs:
                text.SetForegroundColour(kwargs["textcol"])
            if not widgets:
                span = (1, 2)
            else:
                span = (1, 1)
            self.sizer.Add(text, pos=(self.row, col), border=border, flag=wx.ALIGN_CENTRE_VERTICAL | wx.LEFT, span=span)
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
        handler = self._changed_handler(handler)
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
        handler = self._changed_handler(handler)
        ch = wx.Choice(self, choices=choices)
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
        ch.SetSelection(initial)
        return ch

    def button(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        handler = self._changed_handler(handler)
        btn = wx.Button(self, label=label)
        btn.Bind(wx.EVT_BUTTON, handler)
        self.pack("", btn, **kwargs)
        return btn

    def number(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a floating point number
        """
        handler = self._changed_handler(handler)
        num = NumberChooser(self, changed_handler=handler, **kwargs)
        num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def number_list(self, label, handler=None, **kwargs):
        """
        Add a widget to choose a list of floating point numbers
        """
        handler = self._changed_handler(handler)
        num = NumberList(self, changed_handler=handler, **kwargs)
        #num.span = 2
        self.pack(label, num, **kwargs)
        return num

    def integer(self, label, handler=None, pack=True, minval=1, maxval=100, optional=False, initial_on=False, **kwargs):
        """
        Add a widget to choose an integer
        """
        handler = self._changed_handler(handler)
        spin = wx.SpinCtrl(self, min=minval, max=maxval, **kwargs)
        spin.SetValue(kwargs.get("initial", 0))
        spin.Bind(wx.EVT_SPINCTRL, handler)

        if optional:
            cb = wx.CheckBox(self, label=label)
            cb.SetValue(initial_on)
            cb.Bind(wx.EVT_CHECKBOX, handler)
            spin.checkbox = cb
            if pack:
                self.pack("", cb, spin, enable=initial_on)
        elif pack:
            self.pack(label, spin)

        return spin

    def checkbox(self, label, initial=False, handler=None, **kwargs):
        """
        Add a simple on/off option
        """
        handler = self._changed_handler(handler)
        cb = wx.CheckBox(self, label=label)
        #cb.span = 2
        cb.SetValue(initial)
        cb.Bind(wx.EVT_CHECKBOX, handler)
        self.pack("", cb, **kwargs)
        return cb

    def section(self, label):
        """
        Add a section heading
        """
        return self.text(label, bold=True)

    def text(self, label, **kwargs):
        self.pack(label, **kwargs)
        return self.sizer.FindItemAtPosition((self.row-1, 0)).GetWindow()

    def _changed_handler(self, handler):
        def _changed(event):
            if handler:
                handler(event)
            #self.notebook.win.changed()
        return _changed

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
        self.hbox.Add(self.slider, proportion=1, flag=wx.ALIGN_CENTRE_VERTICAL)
        self.hbox.Add(self.spin, proportion=0, flag=wx.ALIGN_CENTRE_VERTICAL)
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

class NumberList(wx.TextCtrl):
    """
    Widget for choosing a list of floating point numbers
    """

    def __init__(self, parent, label=None, minval=0, maxval=1, initial=[0.5], digits=2, changed_handler=None):
        super(NumberList, self).__init__(parent)
        self.minval, self.maxval = minval, maxval
        self.format = "%%.%if" % digits
        self.handler = changed_handler
        self.SetValue(initial)
        self.Bind(wx.EVT_TEXT, self._text_changed)

    def GetValue(self):
        """
        :return: numeric value selected
        """
        return self._validate()

    def SetValue(self, vals):
        """
        Set the numeric value displayed
        """
        try:
            text = ", ".join([self.format % v for v in vals])
        except ValueError:
            text = ""
        wx.TextCtrl.SetValue(self, text)

    def _text_changed(self, event):
        try:
            self._validate()
            self.SetBackgroundColour(wx.NullColour)
        except ValueError as exc:
            self.SetBackgroundColour(wx.Colour(255, 150, 150))
        self.Refresh()

        if self.handler:
            self.handler(event)
        event.Skip()

    def _validate(self):
        try:
            vals = wx.TextCtrl.GetValue(self).replace(",", " ").split()
            floats = []
            for val in vals:
                floats.append(float(val))
            return floats
        except ValueError:
            raise ValueError("%s is not a numeric value" % val)
