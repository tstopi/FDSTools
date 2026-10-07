"""tkinter front end: pick a curve and a fuel mixture, get FDS input lines.

The Design fire tab builds one fire. The Library tab lists saved fires,
overlays the ticked ones and works out a fire that envelopes them all.
"""

import math
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk

from .curves import (GROWTH_TIMES, Q_REF, PowerLawCurve, TabulatedCurve)
from .envelope import fit_power_law, pointwise_max
from .fuels import FUELS
from .library import (LibraryEntry, delete_entry, entry_path, load_library,
                      read_csv_curve, save_entry, user_dir)
from .writer import DesignFire

CUSTOM = "custom"
POWER_LAW = PowerLawCurve.name
TABULATED = TabulatedCurve.name
PALETTE = ["#2e86c1", "#28b463", "#d68910", "#8e44ad", "#17a589",
           "#ba4a00", "#2c3e50", "#d4ac0d", "#a93226", "#5d6d7e"]


def nice_step(span, target=5):
    """Round axis tick step giving about ``target`` ticks over ``span``."""
    raw = span / target
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


class CurvePlot(tk.Canvas):
    """HRR(t) of one or more curves on shared axes.

    ``show`` takes a list of dicts with ``curve`` and optional ``color``,
    ``width``, ``dash``, ``label`` and ``markers`` (draw the &RAMP points).
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 60, 15, 15, 40

    def __init__(self, master, **kw):
        super().__init__(master, background="white", highlightthickness=0, **kw)
        self.series = []
        self.bind("<Configure>", lambda e: self.redraw())

    def show(self, series):
        self.series = list(series)
        self.redraw()

    def redraw(self):
        self.delete("all")
        w, h = self.winfo_width(), self.winfo_height()
        if not self.series or w < 100 or h < 80:
            return
        x0, x1 = self.PAD_L, w - self.PAD_R
        y0, y1 = h - self.PAD_B, self.PAD_T
        t_max = max(s["curve"].duration for s in self.series)
        q_max = max(s["curve"].peak_hrr for s in self.series) * 1.1
        sx = lambda t: x0 + (x1 - x0) * t / t_max
        sy = lambda q: y0 - (y0 - y1) * q / q_max

        grid = "#e4e4e4"
        step = nice_step(t_max)
        t = 0.0
        while t <= t_max + 1e-9:
            self.create_line(sx(t), y0, sx(t), y1, fill=grid)
            self.create_text(sx(t), y0 + 4, text=f"{t:g}", anchor="n")
            t += step
        step = nice_step(q_max)
        q = 0.0
        while q <= q_max + 1e-9:
            self.create_line(x0, sy(q), x1, sy(q), fill=grid)
            self.create_text(x0 - 4, sy(q), text=f"{q:g}", anchor="e")
            q += step
        self.create_rectangle(x0, y1, x1, y0, outline="black")
        self.create_text((x0 + x1) / 2, h - 4, text="Time (s)", anchor="s")
        self.create_text(12, (y0 + y1) / 2, text="HRR (kW)", angle=90)

        n = max(int(x1 - x0), 2)
        legend = []
        for s in self.series:
            c = s["curve"]
            color = s.get("color", "#c0392b")
            t_end = c.duration
            xy = []
            for i in range(n + 1):
                t = t_max * i / n
                if t > t_end:
                    break
                xy += [sx(t), sy(c.hrr(t))]
            xy += [sx(t_end), sy(c.hrr(t_end))]
            self.create_line(*xy, fill=color, width=s.get("width", 2),
                             dash=s.get("dash", ""))
            if s.get("markers"):
                for t, f in c.points():
                    x, y = sx(t), sy(f * c.peak_hrr)
                    self.create_oval(x - 2, y - 2, x + 2, y + 2,
                                     outline="#2c3e50")
            if s.get("label"):
                legend.append((s["label"], color, s.get("dash", "")))
        for i, (label, color, dash) in enumerate(legend):
            y = y1 + 10 + 15 * i
            self.create_line(x0 + 8, y, x0 + 28, y, fill=color, width=2,
                             dash=dash)
            self.create_text(x0 + 32, y, text=label, anchor="w")


class DesignTab(ttk.Frame):
    """One design fire: curve, area and fuel mixture to FDS text."""

    def __init__(self, master):
        super().__init__(master, padding=8)
        self.curve_type = tk.StringVar(value=POWER_LAW)
        self.growth = tk.StringVar(value="fast")
        self.exponent = tk.StringVar()
        self.growth_time = tk.StringVar()
        self.q_ref = tk.StringVar()
        self.t_start = tk.StringVar(value="0")
        self.peak = tk.StringVar(value="2000")
        self.area = tk.StringVar(value="4")
        self.duration = tk.StringVar(value="1200")
        self.decay = tk.BooleanVar(value=False)
        self.decay_start = tk.StringVar(value="600")
        self.decay_time = tk.StringVar(value="300")
        self.decay_exp = tk.StringVar(value="1")
        self.table_info = tk.StringVar(value="No table loaded")
        self.per_fuel = tk.BooleanVar(value=False)
        self.new_fuel = tk.StringVar(value="Polyurethane foam (flexible)")
        self.new_fraction = tk.StringVar(value="1")
        self.status = tk.StringVar()
        self.mixture = []        # [(fuel name, mass fraction)]
        self.table = None        # TabulatedCurve when the type is tabulated
        self._quiet = True

        self._build()
        self._on_growth()
        self.add_fuel()
        for v in (self.exponent, self.growth_time, self.q_ref, self.t_start,
                  self.peak, self.area, self.duration, self.decay_start,
                  self.decay_time, self.decay_exp, self.per_fuel):
            v.trace_add("write", lambda *a: self.refresh())
        self.growth.trace_add("write", lambda *a: self._on_growth())
        for v in (self.curve_type, self.decay):
            v.trace_add("write", lambda *a: self._update_states())
        self._quiet = False
        self._update_states()

    # layout -------------------------------------------------------------
    def _build(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8))
        right = ttk.Frame(self)
        right.grid(row=0, column=1, sticky="nsew")

        box = ttk.LabelFrame(left, text="Heat release rate", padding=6)
        box.grid(row=0, column=0, sticky="ew")
        e = lambda var: ttk.Entry(box, textvariable=var, width=16)
        table = ttk.Frame(box)
        ttk.Label(table, textvariable=self.table_info).pack(side="left")
        self.load_button = ttk.Button(table, text="Load CSV…",
                                      command=self.load_csv)
        self.load_button.pack(side="right")
        rows = [
            ("Curve type", ttk.Combobox(box, textvariable=self.curve_type,
                                        values=[POWER_LAW, TABULATED],
                                        state="readonly", width=14)),
            ("Growth rate", ttk.Combobox(box, textvariable=self.growth,
                                         values=list(GROWTH_TIMES) + [CUSTOM],
                                         state="readonly", width=14)),
            ("Exponent n", e(self.exponent)),
            ("Growth time t_g (s)", e(self.growth_time)),
            ("Reference HRR (kW)", e(self.q_ref)),
            ("Start time (s)", e(self.t_start)),
            ("Peak HRR (kW)", e(self.peak)),
            ("Duration (s)", e(self.duration)),
            ("", ttk.Checkbutton(box, text="Decay", variable=self.decay)),
            ("Decay start (s)", e(self.decay_start)),
            ("Decay duration (s)", e(self.decay_time)),
            ("Decay exponent", e(self.decay_exp)),
            ("Table", table),
            ("Fire area (m²)", e(self.area)),
        ]
        for i, (label, widget) in enumerate(rows):
            ttk.Label(box, text=label).grid(row=i, column=0, sticky="w", pady=2)
            widget.grid(row=i, column=1, sticky="ew", pady=2)
        w = [r[1] for r in rows]
        self.preset_entries = w[2:5]
        self.power_law_widgets = w[1:9]
        self.decay_entries = w[9:12]

        box = ttk.LabelFrame(left, text="Fuel mixture (mass fractions)",
                             padding=6)
        box.grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        left.rowconfigure(1, weight=1)
        self.tree = ttk.Treeview(box, columns=("fuel", "w"), show="headings",
                                 height=5)
        self.tree.heading("fuel", text="Fuel")
        self.tree.heading("w", text="Fraction")
        self.tree.column("fuel", width=190)
        self.tree.column("w", width=70, anchor="e")
        self.tree.grid(row=0, column=0, columnspan=3, sticky="nsew")
        box.rowconfigure(0, weight=1)
        ttk.Combobox(box, textvariable=self.new_fuel, values=list(FUELS),
                     state="readonly", width=24).grid(row=1, column=0,
                                                      pady=(6, 0))
        ttk.Entry(box, textvariable=self.new_fraction, width=6).grid(
            row=1, column=1, padx=4, pady=(6, 0))
        ttk.Button(box, text="Add", command=self.add_fuel, width=7).grid(
            row=1, column=2, pady=(6, 0))
        ttk.Button(box, text="Remove selected", command=self.remove_fuel).grid(
            row=2, column=0, sticky="w", pady=(4, 0))
        ttk.Checkbutton(box, text="One &REAC per fuel", variable=self.per_fuel
                        ).grid(row=3, column=0, columnspan=3, sticky="w",
                               pady=(6, 0))

        right.columnconfigure(0, weight=1)
        right.rowconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)
        self.plot = CurvePlot(right, width=520, height=260)
        self.plot.grid(row=0, column=0, sticky="nsew")
        out = ttk.Frame(right)
        out.grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        out.columnconfigure(0, weight=1)
        out.rowconfigure(0, weight=1)
        self.text = tk.Text(out, height=14, wrap="none", font="TkFixedFont")
        self.text.grid(row=0, column=0, sticky="nsew")
        sb = ttk.Scrollbar(out, orient="vertical", command=self.text.yview)
        sb.grid(row=0, column=1, sticky="ns")
        sbx = ttk.Scrollbar(out, orient="horizontal", command=self.text.xview)
        sbx.grid(row=1, column=0, sticky="ew")
        self.text.configure(yscrollcommand=sb.set, xscrollcommand=sbx.set)

        bar = ttk.Frame(right)
        bar.grid(row=2, column=0, sticky="ew", pady=(6, 0))
        bar.columnconfigure(0, weight=1)
        ttk.Label(bar, textvariable=self.status, foreground="#b03a2e").grid(
            row=0, column=0, sticky="w")
        ttk.Button(bar, text="Copy", command=self.copy).grid(row=0, column=1)
        ttk.Button(bar, text="Save…", command=self.save).grid(
            row=0, column=2, padx=(4, 0))

    # state --------------------------------------------------------------
    def _update_states(self):
        tabulated = self.curve_type.get() == TABULATED
        custom = self.growth.get() == CUSTOM
        for wdg in self.power_law_widgets:
            wdg.state(["disabled" if tabulated else "!disabled"])
        for wdg in self.preset_entries:
            if not tabulated and not custom:
                wdg.state(["disabled"])
        on = not tabulated and self.decay.get()
        for wdg in self.decay_entries:
            wdg.state(["!disabled" if on else "disabled"])
        self.load_button.state(["!disabled" if tabulated else "disabled"])
        self.refresh()

    def _on_growth(self):
        g = self.growth.get()
        if g != CUSTOM:
            quiet, self._quiet = self._quiet, True
            self.exponent.set("2")
            self.growth_time.set(f"{GROWTH_TIMES[g]:g}")
            self.q_ref.set(f"{Q_REF:g}")
            self._quiet = quiet
        self._update_states()

    def _set_table(self, curve):
        self.table = curve
        self.table_info.set(f"{len(curve.table)} points, "
                            f"peak {curve.peak_hrr:g} kW")

    # fuels --------------------------------------------------------------
    def add_fuel(self):
        name = self.new_fuel.get()
        try:
            w = float(self.new_fraction.get())
            if w <= 0:
                raise ValueError
        except ValueError:
            self.status.set("Mass fraction must be a positive number")
            return
        self.mixture = [(n, x) for n, x in self.mixture if n != name]
        self.mixture.append((name, w))
        self._fill_tree()
        self.refresh()

    def remove_fuel(self):
        names = {self.tree.item(i, "values")[0] for i in self.tree.selection()}
        self.mixture = [(n, w) for n, w in self.mixture if n not in names]
        self._fill_tree()
        self.refresh()

    def _fill_tree(self):
        self.tree.delete(*self.tree.get_children())
        total = sum(w for _, w in self.mixture) or 1.0
        for n, w in self.mixture:
            self.tree.insert("", "end", values=(n, f"{w / total:.3f}"))

    # curve in and out ---------------------------------------------------
    def build_curve(self):
        if self.curve_type.get() == TABULATED:
            if self.table is None:
                raise ValueError("Load a CSV table for the tabulated curve")
            return self.table

        def num(var, label):
            try:
                return float(var.get())
            except ValueError:
                raise ValueError(f"{label} must be a number") from None

        decay = {}
        if self.decay.get():
            decay = dict(decay_start=num(self.decay_start, "Decay start"),
                         decay_time=num(self.decay_time, "Decay duration"),
                         decay_exponent=num(self.decay_exp, "Decay exponent"))
        return PowerLawCurve(num(self.exponent, "Exponent"),
                             num(self.growth_time, "Growth time"),
                             num(self.peak, "Peak HRR"),
                             num(self.duration, "Duration"),
                             q_ref=num(self.q_ref, "Reference HRR"),
                             t_start=num(self.t_start, "Start time"), **decay)

    def area_value(self):
        try:
            return float(self.area.get())
        except ValueError:
            raise ValueError("Fire area must be a number") from None

    def build_fire(self):
        curve = self.build_curve()
        if not self.mixture:
            raise ValueError("Add at least one fuel")
        return DesignFire(curve, self.area_value(),
                          [(FUELS[n], w) for n, w in self.mixture],
                          per_fuel_reactions=self.per_fuel.get())

    def set_design(self, curve, fuels=None, area=None):
        """Load a curve, and optionally a fuel mixture and area."""
        self._quiet = True
        try:
            if isinstance(curve, TabulatedCurve):
                self._set_table(curve)
                self.curve_type.set(TABULATED)
            else:
                preset = next((k for k, tg in GROWTH_TIMES.items()
                               if curve.exponent == 2 and curve.q_ref == Q_REF
                               and abs(curve.growth_time - tg) < 1e-9), CUSTOM)
                self.growth.set(preset)
                self.exponent.set(f"{curve.exponent:g}")
                self.growth_time.set(f"{curve.growth_time:.6g}")
                self.q_ref.set(f"{curve.q_ref:g}")
                self.t_start.set(f"{curve.t_start:g}")
                self.peak.set(f"{curve.peak_hrr:.6g}")
                self.duration.set(f"{curve.duration:g}")
                self.decay.set(curve.decay_start is not None)
                if curve.decay_start is not None:
                    self.decay_start.set(f"{curve.decay_start:.6g}")
                    self.decay_time.set(f"{curve.decay_time:.6g}")
                    self.decay_exp.set(f"{curve.decay_exponent:g}")
                self.curve_type.set(POWER_LAW)
            if fuels:
                self.mixture = list(fuels)
                self._fill_tree()
            if area:
                self.area.set(f"{area:g}")
        finally:
            self._quiet = False
        self._update_states()

    def load_csv(self):
        path = filedialog.askopenfilename(
            title="Tabulated HRR: time (s), HRR (kW)",
            filetypes=[("CSV or text", "*.csv *.txt *.dat"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            self._set_table(read_csv_curve(path))
        except (OSError, ValueError) as e:
            messagebox.showerror("Import failed", str(e))
            return
        self.refresh()

    # output -------------------------------------------------------------
    def refresh(self):
        if self._quiet or not hasattr(self, "text"):
            return
        try:
            fire = self.build_fire()
            fds = fire.to_fds()
        except ValueError as e:
            self.status.set(str(e))
            return
        self.status.set("")
        self.plot.show([{"curve": fire.curve, "markers": True}])
        self.text.delete("1.0", "end")
        self.text.insert("1.0", fds)

    def copy(self):
        self.clipboard_clear()
        self.clipboard_append(self.text.get("1.0", "end-1c"))
        self.status.set("Copied to clipboard")

    def save(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".fds",
            filetypes=[("FDS input", "*.fds"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "w") as fh:
                fh.write(self.text.get("1.0", "end-1c"))
        except OSError as e:
            messagebox.showerror("Save failed", str(e))


class LibraryTab(ttk.Frame):
    """Saved design fires, an overlay of the ticked ones and their envelope."""

    KEEP = "Design tab mixture"
    METHODS = (("Pointwise maximum", "max"), ("Fitted power law", "fit"))

    def __init__(self, master, design):
        super().__init__(master, padding=8)
        self.design = design
        self.entries = []
        self.ticked = set()          # indices into self.entries
        self.envelope = None
        self.method = tk.StringVar(value="max")
        self.exponent = tk.StringVar(value="2")
        self.q_ref = tk.StringVar(value=f"{Q_REF:g}")
        self.decay_exp = tk.StringVar(value="1")
        self.reaction = tk.StringVar(value=self.KEEP)
        self.info = tk.StringVar()
        self.result = tk.StringVar()
        self._build()
        for v in (self.method, self.exponent, self.q_ref, self.decay_exp):
            v.trace_add("write", lambda *a: self.update_envelope())
        self.reload()

    def _build(self):
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        left = ttk.Frame(self)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 8))
        left.rowconfigure(0, weight=1)

        box = ttk.LabelFrame(left, text="Library (click ☐ to stack)",
                             padding=6)
        box.grid(row=0, column=0, sticky="nsew")
        box.rowconfigure(0, weight=1)
        cols = ("use", "name", "type", "peak", "where")
        self.tree = ttk.Treeview(box, columns=cols, show="headings",
                                 selectmode="browse", height=10)
        for c, text, width, anchor in [
                ("use", "", 30, "center"), ("name", "Name", 170, "w"),
                ("type", "Type", 85, "w"), ("peak", "Peak kW", 65, "e"),
                ("where", "Source", 60, "w")]:
            self.tree.heading(c, text=text)
            self.tree.column(c, width=width, anchor=anchor,
                             stretch=c == "name")
        self.tree.grid(row=0, column=0, columnspan=3, sticky="nsew")
        sb = ttk.Scrollbar(box, orient="vertical", command=self.tree.yview)
        sb.grid(row=0, column=3, sticky="ns")
        self.tree.configure(yscrollcommand=sb.set)
        self.tree.bind("<Button-1>", self._on_click)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._show_info())
        ttk.Label(box, textvariable=self.info, wraplength=400,
                  justify="left").grid(row=1, column=0, columnspan=4,
                                       sticky="w", pady=(4, 0))
        buttons = [
            ("Save design…", self.save_design),
            ("Import CSV…", self.import_csv),
            ("Load into design", self.load_selected),
            ("Delete", self.delete_selected),
            ("Reload", self.reload),
            ("Clear ticks", self.clear_ticks),
        ]
        for i, (text, cmd) in enumerate(buttons):
            ttk.Button(box, text=text, command=cmd).grid(
                row=2 + i // 3, column=i % 3, sticky="ew", pady=(4, 0),
                padx=(0, 4))

        box = ttk.LabelFrame(left, text="Envelope of ticked fires", padding=6)
        box.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        for i, (text, value) in enumerate(self.METHODS):
            ttk.Radiobutton(box, text=text, value=value,
                            variable=self.method).grid(row=0, column=i,
                                                       sticky="w")
        e = lambda var: ttk.Entry(box, textvariable=var, width=10)
        self.fit_entries = [e(self.exponent), e(self.q_ref), e(self.decay_exp)]
        for i, (label, wdg) in enumerate(zip(
                ("Exponent n", "Reference HRR (kW)", "Decay exponent"),
                self.fit_entries)):
            ttk.Label(box, text=label).grid(row=1 + i, column=0, sticky="w",
                                            pady=2)
            wdg.grid(row=1 + i, column=1, sticky="w", pady=2)
        ttk.Label(box, text="Reaction from").grid(row=4, column=0, sticky="w",
                                                  pady=2)
        self.reaction_box = ttk.Combobox(box, textvariable=self.reaction,
                                         state="readonly", width=28)
        self.reaction_box.grid(row=4, column=1, columnspan=2, sticky="w",
                               pady=2)
        ttk.Label(box, textvariable=self.result, wraplength=400,
                  justify="left").grid(row=5, column=0, columnspan=3,
                                       sticky="w", pady=(4, 0))
        bar = ttk.Frame(box)
        bar.grid(row=6, column=0, columnspan=3, sticky="w", pady=(4, 0))
        ttk.Button(bar, text="Send to design tab",
                   command=self.send_envelope).pack(side="left")
        ttk.Button(bar, text="Save to library…",
                   command=self.save_envelope).pack(side="left", padx=(4, 0))

        self.plot = CurvePlot(self, width=520, height=420)
        self.plot.grid(row=0, column=1, sticky="nsew")

    # library list -------------------------------------------------------
    def reload(self):
        names = {self.entries[i].name for i in self.ticked}
        self.entries, errors = load_library()
        self.ticked = {i for i, e in enumerate(self.entries)
                       if e.name in names}
        self._fill_tree()
        if errors:
            messagebox.showwarning("Some library files were skipped",
                                   "\n".join(errors))

    def _fill_tree(self):
        self.tree.delete(*self.tree.get_children())
        for i, e in enumerate(self.entries):
            self.tree.insert("", "end", iid=str(i), values=(
                "☑" if i in self.ticked else "☐", e.name, e.curve.name,
                f"{e.curve.peak_hrr:g}", "built-in" if e.builtin else "user"))
        self._update_reactions()
        self.update_envelope()

    def _on_click(self, event):
        if self.tree.identify_region(event.x, event.y) != "cell":
            return
        if self.tree.identify_column(event.x) != "#1":
            return
        row = self.tree.identify_row(event.y)
        if not row:
            return
        i = int(row)
        self.ticked ^= {i}
        self.tree.set(row, "use", "☑" if i in self.ticked else "☐")
        self._update_reactions()
        self.update_envelope()

    def clear_ticks(self):
        self.ticked.clear()
        self._fill_tree()

    def selected(self):
        sel = self.tree.selection()
        return self.entries[int(sel[0])] if sel else None

    def _show_info(self):
        e = self.selected()
        if e is None:
            self.info.set("")
            return
        fuels = ", ".join(f"{n} {w:g}" for n, w in e.fuels) or "no fuel set"
        area = f", area {e.area:g} m²" if e.area else ""
        self.info.set(f"{e.description}\nFuel: {fuels}{area}\n"
                      f"File: {e.path}")

    def _update_reactions(self):
        names = [self.entries[i].name for i in sorted(self.ticked)
                 if self.entries[i].fuels]
        self.reaction_box["values"] = [self.KEEP] + names
        if self.reaction.get() not in self.reaction_box["values"]:
            self.reaction.set(self.KEEP)

    # envelope -----------------------------------------------------------
    def update_envelope(self):
        fit = self.method.get() == "fit"
        for wdg in self.fit_entries:
            wdg.state(["!disabled" if fit else "disabled"])
        curves = [self.entries[i].curve for i in sorted(self.ticked)]
        self.envelope = None
        series = [{"curve": c, "color": PALETTE[k % len(PALETTE)],
                   "label": self.entries[i].name}
                  for k, (i, c) in enumerate(zip(sorted(self.ticked), curves))]
        try:
            if not curves:
                raise ValueError("Tick library fires to stack them")
            if fit:
                nums = []
                for var, label in ((self.exponent, "Exponent"),
                                   (self.q_ref, "Reference HRR"),
                                   (self.decay_exp, "Decay exponent")):
                    try:
                        nums.append(float(var.get()))
                    except ValueError:
                        raise ValueError(f"{label} must be a number") from None
                self.envelope = fit_power_law(curves, *nums)
            else:
                self.envelope = pointwise_max(curves)
        except ValueError as e:
            self.result.set(str(e))
        else:
            env = self.envelope
            text = f"Envelope peak {env.peak_hrr:.6g} kW"
            if fit:
                text += (f", t_g {env.growth_time:.4g} s "
                         f"(alpha {env.alpha:.4g} kW/s^{env.exponent:g})")
                if env.t_start:
                    text += f", starts at {env.t_start:g} s"
                if env.decay_start is not None:
                    text += (f", decay {env.decay_start:.4g} to "
                             f"{env.t_end:.4g} s")
                else:
                    text += ", no decay"
            else:
                text += f", {len(env.table)} points"
            self.result.set(text)
            series.append({"curve": env, "color": "black", "width": 3,
                           "dash": (6, 4), "label": "Envelope"})
        self.plot.show(series)

    def _reaction_entry(self):
        name = self.reaction.get()
        return next((self.entries[i] for i in self.ticked
                     if self.entries[i].name == name), None)

    def send_envelope(self):
        if self.envelope is None:
            messagebox.showinfo("No envelope", self.result.get())
            return
        src = self._reaction_entry()
        self.design.set_design(self.envelope, src.fuels if src else None)
        self.master.select(self.design)

    def save_envelope(self):
        if self.envelope is None:
            messagebox.showinfo("No envelope", self.result.get())
            return
        src = self._reaction_entry()
        names = ", ".join(self.entries[i].name for i in sorted(self.ticked))
        how = {v: t for t, v in self.METHODS}
        self._save(self.envelope, src.fuels if src else self.design.mixture,
                   f"{how[self.method.get()]} envelope of: {names}")

    # entries ------------------------------------------------------------
    def _save(self, curve, fuels, description=""):
        name = simpledialog.askstring("Save to library", "Name:", parent=self)
        if not name:
            return
        desc = simpledialog.askstring("Save to library", "Description:",
                                      initialvalue=description, parent=self)
        if desc is None:
            return
        try:
            area = self.design.area_value()
        except ValueError:
            area = None
        try:
            entry = LibraryEntry(name, curve, list(fuels), area, desc)
            path = entry_path(name)
            if path.exists() and not messagebox.askyesno(
                    "Overwrite?", f"{path} exists. Replace it?"):
                return
            save_entry(entry)
        except (OSError, ValueError) as e:
            messagebox.showerror("Save failed", str(e))
            return
        self.reload()

    def save_design(self):
        try:
            curve = self.design.build_curve()
        except ValueError as e:
            messagebox.showerror("Design fire is not valid", str(e))
            return
        self._save(curve, self.design.mixture)

    def import_csv(self):
        path = filedialog.askopenfilename(
            title="Tabulated HRR: time (s), HRR (kW)",
            filetypes=[("CSV or text", "*.csv *.txt *.dat"),
                       ("All files", "*.*")])
        if not path:
            return
        try:
            curve = read_csv_curve(path)
        except (OSError, ValueError) as e:
            messagebox.showerror("Import failed", str(e))
            return
        self._save(curve, self.design.mixture, f"Imported from {path}")

    def load_selected(self):
        e = self.selected()
        if e is None:
            return
        self.design.set_design(e.curve, e.fuels or None, e.area)
        self.master.select(self.design)

    def delete_selected(self):
        e = self.selected()
        if e is None:
            return
        if e.builtin:
            messagebox.showinfo("Built-in entry",
                                "Built-in entries cannot be deleted. User "
                                f"entries live in {user_dir()}.")
            return
        if not messagebox.askyesno("Delete", f"Delete '{e.name}'?"):
            return
        try:
            delete_entry(e)
        except OSError as err:
            messagebox.showerror("Delete failed", str(err))
        self.ticked.clear()
        self.reload()


class App(ttk.Notebook):
    def __init__(self, master):
        super().__init__(master)
        master.title("FDS Design Fire Generator")
        self.grid(sticky="nsew")
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)
        self.design = DesignTab(self)
        self.library = LibraryTab(self, self.design)
        self.add(self.design, text="Design fire")
        self.add(self.library, text="Library")


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
