"""tkinter front end: pick a curve and a fuel mixture, get FDS input lines."""

import math
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .curves import CURVE_TYPES, GROWTH_RATES
from .fuels import FUELS
from .writer import DesignFire

CUSTOM = "custom"


def nice_step(span, target=5):
    """Round axis tick step giving about ``target`` ticks over ``span``."""
    raw = span / target
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


class CurvePlot(tk.Canvas):
    PAD_L, PAD_R, PAD_T, PAD_B = 60, 15, 15, 40

    def __init__(self, master, **kw):
        super().__init__(master, background="white", highlightthickness=0, **kw)
        self.curve = None
        self.bind("<Configure>", lambda e: self.redraw())

    def show(self, curve):
        self.curve = curve
        self.redraw()

    def redraw(self):
        self.delete("all")
        w, h = self.winfo_width(), self.winfo_height()
        if self.curve is None or w < 100 or h < 80:
            return
        x0, x1 = self.PAD_L, w - self.PAD_R
        y0, y1 = h - self.PAD_B, self.PAD_T
        t_max = self.curve.duration
        q_max = self.curve.peak_hrr * 1.1
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
        xy = []
        for i in range(n + 1):
            t = t_max * i / n
            xy += [sx(t), sy(self.curve.hrr(t))]
        self.create_line(*xy, fill="#c0392b", width=2)
        # ramp points as written to FDS
        for t, f in self.curve.points():
            x, y = sx(t), sy(f * self.curve.peak_hrr)
            self.create_oval(x - 2, y - 2, x + 2, y + 2, outline="#2c3e50")


class App(ttk.Frame):
    def __init__(self, master):
        super().__init__(master, padding=8)
        self.master = master
        master.title("FDS Design Fire Generator")
        self.grid(sticky="nsew")
        master.columnconfigure(0, weight=1)
        master.rowconfigure(0, weight=1)

        self.curve_type = tk.StringVar(value=next(iter(CURVE_TYPES)))
        self.growth = tk.StringVar(value="fast")
        self.alpha = tk.StringVar()
        self.peak = tk.StringVar(value="2000")
        self.area = tk.StringVar(value="4")
        self.duration = tk.StringVar(value="1200")
        self.per_fuel = tk.BooleanVar(value=False)
        self.new_fuel = tk.StringVar(value="Polyurethane foam (flexible)")
        self.new_fraction = tk.StringVar(value="1")
        self.status = tk.StringVar()
        self.mixture = []        # [(fuel name, mass fraction)]

        self._build()
        self._on_growth()
        self.add_fuel()
        for v in (self.curve_type, self.alpha, self.peak, self.area,
                  self.duration, self.per_fuel):
            v.trace_add("write", lambda *a: self.refresh())
        self.growth.trace_add("write", lambda *a: self._on_growth())
        self.refresh()

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
        rows = [
            ("Curve type", ttk.Combobox(box, textvariable=self.curve_type,
                                        values=list(CURVE_TYPES),
                                        state="readonly", width=14)),
            ("Growth rate", ttk.Combobox(box, textvariable=self.growth,
                                         values=list(GROWTH_RATES) + [CUSTOM],
                                         state="readonly", width=14)),
            ("alpha (kW/s²)", ttk.Entry(box, textvariable=self.alpha, width=16)),
            ("Peak HRR (kW)", ttk.Entry(box, textvariable=self.peak, width=16)),
            ("Fire area (m²)", ttk.Entry(box, textvariable=self.area, width=16)),
            ("Duration (s)", ttk.Entry(box, textvariable=self.duration, width=16)),
        ]
        for i, (label, widget) in enumerate(rows):
            ttk.Label(box, text=label).grid(row=i, column=0, sticky="w", pady=2)
            widget.grid(row=i, column=1, sticky="ew", pady=2)
        self.alpha_entry = rows[2][1]

        box = ttk.LabelFrame(left, text="Fuel mixture (mass fractions)",
                             padding=6)
        box.grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        left.rowconfigure(1, weight=1)
        self.tree = ttk.Treeview(box, columns=("fuel", "w"), show="headings",
                                 height=6)
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

    # actions ------------------------------------------------------------
    def _on_growth(self):
        g = self.growth.get()
        if g == CUSTOM:
            self.alpha_entry.state(["!disabled"])
        else:
            self.alpha_entry.state(["!disabled"])
            self.alpha.set(f"{GROWTH_RATES[g]:.5g}")
            self.alpha_entry.state(["disabled"])
        self.refresh()

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

    def build_fire(self):
        def num(var, label):
            try:
                return float(var.get())
            except ValueError:
                raise ValueError(f"{label} must be a number") from None

        curve_cls = CURVE_TYPES[self.curve_type.get()]
        curve = curve_cls(num(self.alpha, "alpha"), num(self.peak, "Peak HRR"),
                          num(self.duration, "Duration"))
        if not self.mixture:
            raise ValueError("Add at least one fuel")
        return DesignFire(curve, num(self.area, "Fire area"),
                          [(FUELS[n], w) for n, w in self.mixture],
                          per_fuel_reactions=self.per_fuel.get())

    def refresh(self):
        if not hasattr(self, "text"):
            return
        try:
            fire = self.build_fire()
            fds = fire.to_fds()
        except ValueError as e:
            self.status.set(str(e))
            return
        self.status.set("")
        self.plot.show(fire.curve)
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


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
