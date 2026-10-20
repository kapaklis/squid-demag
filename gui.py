"""GUI for demag.py: load a SQUID loop, set the specimen shape, see and save the corrected loop."""
import os
import sys
import tkinter as tk
from tkinter import filedialog, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from demag import N_along, N_in_plane, correct, diagonal_angle, find_column, read_table, write_csv


class App(tk.Tk):
    def __init__(self, path=None):
        super().__init__()
        self.title("SQUID demagnetization correction")
        self.header, self.table, self.corrected = [], [], []

        side = ttk.Frame(self, padding=10)
        side.pack(side="left", fill="y")
        row = 0

        def label(text):
            nonlocal row
            ttk.Label(side, text=text).grid(row=row, column=0, sticky="w", pady=2)

        ttk.Button(side, text="Open data file…", command=self.open).grid(row=row, column=0, columnspan=2, sticky="ew")
        row += 1
        self.file_label = ttk.Label(side, text="no file", width=34)
        self.file_label.grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 8))
        row += 1

        self.field_col, self.moment_col = tk.StringVar(), tk.StringVar()
        self.combos = []
        for text, var in [("Field column (Oe)", self.field_col), ("Moment column (emu)", self.moment_col)]:
            label(text)
            row += 1
            cb = ttk.Combobox(side, textvariable=var, state="readonly", width=32)
            cb.grid(row=row, column=0, columnspan=2, sticky="ew")
            self.combos.append(cb)
            row += 1

        ttk.Separator(side).grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1
        label("Prism sides (mm)")
        row += 1
        self.dims = [tk.StringVar(value=v) for v in ("5", "5", "0.5")]
        for name, var in zip("ABC", self.dims):
            label(f"  {name}")
            ttk.Entry(side, textvariable=var, width=12).grid(row=row, column=1, sticky="w")
            row += 1

        label("Field applied along")
        row += 1
        self.axis = tk.IntVar(value=2)
        axes = ttk.Frame(side)
        axes.grid(row=row, column=0, columnspan=2, sticky="w")
        for i, name in enumerate("ABC"):
            ttk.Radiobutton(axes, text=name, value=i, variable=self.axis).pack(side="left", padx=(0, 10))
        row += 1
        plane = ttk.Frame(side)
        plane.grid(row=row, column=0, columnspan=2, sticky="w")
        ttk.Radiobutton(plane, text="in A–B plane, φ from A (°)", value=3, variable=self.axis).pack(side="left")
        self.angle = tk.StringVar(value="45")
        ttk.Entry(plane, textvariable=self.angle, width=7).pack(side="left", padx=4)
        ttk.Button(plane, text="diagonal", width=7, command=self.set_diagonal).pack(side="left")
        row += 1

        ttk.Separator(side).grid(row=row, column=0, columnspan=2, sticky="ew", pady=8)
        row += 1
        label("Overrides (leave blank to use the prism)")
        row += 1
        self.N_override, self.V_override = tk.StringVar(), tk.StringVar()
        for text, var in [("  N", self.N_override), ("  Volume (cm³)", self.V_override)]:
            label(text)
            ttk.Entry(side, textvariable=var, width=12).grid(row=row, column=1, sticky="w")
            row += 1

        self.result = ttk.Label(side, text="", justify="left")
        self.result.grid(row=row, column=0, columnspan=2, sticky="w", pady=8)
        row += 1
        self.save_button = ttk.Button(side, text="Save corrected CSV…", command=self.save, state="disabled")
        self.save_button.grid(row=row, column=0, columnspan=2, sticky="ew")

        plot = ttk.Frame(self)
        plot.pack(side="right", fill="both", expand=True)
        self.fig = Figure(figsize=(7, 5), layout="constrained")
        self.ax = self.fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(self.fig, plot)
        NavigationToolbar2Tk(self.canvas, plot)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        for var in [self.field_col, self.moment_col, *self.dims, self.axis, self.angle, self.N_override, self.V_override]:
            var.trace_add("write", lambda *_: self.update_loop())
        if path:
            self.load(path)
        self.update_loop()

    def open(self):
        path = filedialog.askopenfilename(filetypes=[("MPMS data or CSV", "*.dat *.csv *.txt"), ("All files", "*")])
        if path:
            self.load(path)

    def load(self, path):
        try:
            self.header, self.table = read_table(path)
        except (OSError, ValueError) as e:
            self.result.config(text=f"Could not read file:\n{e}")
            return
        self.path = path
        self.file_label.config(text=os.path.basename(path))
        for cb in self.combos:
            cb.config(values=self.header)
        for var, key in [(self.field_col, "Magnetic Field"), (self.moment_col, "Moment")]:
            try:
                var.set(self.header[find_column(self.header, key)])
            except ValueError:
                var.set("")

    def set_diagonal(self):
        try:
            dims = [float(v.get()) for v in self.dims]
        except ValueError:
            return
        self.angle.set(f"{diagonal_angle(dims):.4g}")
        self.axis.set(3)

    def update_loop(self):
        self.corrected = []
        self.save_button.config(state="disabled")
        try:
            dims = [float(v.get()) for v in self.dims]
            if min(dims) <= 0:
                raise ValueError("sides must be positive")
            if self.N_override.get().strip():
                N = float(self.N_override.get())
            elif self.axis.get() == 3:
                N = N_in_plane(dims, float(self.angle.get()))
            else:
                N = N_along(dims, self.axis.get())
            V = float(self.V_override.get()) if self.V_override.get().strip() else dims[0] * dims[1] * dims[2] / 1000
            if not 0 <= N <= 1 or V <= 0:
                raise ValueError("need 0 ≤ N ≤ 1 and volume > 0")
        except ValueError as e:
            self.result.config(text=f"Check the inputs: {e}")
            return
        status = f"N = {N:.6f}\nV = {V:.6g} cm³"

        rows = []
        if self.field_col.get() and self.moment_col.get():
            hi, mi = self.header.index(self.field_col.get()), self.header.index(self.moment_col.get())
            for r in self.table:
                try:
                    rows.append((float(r[hi]), float(r[mi])))
                except (ValueError, IndexError):
                    continue  # MPMS leaves failed measurements blank
        self.ax.clear()
        if rows:
            self.corrected = correct(rows, N, V)
            H, _, M, Hi = zip(*self.corrected)
            self.ax.plot(H, M, ".-", color="0.6", ms=3, lw=0.8, label="applied field")
            self.ax.plot(Hi, M, ".-", color="C0", ms=3, lw=0.8, label="internal field")
            self.ax.legend()
            self.save_button.config(state="normal")
            status += f"\n{len(rows)} points"
        elif self.header:
            status += "\nno numeric data in these columns"
        self.ax.set_xlabel("H (Oe)")
        self.ax.set_ylabel("M (emu/cm³)")
        self.ax.axhline(0, color="0.85", lw=0.5)
        self.ax.axvline(0, color="0.85", lw=0.5)
        self.canvas.draw_idle()
        self.result.config(text=status)

    def save(self):
        stem = os.path.splitext(os.path.basename(self.path))[0]
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile=f"{stem}_corrected.csv",
                                            initialdir=os.path.dirname(self.path))
        if path:
            write_csv(path, self.corrected)
            self.result.config(text=self.result.cget("text") + f"\nsaved {os.path.basename(path)}")


if __name__ == "__main__":
    App(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()
