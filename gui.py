"""GUI for demag.py: load a SQUID loop, set the background and specimen shape, see and save the corrected loop."""
import os
import sys
import tkinter as tk
import tkinter.font as tkfont
from math import cos, radians, sin
from tkinter import filedialog, ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from demag import (SHAPES, N_along, N_in_plane, correct, diagonal_angle, find_column, fit_background, pick_moment,
                   read_info, read_loop, read_table, reference_background, unit_warning, volume, write_csv)

HERE = os.path.dirname(os.path.abspath(__file__))
BG, CARD, BORDER, HOVER = "#f3f4f7", "#ffffff", "#e3e5ea", "#f1f2f6"
TEXT, MUTED = "#1d2129", "#6b7280"
ACCENT, ACCENT_DARK, ACCENT_SOFT = "#4f46e5", "#4338ca", "#e9e9fd"
FIELD = "#e5484d"
SIZE_LABELS = {"prism": "Sides (mm)", "ellipsoid": "Axis lengths (mm)",
               "cylinder": "Diameter A = B, length C (mm)  ·  axis along C"}
FILETYPES = [("MPMS data or CSV", "*.dat *.DAT *.csv *.CSV *.txt"), ("All files", "*")]
INFO_KEYS = {"SAMPLE_MATERIAL": "{}", "SAMPLE_COMMENT": "{}", "SAMPLE_SIZE": "size {}", "SAMPLE_SHAPE": "shape {}",
             "SAMPLE_MASS": "{} g"}


def apply_style(root):
    """Flat light theme on top of ttk's 'clam', which (unlike 'aqua') lets every colour be set."""
    family = tkfont.nametofont("TkDefaultFont").actual("family")
    body, small = (family, 13), (family, 11)
    s = ttk.Style(root)
    s.theme_use("clam")
    root.configure(bg=BG)
    s.configure(".", background=CARD, foreground=TEXT, font=body, bordercolor=BORDER, lightcolor=CARD,
                darkcolor=CARD, troughcolor=BG, focuscolor=ACCENT, selectbackground=ACCENT, selectforeground="white")
    s.configure("App.TFrame", background=BG)
    s.configure("Card.TFrame", background=CARD, relief="solid", borderwidth=1, bordercolor=BORDER)
    s.configure("App.TLabel", background=BG)
    s.configure("Title.TLabel", background=BG, font=(family, 18, "bold"))
    s.configure("Subtitle.TLabel", background=BG, foreground=MUTED, font=small)
    s.configure("Heading.TLabel", font=(family, 13, "bold"))
    s.configure("Muted.TLabel", foreground=MUTED, font=small)
    s.configure("Stat.TLabel", font=(family, 19, "bold"))
    s.configure("StatSmall.TLabel", font=(family, 13, "bold"))
    s.configure("Error.TLabel", foreground=FIELD, font=small)
    s.configure("Footnote.TLabel", background=BG, foreground="#a0a5b0", font=(family, 10))

    for w in ("TEntry", "TCombobox"):
        s.configure(w, fieldbackground="white", bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    padding=(6, 4), insertcolor=TEXT, arrowcolor=MUTED)
        s.map(w, bordercolor=[("focus", ACCENT)], lightcolor=[("focus", ACCENT)],
              fieldbackground=[("disabled", HOVER), ("readonly", "white")], foreground=[("disabled", MUTED)],
              selectbackground=[("readonly", "white")],
              selectforeground=[("readonly", TEXT)])
    root.option_add("*TCombobox*Listbox.background", "white")
    root.option_add("*TCombobox*Listbox.selectBackground", ACCENT)
    root.option_add("*TCombobox*Listbox.font", body)

    s.configure("TButton", background=CARD, bordercolor=BORDER, lightcolor=CARD, darkcolor=CARD,
                padding=(10, 5), relief="flat")
    s.map("TButton", background=[("active", HOVER)])
    s.configure("Accent.TButton", background=ACCENT, foreground="white", bordercolor=ACCENT, lightcolor=ACCENT,
                darkcolor=ACCENT, padding=(14, 7), font=(family, 13, "bold"))
    s.map("Accent.TButton", background=[("disabled", "#c9cbf3"), ("active", ACCENT_DARK)],
          bordercolor=[("disabled", "#c9cbf3")], lightcolor=[("disabled", "#c9cbf3")],
          darkcolor=[("disabled", "#c9cbf3")], foreground=[("disabled", "white")])

    # segmented control: radio buttons drawn as buttons
    s.configure("Seg.Toolbutton", background=HOVER, bordercolor=BORDER, lightcolor=HOVER, darkcolor=HOVER,
                padding=(10, 5), anchor="center", relief="flat")
    s.map("Seg.Toolbutton", background=[("selected", ACCENT), ("active", ACCENT_SOFT)],
          foreground=[("selected", "white")], lightcolor=[("selected", ACCENT)], darkcolor=[("selected", ACCENT)],
          bordercolor=[("selected", ACCENT)])

    s.configure("TRadiobutton", background=CARD, indicatorbackground="white", indicatorforeground="white",
                upperbordercolor=BORDER, lowerbordercolor=BORDER, padding=(0, 3))
    s.map("TRadiobutton", indicatorbackground=[("selected", ACCENT)], background=[("active", CARD)],
          upperbordercolor=[("selected", ACCENT)], lowerbordercolor=[("selected", ACCENT)])
    s.configure("TSeparator", background=BORDER)


def length(mm):
    """Side length with a readable unit: 5 mm, 200 µm, 50 nm."""
    if mm >= 0.1:
        return f"{mm:g} mm"
    if mm >= 1e-4:
        return f"{mm * 1e3:.3g} µm"
    return f"{mm * 1e6:.3g} nm"


class App(tk.Tk):
    def __init__(self, path=None):
        super().__init__()
        self.title("SQUID Demag")
        self.geometry("1400x860")
        self.minsize(1150, 760)
        try:
            self.icon = tk.PhotoImage(file=os.path.join(HERE, "icon.png"))
            self.iconphoto(True, self.icon)
        except tk.TclError:
            pass  # icon is cosmetic
        apply_style(self)
        self.header, self.table, self.corrected, self.ref_path, self.path = [], [], [], None, None

        # header bar
        head = ttk.Frame(self, style="App.TFrame", padding=(20, 16, 20, 8))
        head.pack(fill="x")
        titles = ttk.Frame(head, style="App.TFrame")
        titles.pack(side="left")
        ttk.Label(titles, text="SQUID loop correction", style="Title.TLabel").pack(anchor="w")
        ttk.Label(titles, text="Demagnetizing field and substrate background", style="Subtitle.TLabel").pack(anchor="w")
        ttk.Button(head, text="Open loop…", style="Accent.TButton", command=self.open).pack(side="right")
        self.file_label = ttk.Label(head, text="No file loaded", style="Subtitle.TLabel")
        self.file_label.pack(side="right", padx=12)

        ttk.Label(self, text="© 2026 Vassilios Kapaklis  ·  vassilios.kapaklis@physics.uu.se  ·  MIT License", style="Footnote.TLabel").pack(
            side="bottom", anchor="e", padx=22, pady=(0, 8))
        body = ttk.Frame(self, style="App.TFrame", padding=(20, 8, 20, 6))
        body.pack(fill="both", expand=True)
        body.columnconfigure(1, weight=1)
        body.rowconfigure(0, weight=1)
        left = ttk.Frame(body, style="App.TFrame", width=310)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 16))
        right = ttk.Frame(body, style="App.TFrame", width=300)
        right.grid(row=0, column=2, sticky="ns", padx=(16, 0))

        self.build_data_card(left)
        self.build_sample_card(left)
        self.build_background_card(left)
        self.build_plot_card(body)
        self.build_geometry_card(right)
        self.build_results_card(right)

        for var in [self.field_col, self.moment_col, self.shape, *self.dims, self.axis, self.angle, self.N_override,
                    self.V_override, self.bg_mode, self.bg_hmin, self.bg_slope, self.bg_scale]:
            var.trace_add("write", lambda *_: self.update_loop())
        self.update_loop()
        if path:
            self.load(path)  # after the first draw, so a load error isn't overwritten

    # ---- layout -------------------------------------------------------------------------------------------------

    def card(self, parent, title):
        c = ttk.Frame(parent, style="Card.TFrame", padding=16)
        c.pack(fill="x", pady=(0, 14))
        ttk.Label(c, text=title, style="Heading.TLabel").grid(row=0, column=0, columnspan=6, sticky="w", pady=(0, 10))
        return c

    def build_data_card(self, parent):
        c = self.card(parent, "Data")
        self.field_col, self.moment_col = tk.StringVar(), tk.StringVar()
        self.combos = []
        for r, (text, var) in enumerate([("Field column (Oe)", self.field_col), ("Moment column (emu)", self.moment_col)]):
            ttk.Label(c, text=text, style="Muted.TLabel").grid(row=1 + 2 * r, column=0, sticky="w")
            cb = ttk.Combobox(c, textvariable=var, state="readonly", width=30)
            cb.grid(row=2 + 2 * r, column=0, sticky="ew", pady=(2, 8))
            self.combos.append(cb)
        self.info_label = ttk.Label(c, text="", style="Muted.TLabel", wraplength=260, justify="left")
        self.info_label.grid(row=5, column=0, sticky="w")

    def build_sample_card(self, parent):
        c = self.card(parent, "Sample")
        self.shape = tk.StringVar(value="prism")
        shapes = ttk.Frame(c)
        shapes.grid(row=1, column=0, columnspan=6, sticky="w", pady=(0, 10))
        for value in SHAPES:
            ttk.Radiobutton(shapes, text=value.capitalize(), value=value, variable=self.shape,
                            style="Seg.Toolbutton").pack(side="left")
        self.size_label = ttk.Label(c, text=SIZE_LABELS["prism"], style="Muted.TLabel")
        self.size_label.grid(row=2, column=0, columnspan=6, sticky="w")
        self.dims = [tk.StringVar(value=v) for v in ("5", "5", "0.5")]
        self.dim_entries = []
        for i, (name, var) in enumerate(zip("ABC", self.dims)):
            ttk.Label(c, text=name).grid(row=3, column=2 * i, sticky="e", padx=(0 if i == 0 else 8, 4))
            e = ttk.Entry(c, textvariable=var, width=6)
            e.grid(row=3, column=2 * i + 1, sticky="w", pady=(2, 10))
            self.dim_entries.append(e)

        ttk.Label(c, text="Field along", style="Muted.TLabel").grid(row=4, column=0, columnspan=6, sticky="w")
        self.axis = tk.IntVar(value=2)
        seg = ttk.Frame(c)
        seg.grid(row=5, column=0, columnspan=6, sticky="w", pady=(2, 8))
        for value, text in enumerate(["A", "B", "C", "A–B plane"]):
            ttk.Radiobutton(seg, text=text, value=value, variable=self.axis, style="Seg.Toolbutton").pack(side="left")

        plane = ttk.Frame(c)
        plane.grid(row=6, column=0, columnspan=6, sticky="w", pady=(0, 10))
        ttk.Label(plane, text="φ from A").pack(side="left")
        self.angle = tk.StringVar(value="45")
        ttk.Entry(plane, textvariable=self.angle, width=6).pack(side="left", padx=(6, 2))
        ttk.Label(plane, text="°").pack(side="left")
        ttk.Button(plane, text="Diagonal", command=self.set_diagonal).pack(side="left", padx=(10, 0))

        ttk.Label(c, text="Overrides  ·  blank = from the shape", style="Muted.TLabel").grid(
            row=7, column=0, columnspan=6, sticky="w")
        self.N_override, self.V_override = tk.StringVar(), tk.StringVar()
        over = ttk.Frame(c)
        over.grid(row=8, column=0, columnspan=6, sticky="w", pady=(2, 0))
        ttk.Label(over, text="N").pack(side="left")
        ttk.Entry(over, textvariable=self.N_override, width=7).pack(side="left", padx=(4, 12))
        ttk.Label(over, text="V (cm³)").pack(side="left")
        ttk.Entry(over, textvariable=self.V_override, width=9).pack(side="left", padx=4)

    def build_background_card(self, parent):
        c = self.card(parent, "Substrate background")
        c.columnconfigure(0, weight=1)
        self.bg_mode = tk.StringVar(value="none")
        self.bg_hmin, self.bg_slope, self.bg_scale = tk.StringVar(value="10000"), tk.StringVar(), tk.StringVar(value="1")
        for r, (value, text, var) in enumerate([("none", "None", None),
                                                ("fit", "Linear fit,  |H| ≥ (Oe)", self.bg_hmin),
                                                ("slope", "Slope χ (emu/Oe)", self.bg_slope),
                                                ("file", "Bare substrate  ×", self.bg_scale)]):
            ttk.Radiobutton(c, text=text, value=value, variable=self.bg_mode).grid(row=1 + r, column=0, sticky="w")
            if var:
                ttk.Entry(c, textvariable=var, width=9).grid(row=1 + r, column=1, sticky="e", pady=2)
        ref = ttk.Frame(c)
        ref.grid(row=5, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        ttk.Button(ref, text="Substrate file…", command=self.open_reference).pack(side="left")
        self.ref_label = ttk.Label(ref, text="none", style="Muted.TLabel")
        self.ref_label.pack(side="left", padx=8)

    def build_plot_card(self, parent):
        c = ttk.Frame(parent, style="Card.TFrame", padding=12)
        c.grid(row=0, column=1, sticky="nsew")
        self.fig = Figure(figsize=(7, 5), layout="constrained", facecolor=CARD)
        self.ax = self.fig.add_subplot()
        self.canvas = FigureCanvasTkAgg(self.fig, c)
        bar = NavigationToolbar2Tk(self.canvas, c, pack_toolbar=False)
        for w in [bar, *bar.winfo_children()]:
            try:
                w.configure(background=CARD, highlightthickness=0)
            except tk.TclError:
                pass
        bar.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().configure(highlightthickness=0)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def build_geometry_card(self, parent):
        c = self.card(parent, "Geometry")
        self.geo_fig = Figure(figsize=(3.0, 2.7), facecolor=CARD)
        self.geo_ax = self.geo_fig.add_axes([0, 0, 1, 1], projection="3d", facecolor=CARD)
        self.geo_canvas = FigureCanvasTkAgg(self.geo_fig, c)
        self.geo_canvas.get_tk_widget().configure(highlightthickness=0)
        self.geo_canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew")
        self.geo_note = ttk.Label(c, text="", style="Muted.TLabel")
        self.geo_note.grid(row=2, column=0, sticky="w", pady=(6, 0))

    def build_results_card(self, parent):
        c = self.card(parent, "Results")
        c.columnconfigure((0, 1), weight=1, uniform="stat")
        self.stats = {}
        for i, (key, label, style) in enumerate([("N", "Demag factor N", "Stat.TLabel"),
                                                 ("V", "Volume (cm³)", "Stat.TLabel"),
                                                 ("bg", "Background", "StatSmall.TLabel"),
                                                 ("pts", "Points", "Stat.TLabel")]):
            tile = ttk.Frame(c)
            tile.grid(row=1 + i // 2, column=i % 2, sticky="w", pady=(0, 10))
            ttk.Label(tile, text=label, style="Muted.TLabel").pack(anchor="w")
            self.stats[key] = tk.StringVar(value="—")
            ttk.Label(tile, textvariable=self.stats[key], style=style).pack(anchor="w")
        self.message = ttk.Label(c, text="", style="Error.TLabel", wraplength=250, justify="left")
        self.message.grid(row=3, column=0, columnspan=2, sticky="w")
        self.save_button = ttk.Button(parent, text="Save corrected CSV…", style="Accent.TButton", command=self.save)
        self.save_button.pack(fill="x")

    # ---- actions ------------------------------------------------------------------------------------------------

    def open(self):
        path = filedialog.askopenfilename(filetypes=FILETYPES)
        if path:
            self.load(path)

    def load(self, path):
        try:
            self.header, self.table = read_table(path)
        except (OSError, ValueError) as e:
            self.message.config(text=f"Could not read file: {e}", style="Error.TLabel")
            return
        self.path = path
        self.file_label.config(text=os.path.basename(path))
        for cb in self.combos:
            cb.config(values=self.header)
        info = read_info(path)  # MPMS3 sample notes, to help fill in the dimensions
        self.info_label.config(text="  ·  ".join(
            INFO_KEYS[k].format(info[k]) for k in INFO_KEYS if info.get(k)))
        for var, pick in [(self.field_col, lambda: find_column(self.header, "Magnetic Field")),
                          (self.moment_col, lambda: pick_moment(self.header, self.table))]:
            try:
                var.set(self.header[pick()])
            except ValueError:
                var.set("")

    def open_reference(self):
        path = filedialog.askopenfilename(filetypes=FILETYPES)
        if path:
            self.ref_path = path
            self.ref_label.config(text=os.path.basename(path))
            self.bg_mode.set("file")
            self.update_loop()

    def background(self, rows):
        """(function of H, short description) for the chosen background."""
        mode = self.bg_mode.get()
        if mode == "fit":
            chi = fit_background(rows, float(self.bg_hmin.get()))
            return (lambda H: chi * H), f"χ = {chi:.3g} emu/Oe (fit)"
        if mode == "slope":
            chi = float(self.bg_slope.get())
            return (lambda H: chi * H), f"χ = {chi:.3g} emu/Oe"
        if mode == "file":
            if not self.ref_path:
                raise ValueError("open a substrate file")
            # the reference picks its own moment column, so a VSM substrate run can correct a DC-scan sample
            ref = read_loop(self.ref_path, self.field_col.get())[2]
            return reference_background(ref, float(self.bg_scale.get())), f"substrate × {float(self.bg_scale.get()):g}"
        return (lambda H: 0.0), "none"

    def set_diagonal(self):
        try:
            dims = [float(v.get()) for v in self.dims]
        except ValueError:
            return
        self.angle.set(f"{diagonal_angle(dims):.4g}")
        self.axis.set(3)

    def save(self):
        stem = os.path.splitext(os.path.basename(self.path))[0]
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile=f"{stem}_corrected.csv",
                                            initialdir=os.path.dirname(self.path))
        if path:
            write_csv(path, self.corrected)
            self.message.config(text=f"Saved {os.path.basename(path)}", style="Muted.TLabel")

    # ---- drawing ------------------------------------------------------------------------------------------------

    def update_loop(self):
        self.corrected = []
        self.save_button.state(["disabled"])
        self.message.config(text="", style="Error.TLabel")
        shape = self.shape.get()
        self.size_label.config(text=SIZE_LABELS[shape])
        self.dim_entries[1].state(["disabled"] if shape == "cylinder" else ["!disabled"])
        if shape == "cylinder" and self.dims[1].get() != self.dims[0].get():
            self.dims[1].set(self.dims[0].get())  # B mirrors the diameter; the trace redraws
            return
        try:
            dims = [float(v.get()) for v in self.dims]
            if min(dims) <= 0:
                raise ValueError("sides must be positive")
            phi = float(self.angle.get()) if self.axis.get() == 3 else None
            if self.N_override.get().strip():
                N = float(self.N_override.get())
            elif phi is not None:
                N = N_in_plane(dims, phi, shape)
            else:
                N = N_along(dims, self.axis.get(), shape)
            V = float(self.V_override.get()) if self.V_override.get().strip() else volume(dims, shape)
            if not 0 <= N <= 1 or V <= 0:
                raise ValueError("need 0 ≤ N ≤ 1 and volume > 0")
        except ValueError as e:
            self.message.config(text=f"Check the inputs: {e}")
            return
        self.draw_geometry(dims, phi)
        self.stats["N"].set(f"{N:.4f}")
        self.stats["V"].set(f"{V:.3g}")

        rows = []
        if self.field_col.get() and self.moment_col.get():
            hi, mi = self.header.index(self.field_col.get()), self.header.index(self.moment_col.get())
            for r in self.table:
                try:
                    rows.append((float(r[hi]), float(r[mi])))
                except (ValueError, IndexError):
                    continue  # MPMS leaves failed measurements blank
        self.stats["pts"].set(str(len(rows)) if rows else "—")

        ax = self.ax
        ax.clear()
        ax.set_facecolor(CARD)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(BORDER)
        ax.tick_params(colors=MUTED, labelsize=10)
        ax.grid(color="#eef0f3", lw=0.8)
        ax.set_axisbelow(True)
        ax.set_xlabel("H (Oe)", color=MUTED)
        ax.set_ylabel("M (emu/cm³)", color=MUTED)
        if rows:
            try:
                background, note = self.background(rows)
                self.corrected = correct(rows, N, V, background)
            except (OSError, ValueError) as e:
                self.message.config(text=f"Background: {e}")
                self.stats["bg"].set("—")
                self.canvas.draw_idle()
                return
            self.stats["bg"].set(note)
            warning = unit_warning(self.field_col.get(), self.moment_col.get())
            if warning:
                self.message.config(text=warning.capitalize())
            H, m, _, M, Hi = zip(*self.corrected)
            if self.bg_mode.get() != "none":
                ax.plot(H, [x / V for x in m], "--", color="#c3c7cf", lw=1, label="raw")
            ax.plot(H, M, color="#9aa0ab", lw=1.2, label="vs applied field")
            ax.plot(Hi, M, color=ACCENT, lw=1.8, label="vs internal field")
            ax.axhline(0, color=BORDER, lw=0.8)
            ax.axvline(0, color=BORDER, lw=0.8)
            ax.legend(frameon=False, labelcolor=TEXT, fontsize=10)
            self.save_button.state(["!disabled"])
        else:
            self.stats["bg"].set("—")
            ax.text(0.5, 0.5, "Open a loop to start" if not self.header else "No numeric data in these columns",
                    transform=ax.transAxes, ha="center", va="center", color=MUTED, fontsize=13)
            ax.set_xticks([])
            ax.set_yticks([])
        self.canvas.draw_idle()

    def draw_geometry(self, dims, phi):
        """3D sketch of the specimen with the field arrow; drag to rotate."""
        ax = self.geo_ax
        view = (ax.elev, ax.azim) if ax.has_data() else (30, -58)
        ax.clear()
        ax.set_axis_off()
        rel = np.array(dims) / max(dims)
        size = np.maximum(rel, 0.08)  # keep thin films visible
        x, y, z = size / 2
        shape = self.shape.get()
        t = np.linspace(0, 2 * np.pi, 96)
        if shape == "prism":
            v = np.array([[i, j, k] for i in (-x, x) for j in (-y, y) for k in (-z, z)])
            faces = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]
            ax.add_collection3d(Poly3DCollection([v[f] for f in faces], facecolor=ACCENT, alpha=0.16,
                                                 edgecolor=ACCENT, linewidths=1.2))
        elif shape == "ellipsoid":
            u, w = np.meshgrid(t, np.linspace(0, np.pi, 32))
            ax.plot_surface(x * np.cos(u) * np.sin(w), y * np.sin(u) * np.sin(w), z * np.cos(w),
                            color=ACCENT, alpha=0.14, linewidth=0, shade=False)
            zero = 0 * t
            for curve in [(x * np.cos(t), y * np.sin(t), zero), (x * np.cos(t), zero, z * np.sin(t)),
                          (zero, y * np.cos(t), z * np.sin(t))]:
                ax.plot(*curve, color=ACCENT, lw=1)
        else:
            u, h = np.meshgrid(t, [-z, z])
            ax.plot_surface(x * np.cos(u), y * np.sin(u), h, color=ACCENT, alpha=0.14, linewidth=0, shade=False)
            ax.add_collection3d(Poly3DCollection(
                [np.column_stack([x * np.cos(t), y * np.sin(t), np.full_like(t, zz)]) for zz in (-z, z)],
                facecolor=ACCENT, alpha=0.16, edgecolor=ACCENT, linewidths=1.2))
            for zz in (-z, z):
                ax.plot(x * np.cos(t), y * np.sin(t), np.full_like(t, zz), color=ACCENT, lw=1.2)
            for side in (-1, 1):  # silhouette edges
                ax.plot([side * x * np.cos(np.radians(32))] * 2, [side * y * np.sin(np.radians(32))] * 2, [-z, z],
                        color=ACCENT, lw=1)
        kw = dict(color=TEXT, fontsize=9, ha="center", va="center")
        ax.text(0, -y - 0.2, -z, f"{'⌀ A = B' if shape == 'cylinder' else 'A'}  {length(dims[0])}", **kw)
        if shape != "cylinder":
            ax.text(x + 0.2, 0, -z, f"B  {length(dims[1])}", **kw)
        ax.text(-x, -y, z + 0.13, f"C  {length(dims[2])}", **kw)

        axis = self.axis.get()
        d = np.array([cos(radians(phi)), sin(radians(phi)), 0]) if phi is not None else np.eye(3)[axis]
        ax.quiver(*(-0.72 * d), *(1.44 * d), color=FIELD, linewidth=2.4, arrow_length_ratio=0.14)
        ax.text(*(0.86 * d + [0, 0, 0.06]), "H", color=FIELD, fontsize=12, fontweight="bold", ha="center")
        if phi is not None:
            t = np.linspace(0, radians(phi), 30)
            ax.plot(0.34 * np.cos(t), 0.34 * np.sin(t), 0 * t, color=FIELD, lw=1)
            ax.plot([0, 0.55], [0, 0], [0, 0], color=MUTED, lw=0.8, ls=":")
            ax.text(0.44 * cos(radians(phi / 2)), 0.44 * sin(radians(phi / 2)), 0, "φ", color=FIELD, fontsize=10)
        lim = 0.65
        ax.set(xlim=(-lim, lim), ylim=(-lim, lim), zlim=(-lim, lim))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(*view)
        self.geo_note.config(text="Drag to rotate" + ("  ·  not to scale" if (size != rel).any() else ""))
        self.geo_canvas.draw_idle()


if __name__ == "__main__":
    App(sys.argv[1] if len(sys.argv) > 1 else None).mainloop()
