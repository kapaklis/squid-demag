"""Draw icon.png (1024 px) and icon.icns: a sheared raw loop and a sharp corrected loop on an indigo tile."""
import os
import subprocess
import tempfile

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))

fig = Figure(figsize=(10.24, 10.24), dpi=100)
fig.patch.set_alpha(0)
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, 1024), ylim=(0, 1024))
ax.set_axis_off()

# macOS icon grid: 824 px tile, 185 px corner radius, centred on the 1024 canvas
tile = FancyBboxPatch((100 + 185, 100 + 185), 824 - 370, 824 - 370, boxstyle="round,pad=185",
                      facecolor="none", edgecolor="none")
ax.add_patch(tile)
grad = np.linspace(0, 1, 256)[:, None]
img = ax.imshow(grad, extent=(100, 924, 100, 924), origin="lower", aspect="auto",
                cmap=LinearSegmentedColormap.from_list("tile", ["#312e81", "#4f46e5", "#7c7cf6"]))
img.set_clip_path(tile)

# faint axes
ax.plot([512, 512], [215, 809], color="white", alpha=0.18, lw=6, solid_capstyle="round")
ax.plot([215, 809], [512, 512], color="white", alpha=0.18, lw=6, solid_capstyle="round")

h = np.linspace(-1, 1, 400)
cx, cy, sx, sy = 512, 512, 250, 205
for shear, alpha, lw in [(0.5, 0.32, 24), (0.0, 1.0, 32)]:
    for s in (1, -1):
        m = np.tanh((h - s * 0.4) / 0.2)  # ascending branch switches at +Hc, descending at -Hc
        H = h + shear * m * 0.45  # demagnetizing shear tilts the raw loop
        ax.plot(cx + sx * H, cy + sy * m, color="white", alpha=alpha, lw=lw, solid_capstyle="round",
                solid_joinstyle="round")

png = os.path.join(HERE, "icon.png")
fig.savefig(png, dpi=100, transparent=True)

with tempfile.TemporaryDirectory() as d:
    iconset = os.path.join(d, "icon.iconset")
    os.mkdir(iconset)
    for size in (16, 32, 128, 256, 512):
        for scale in (1, 2):
            px = size * scale
            name = f"icon_{size}x{size}{'@2x' if scale == 2 else ''}.png"
            subprocess.run(["sips", "-z", str(px), str(px), png, "--out", os.path.join(iconset, name)],
                           check=True, capture_output=True)
    subprocess.run(["iconutil", "-c", "icns", iconset, "-o", os.path.join(HERE, "icon.icns")], check=True)
print("wrote icon.png and icon.icns")
