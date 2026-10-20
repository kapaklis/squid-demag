"""Demagnetization and substrate-background correction of SQUID M(H) loops (CGS: Oe, emu).

M = (moment - background(H_applied)) / volume,  H_int = H_applied - 4*pi*N*M.
N for a rectangular prism from A. Aharoni, J. Appl. Phys. 83, 3432 (1998); for an ellipsoid from
J. A. Osborn, Phys. Rev. 67, 351 (1945); for a cylinder (magnetometric) from R. I. Joseph, J. Appl. Phys. 37,
4639 (1966).
"""
import argparse
import csv
import sys
import numpy as np
from math import atan, atan2, cos, degrees, isclose, log, pi, radians, sin, sqrt
from scipy.special import ellipe, elliprd

SHAPES = ("prism", "ellipsoid", "cylinder")
VOLUME_FACTOR = {"prism": 1, "ellipsoid": pi / 6, "cylinder": pi / 4}  # V = factor * A * B * C


def _ln_ratio(s, x):
    # ln((s-x)/(s+x)) written as ln((s^2-x^2)/(s+x)^2) so thin shapes don't lose s-x to cancellation
    return log((s * s - x * x) / (s + x) ** 2)


def prism_N(a, b, c):
    """Demagnetizing factor along c for a prism with sides a, b, c (any length unit)."""
    a, b, c = a / 2, b / 2, c / 2  # Aharoni uses half-sides
    r = sqrt(a * a + b * b + c * c)
    ab, bc, ac = sqrt(a * a + b * b), sqrt(b * b + c * c), sqrt(a * a + c * c)
    pi_N = (
        (b * b - c * c) / (2 * b * c) * _ln_ratio(r, a)
        + (a * a - c * c) / (2 * a * c) * _ln_ratio(r, b)
        - b / (2 * c) * _ln_ratio(ab, a)
        - a / (2 * c) * _ln_ratio(ab, b)
        + c / (2 * a) * _ln_ratio(bc, b)
        + c / (2 * b) * _ln_ratio(ac, a)
        + 2 * atan(a * b / (c * r))
        + (a ** 3 + b ** 3 - 2 * c ** 3) / (3 * a * b * c)
        + (a * a + b * b - 2 * c * c) / (3 * a * b * c) * r
        + c / (a * b) * (ac + bc)
        - (ab ** 3 + bc ** 3 + ac ** 3) / (3 * a * b * c)
    )
    return pi_N / pi


def ellipsoid_N(a, b, c):
    """Demagnetizing factor along a for an ellipsoid with axes a, b, c (Carlson form of Osborn's integral)."""
    return a * b * c / 3 * float(elliprd(b * b, c * c, a * a))


def cylinder_N(length, diameter):
    """Axial magnetometric demagnetizing factor of a cylinder."""
    t = length / diameter
    m, m1 = 1 / (1 + t * t), t * t / (1 + t * t)  # m1 = 1 - m, kept exact for thin discs
    # t^2 K + (1 - t^2) E with K - E = (m/3) R_D(0, 1-m, 1), which avoids K - E cancelling for long rods
    x = ellipe(m) + m1 / 3 * elliprd(0, m1, 1)
    return 1 + 4 / (3 * pi * t) * (1 - sqrt(1 + t * t) * x)


def factors(dims, shape="prism"):
    """(N_A, N_B, N_C) for a body of extent dims=(A, B, C) along its axes.

    prism: sides; ellipsoid: full axis lengths; cylinder: diameter A = B, length C (axis along C).
    """
    A, B, C = dims
    if shape == "prism":
        return tuple(prism_N(*(d for i, d in enumerate(dims) if i != ax), dims[ax]) for ax in range(3))
    if shape == "ellipsoid":
        return ellipsoid_N(A, B, C), ellipsoid_N(B, C, A), ellipsoid_N(C, A, B)
    if shape == "cylinder":
        if not isclose(A, B, rel_tol=1e-9):
            raise ValueError("a cylinder needs A = B (the diameter)")
        Nc = cylinder_N(C, A)
        return (1 - Nc) / 2, (1 - Nc) / 2, Nc  # transverse factors follow from the trace being 1
    raise ValueError(f"unknown shape '{shape}'")


def volume(dims, shape="prism"):
    """Volume in cm^3 from extents in mm."""
    return VOLUME_FACTOR[shape] * dims[0] * dims[1] * dims[2] / 1000


def N_along(dims, axis, shape="prism"):
    """Demagnetizing factor for the field along axis 0, 1 or 2 (A, B, C)."""
    return factors(dims, shape)[axis]


def diagonal_angle(dims):
    """Angle in degrees from A of the A-B diagonal."""
    return degrees(atan2(dims[1], dims[0]))


def N_in_plane(dims, phi_deg, shape="prism"):
    """Factor for the field in the A-B plane at phi degrees from A, assuming M stays parallel to the field.

    The tensor is diag(N_A, N_B, N_C), so N_eff = N_A cos^2 phi + N_B sin^2 phi. When N_A != N_B there is
    also a transverse field (N_A - N_B) sin phi cos phi M that this ignores; it vanishes when A = B.
    """
    c, s = cos(radians(phi_deg)), sin(radians(phi_deg))
    NA, NB, _ = factors(dims, shape)
    return NA * c * c + NB * s * s


def read_table(path):
    """Return (header, rows of strings) from an MPMS .dat or a plain CSV."""
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    if "[Data]" in lines:  # MPMS .dat: header block, then CSV after [Data]
        lines = lines[lines.index("[Data]") + 1:]
    reader = csv.reader(lines)
    header = next(reader, None)
    if not header:
        raise ValueError(f"{path}: no header row")
    return header, list(reader)


def find_column(header, key):
    """Index of the first column whose name contains key (case-insensitive)."""
    hits = [i for i, h in enumerate(header) if key.lower() in h.lower()]
    if not hits:
        raise ValueError(f"No column matching '{key}'. Columns: {header}")
    return hits[0]


def read_loop(path, field_key, moment_key):
    """Return (field column name, moment column name, rows of (H, m)) from an MPMS .dat or a plain CSV."""
    header, table = read_table(path)
    hi, mi = find_column(header, field_key), find_column(header, moment_key)
    rows = []
    for row in table:
        try:
            rows.append((float(row[hi]), float(row[mi])))
        except (ValueError, IndexError):
            continue  # MPMS leaves failed measurements blank
    if not rows:
        raise ValueError("No numeric (field, moment) rows found.")
    return header[hi], header[mi], rows


def fit_background(rows, h_min):
    """Slope chi (emu/Oe) of m = chi*H + Ms*sign(H) + c fitted to the points with |H| >= h_min (saturated)."""
    pts = [(H, m) for H, m in rows if abs(H) >= h_min]
    if len(pts) < 3:
        raise ValueError(f"only {len(pts)} points with |H| >= {h_min} Oe; lower the threshold")
    H, m = np.array(pts).T
    # lstsq copes with points on one side only, where sign(H) and c are degenerate but chi is not
    return float(np.linalg.lstsq(np.column_stack([H, np.sign(H), np.ones_like(H)]), m, rcond=None)[0][0])


def reference_background(ref_rows, scale=1.0):
    """Background function from a measured bare substrate: its moment at H (branches averaged), times scale."""
    H, m = np.array(ref_rows).T
    Hu, idx = np.unique(H, return_inverse=True)
    mu = np.bincount(idx, m) / np.bincount(idx)
    tol = 0.01 * max(abs(Hu[0]), abs(Hu[-1]))  # MPMS end points scatter a little around the set field

    def bg(h):
        if not Hu[0] - tol <= h <= Hu[-1] + tol:
            raise ValueError(f"field {h} Oe is outside the reference range [{Hu[0]}, {Hu[-1]}] Oe")
        return scale * float(np.interp(h, Hu, mu))
    return bg


def correct(rows, N, V, background=lambda H: 0.0):
    """Rows of (H, m) -> rows of (H, m, background, M, H_int)."""
    out = []
    for H, m in rows:
        b = background(H)
        M = (m - b) / V
        out.append((H, m, b, M, H - 4 * pi * N * M))
    return out


def write_csv(path_or_file, corrected):
    out = open(path_or_file, "w", newline="") if isinstance(path_or_file, str) else path_or_file
    w = csv.writer(out)
    w.writerow(["H_applied (Oe)", "moment (emu)", "background (emu)", "M (emu/cm3)", "H_int (Oe)"])
    w.writerows(corrected)
    if out is not path_or_file:
        out.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="MPMS .dat or CSV file with a header row")
    p.add_argument("-o", "--output", help="output CSV (default: stdout)")
    p.add_argument("--dims", nargs=3, type=float, metavar=("A", "B", "C"),
                   help="extent in mm: prism sides (film: lateral, lateral, thickness), ellipsoid axis lengths, "
                        "or cylinder diameter, diameter, length")
    p.add_argument("--shape", choices=SHAPES, default="prism", help="specimen shape (default prism)")
    p.add_argument("--axis", choices="ABC", default="C", help="side the field is applied along (default C)")
    p.add_argument("--angle", metavar="DEG",
                   help="field in the A-B plane at DEG degrees from A, or 'diag' for the A-B diagonal; overrides --axis")
    p.add_argument("--N", type=float, help="use this demagnetizing factor instead of the shape formula")
    p.add_argument("--volume", type=float, help="magnetic volume in cm^3 (default: from --dims)")
    p.add_argument("--field-col", default="Magnetic Field", help="substring of the field column name (Oe)")
    p.add_argument("--moment-col", default="Moment", help="substring of the moment column name (emu)")
    bg = p.add_mutually_exclusive_group()
    bg.add_argument("--bg-fit", type=float, metavar="HMIN",
                    help="subtract a linear background fitted to the points with |H| >= HMIN Oe")
    bg.add_argument("--bg-slope", type=float, metavar="CHI", help="subtract CHI*H, CHI in emu/Oe; write a negative value as --bg-slope=-3e-8")
    bg.add_argument("--bg-file", metavar="PATH", help="subtract a measured bare-substrate loop (same column names)")
    p.add_argument("--bg-scale", type=float, default=1.0,
                   help="multiply the --bg-file moment by this, e.g. sample/reference substrate mass (default 1)")
    args = p.parse_args(argv)

    if args.dims and min(args.dims) <= 0:
        p.error("--dims must be positive")
    if args.angle is not None:
        if not args.dims:
            p.error("--angle needs --dims")
        try:
            args.angle = diagonal_angle(args.dims) if args.angle == "diag" else float(args.angle)
        except ValueError:
            p.error("--angle takes degrees or 'diag'")
    if args.N is not None:
        N = args.N
    elif not args.dims:
        p.error("give --dims or --N")
    else:
        try:
            if args.angle is not None:
                N = N_in_plane(args.dims, args.angle, args.shape)
            else:
                N = N_along(args.dims, "ABC".index(args.axis), args.shape)
        except ValueError as e:
            p.error(str(e))
    V = args.volume or (volume(args.dims, args.shape) if args.dims else None)
    if not V or V <= 0:
        p.error("give a positive --volume or --dims")
    if not 0 <= N <= 1:
        p.error(f"N must lie in [0, 1], got {N}")

    try:
        fcol, mcol, rows = read_loop(args.input, args.field_col, args.moment_col)
        background, note = lambda H: 0.0, "none"
        if args.bg_fit is not None:
            chi = fit_background(rows, args.bg_fit)
            background, note = (lambda H: chi * H), f"fit chi = {chi:.6g} emu/Oe"
        elif args.bg_slope is not None:
            background, note = (lambda H: args.bg_slope * H), f"chi = {args.bg_slope:.6g} emu/Oe"
        elif args.bg_file:
            background = reference_background(read_loop(args.bg_file, args.field_col, args.moment_col)[2], args.bg_scale)
            note = f"{args.bg_file} x {args.bg_scale:g}"
        print(f"field: '{fcol}'  moment: '{mcol}'  N = {N:.6f}  V = {V:.6g} cm^3  rows = {len(rows)}  "
              f"background: {note}", file=sys.stderr)
        corrected = correct(rows, N, V, background)
    except (OSError, ValueError) as e:
        sys.exit(str(e))
    write_csv(args.output or sys.stdout, corrected)


if __name__ == "__main__":
    main()
