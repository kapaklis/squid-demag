"""Demagnetization correction of SQUID M(H) loops (CGS: Oe, emu).

H_int = H_applied - 4*pi*N*M, with M = moment / volume.
N for a rectangular prism from A. Aharoni, J. Appl. Phys. 83, 3432 (1998).
"""
import argparse
import csv
import sys
from math import atan, atan2, cos, degrees, log, pi, radians, sin, sqrt


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


def N_along(dims, axis):
    """Demagnetizing factor of a prism with sides dims=(A, B, C) for the field along axis 0, 1 or 2."""
    return prism_N(*(d for i, d in enumerate(dims) if i != axis), dims[axis])


def diagonal_angle(dims):
    """Angle in degrees from A of the A-B diagonal."""
    return degrees(atan2(dims[1], dims[0]))


def N_in_plane(dims, phi_deg):
    """Factor for the field in the A-B plane at phi degrees from A, assuming M stays parallel to the field.

    The prism's tensor is diag(N_A, N_B, N_C), so N_eff = N_A cos^2 phi + N_B sin^2 phi. When A != B there is
    also a transverse field (N_A - N_B) sin phi cos phi M that this ignores; it vanishes for a square.
    """
    c, s = cos(radians(phi_deg)), sin(radians(phi_deg))
    return N_along(dims, 0) * c * c + N_along(dims, 1) * s * s


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


def correct(rows, N, V):
    """Rows of (H, m) -> rows of (H, m, M, H_int)."""
    return [(H, m, m / V, H - 4 * pi * N * m / V) for H, m in rows]


def write_csv(path_or_file, corrected):
    out = open(path_or_file, "w", newline="") if isinstance(path_or_file, str) else path_or_file
    w = csv.writer(out)
    w.writerow(["H_applied (Oe)", "moment (emu)", "M (emu/cm3)", "H_int (Oe)"])
    w.writerows(corrected)
    if out is not path_or_file:
        out.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="MPMS .dat or CSV file with a header row")
    p.add_argument("-o", "--output", help="output CSV (default: stdout)")
    p.add_argument("--dims", nargs=3, type=float, metavar=("A", "B", "C"),
                   help="prism sides in mm (film: lateral, lateral, thickness)")
    p.add_argument("--axis", choices="ABC", default="C", help="side the field is applied along (default C)")
    p.add_argument("--angle", metavar="DEG",
                   help="field in the A-B plane at DEG degrees from A, or 'diag' for the A-B diagonal; overrides --axis")
    p.add_argument("--N", type=float, help="use this demagnetizing factor instead of the prism formula")
    p.add_argument("--volume", type=float, help="magnetic volume in cm^3 (default: from --dims)")
    p.add_argument("--field-col", default="Magnetic Field", help="substring of the field column name (Oe)")
    p.add_argument("--moment-col", default="Moment", help="substring of the moment column name (emu)")
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
    elif args.angle is not None:
        N = N_in_plane(args.dims, args.angle)
    else:
        N = N_along(args.dims, "ABC".index(args.axis))
    V = args.volume or (args.dims[0] * args.dims[1] * args.dims[2] / 1000 if args.dims else None)
    if not V or V <= 0:
        p.error("give a positive --volume or --dims")
    if not 0 <= N <= 1:
        p.error(f"N must lie in [0, 1], got {N}")

    try:
        fcol, mcol, rows = read_loop(args.input, args.field_col, args.moment_col)
    except ValueError as e:
        sys.exit(str(e))
    print(f"field: '{fcol}'  moment: '{mcol}'  N = {N:.6f}  V = {V:.6g} cm^3  rows = {len(rows)}", file=sys.stderr)
    write_csv(args.output or sys.stdout, correct(rows, N, V))


if __name__ == "__main__":
    main()
