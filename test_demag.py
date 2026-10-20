"""Run: python3 test_demag.py"""
import csv
import os
import tempfile
from math import isclose, pi

from demag import N_along, N_in_plane, diagonal_angle, main, prism_N

# cube, and the three factors of any prism sum to 1
assert isclose(prism_N(1, 1, 1), 1 / 3, abs_tol=1e-12)
for a, b, c in [(2, 3, 5), (5, 5, 0.5), (10, 4, 1e-6)]:
    assert isclose(prism_N(b, c, a) + prism_N(c, a, b) + prism_N(a, b, c), 1, abs_tol=1e-7), (a, b, c)  # float cancellation ~1e-9 at aspect 1e7

# limits: thin film out of plane -> 1, long bar along its axis -> 0, across it -> 1/2
assert isclose(prism_N(5, 5, 50e-6), 1, abs_tol=1e-4)  # 5x5 mm film, 50 nm thick
assert prism_N(1, 1, 1e4) < 1e-3
assert isclose(prism_N(1, 1e4, 1), 0.5, abs_tol=1e-3)

# field along any side: the three directions sum to 1, a film is ~0 in plane
dims = (5, 3, 0.1)
assert isclose(sum(N_along(dims, i) for i in range(3)), 1, abs_tol=1e-12)
assert N_along(dims, 0) < N_along(dims, 1) < N_along(dims, 2)

# in-plane angle: 0 and 90 deg reduce to A and B; a square is the same at any angle
assert isclose(N_in_plane(dims, 0), N_along(dims, 0)) and isclose(N_in_plane(dims, 90), N_along(dims, 1))
assert isclose(N_in_plane((4, 4, 1), 30), N_along((4, 4, 1), 0))
assert isclose(diagonal_angle((4, 4, 1)), 45) and isclose(diagonal_angle((3, 4, 1)), 53.13010235415598)

# end to end on a fake MPMS .dat: H_int = H - 4 pi N m / V
with tempfile.TemporaryDirectory() as d:
    src, dst = os.path.join(d, "loop.dat"), os.path.join(d, "out.csv")
    with open(src, "w") as f:
        f.write("[Header]\nTITLE,x\n[Data]\nTime,Magnetic Field (Oe),DC Moment Free Ctr (emu)\n"
                "0,1000,0.01\n1,500,\n2,-1000,-0.01\n")
    main([src, "-o", dst, "--dims", "2", "2", "2"])
    rows = list(csv.reader(open(dst)))[1:]
    assert len(rows) == 2  # blank moment row skipped
    H, m, M, Hi = map(float, rows[0])
    assert isclose(M, 0.01 / 0.008) and isclose(Hi, 1000 - 4 * pi / 3 * M)
    main([src, "-o", dst, "--dims", "4", "1", "2", "--axis", "A"])
    H, m, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(Hi, 1000 - 4 * pi * prism_N(1, 2, 4) * M)
    main([src, "-o", dst, "--dims", "3", "4", "1", "--angle", "diag"])
    H, m, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(Hi, 1000 - 4 * pi * (0.36 * N_along((3, 4, 1), 0) + 0.64 * N_along((3, 4, 1), 1)) * M)

print("all checks passed")
