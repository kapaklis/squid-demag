"""Run: python3 test_demag.py"""
import csv
import os
import tempfile
from math import isclose, log, pi, sqrt, tanh

from demag import pick_moment, read_info, read_loop, unit_warning, cylinder_N, ellipsoid_N, factors, volume, N_along, N_in_plane, diagonal_angle, fit_background, main, prism_N, reference_background

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

# ellipsoid: sphere 1/3, trace 1, prolate spheroid against its closed form
assert isclose(ellipsoid_N(2, 2, 2), 1 / 3)
assert isclose(sum(factors((5, 3, 0.2), "ellipsoid")), 1)
r = 4  # long axis / short axis
assert isclose(ellipsoid_N(r, 1, 1), (r / sqrt(r * r - 1) * log(r + sqrt(r * r - 1)) - 1) / (r * r - 1))

# cylinder: L = D gives 0.3116 (Chen, Brug & Goldfarb 1991), thin disc -> 1, long rod -> 0, transverse = (1-N)/2
assert isclose(cylinder_N(1, 1), 0.3116, abs_tol=1e-4)
for tau in (1e-3, 1e-5, 1e-7, 1e-9):  # thin disc: 1 - N -> (2 tau / pi)(ln(4 / tau) - 1/2)
    assert isclose(cylinder_N(tau, 1), 1 - 2 * tau / pi * (log(4 / tau) - 0.5), abs_tol=max(1e-3 * tau, 1e-7)), tau
for tau in (1e2, 1e4, 1e6):  # long rod: N -> 4 / (3 pi tau)
    assert isclose(cylinder_N(tau, 1), 4 / (3 * pi * tau), rel_tol=2 / tau), tau  # next term is -1/(2 tau^2)
NA, NB, NC = factors((2, 2, 3), "cylinder")
assert NA == NB and isclose(NA + NB + NC, 1)
try:
    factors((2, 3, 3), "cylinder")
    raise AssertionError("expected A != B to be refused")
except ValueError:
    pass
assert isclose(volume((2, 2, 1), "cylinder"), pi / 1000) and isclose(volume((2, 2, 2), "ellipsoid"), pi / 6 * 8 / 1000)

# background fit recovers a diamagnetic slope under a saturating loop, on both branches or one side only
chi = -2e-7
loop = [(H, s * 0.05 * tanh((H + s * 400) / 300) + chi * H) for s in (1, -1) for H in range(-50000, 50001, 1000)]
assert isclose(fit_background(loop, 20000), chi, rel_tol=1e-6)
assert isclose(fit_background([r for r in loop if r[0] > 0], 20000), chi, rel_tol=1e-6)
try:
    fit_background(loop, 1e6)
    raise AssertionError("expected too few points")
except ValueError:
    pass

# reference substrate: branches averaged, interpolated, scaled, refuses fields well outside its range
bg = reference_background([(-1000, 2e-4), (0, 0.0), (1000, -2e-4), (1000, -4e-4)], scale=2)
assert isclose(bg(500), 2 * -1.5e-4) and isclose(bg(1005), 2 * -3e-4)
try:
    bg(2000)
    raise AssertionError("expected out of range")
except ValueError:
    pass

# end to end on a fake MPMS .dat: H_int = H - 4 pi N m / V
with tempfile.TemporaryDirectory() as d:
    src, dst = os.path.join(d, "loop.dat"), os.path.join(d, "out.csv")
    with open(src, "w") as f:
        f.write("[Header]\nTITLE,x\n[Data]\nTime,Magnetic Field (Oe),DC Moment Free Ctr (emu)\n"
                "0,1000,0.01\n1,500,\n2,-1000,-0.01\n")
    main([src, "-o", dst, "--dims", "2", "2", "2"])
    rows = list(csv.reader(open(dst)))[1:]
    assert len(rows) == 2  # blank moment row skipped
    H, m, b, M, Hi = map(float, rows[0])
    assert isclose(M, 0.01 / 0.008) and isclose(Hi, 1000 - 4 * pi / 3 * M)
    main([src, "-o", dst, "--dims", "4", "1", "2", "--axis", "A"])
    H, m, b, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(Hi, 1000 - 4 * pi * prism_N(1, 2, 4) * M)
    main([src, "-o", dst, "--dims", "3", "4", "1", "--angle", "diag"])
    H, m, b, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(Hi, 1000 - 4 * pi * (0.36 * N_along((3, 4, 1), 0) + 0.64 * N_along((3, 4, 1), 1)) * M)
    main([src, "-o", dst, "--dims", "2", "2", "2", "--shape", "ellipsoid"])
    H, m, b, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(M, 0.01 / (pi / 6 * 0.008)) and isclose(Hi, 1000 - 4 * pi / 3 * M)
    main([src, "-o", dst, "--dims", "2", "2", "2", "--bg-slope=-1e-6"])
    H, m, b, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(b, -1e-3) and isclose(M, (0.01 + 1e-3) / 0.008) and isclose(Hi, 1000 - 4 * pi / 3 * M)
    main([src, "-o", dst, "--dims", "2", "2", "2", "--bg-file", src, "--bg-scale", "0.5"])
    H, m, b, M, Hi = map(float, list(csv.reader(open(dst)))[1])
    assert isclose(b, 0.005) and isclose(M, 0.005 / 0.008)

    # MPMS3 layout: INFO header (values may hold commas), several moment columns, DC scans leave 'Moment (emu)' empty
    mp3 = os.path.join(d, "dc.dat")
    with open(mp3, "w") as f:
        f.write("[Header]\n; MPMS3 Data File (default extension .dat)\n"
                "INFO,MPMS3 Measurement Release 1.2.2 Build 0, MultiVu Release 2.3.6.1,APPNAME\n"
                "INFO,4x4 mm^2,SAMPLE_SIZE\nINFO,,SAMPLE_MASS\n[Data]\n"
                "Comment,Magnetic Field (Oe),Moment (emu),M. Std. Err. (emu),AC Moment (emu),"
                "DC Moment Fixed Ctr (emu),DC Moment Err Fixed Ctr (emu),DC Moment Free Ctr (emu)\n"
                ",1000,,,,0.011,1e-6,0.012\n,-1000,,,,-0.011,1e-6,-0.012\n")
    info = read_info(mp3)
    assert info["APPNAME"] == "MPMS3 Measurement Release 1.2.2 Build 0, MultiVu Release 2.3.6.1"
    assert info["SAMPLE_SIZE"] == "4x4 mm^2" and info["SAMPLE_MASS"] == ""
    fcol, mcol, rows = read_loop(mp3, "Magnetic Field")
    assert mcol == "DC Moment Free Ctr (emu)" and rows == [(1000, 0.012), (-1000, -0.012)]
    assert read_loop(mp3, "Magnetic Field", "Fixed Ctr")[1] == "DC Moment Fixed Ctr (emu)"
    assert pick_moment(["Moment (emu)", "DC Moment Free Ctr (emu)"], [["1", ""], ["2", ""]]) == 0  # VSM file
    assert read_info(src) == {}
assert unit_warning("Magnetic Field (Oe)", "Moment (emu)") == ""
assert "Moment (A m2)" in unit_warning("Magnetic Field (Oe)", "Moment (A m2)")

print("all checks passed")
