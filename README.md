# squid-demag

Corrects SQUID magnetometry hysteresis loops for a substrate background and for the demagnetizing field
of the specimen:

    M     = (moment - background(H_applied)) / V      (emu/cm^3)
    H_int = H_applied - 4*pi*N*M       (Oe)

`N` comes from the specimen shape (`--shape`, or the Prism | Ellipsoid | Cylinder selector in the GUI):

| Shape | `--dims A B C` | N from | Volume |
|---|---|---|---|
| `prism` (default): crystal, plate, film | sides | Aharoni, J. Appl. Phys. **83**, 3432 (1998) | ABC |
| `ellipsoid` | full axis lengths | Osborn, Phys. Rev. **67**, 351 (1945) | πABC/6 |
| `cylinder` | diameter, diameter, length (axis along C) | Joseph, J. Appl. Phys. **37**, 4639 (1966) | πA²C/4 |

Prism and cylinder factors are magnetometric (averaged over the body), which is what a magnetometer measuring the
total moment needs; only an ellipsoid magnetizes uniformly, so its correction is exact. For a cylinder the
transverse factor is (1 − N_axis)/2. For any other shape, pass `N` yourself.

## GUI

```bash
python3 gui.py            # or: python3 gui.py loop.dat
```

To get a double-clickable app with its icon, run `./make_app.sh`. It builds `SQUID Demag.app` here; drag it to
/Applications or the Dock. The app runs `gui.py` from this folder with the `python3` found at build time, so rebuild
it if you move the folder or change Python. `python3 make_icon.py` redraws `icon.png` and `icon.icns`.

The Geometry panel draws the prism with the field arrow (and the angle φ for in-plane fields); drag it to rotate.
Very thin sides are drawn thicker than scale so films stay visible.

Open a data file, pick the field and moment columns, enter the prism sides A, B, C (mm) and the side the
field is applied along. The plot shows M against the applied and the internal field and updates as you type.
Fill in N or the volume to override the prism values. Pick a substrate background (below); the dashed curve is the
raw loop. Save writes the corrected CSV.

## Command line

Rectangular specimen, sides in mm, field along C (the default):

```bash
python3 demag.py loop.dat --dims 5 5 0.5 -o loop_corrected.csv
```

Field along side A instead (`--axis A|B|C`):

```bash
python3 demag.py loop.dat --dims 5 5 0.5 --axis A -o loop_corrected.csv
```

Field in the A–B plane at an angle from A, or along the A–B diagonal (`--angle DEG|diag`). In the GUI,
choose "in A–B plane" and type φ, or press "diagonal":

```bash
python3 demag.py loop.dat --dims 3 4 1 --angle diag -o loop_corrected.csv
```

This uses N_eff = N_A cos²φ + N_B sin²φ and assumes M stays parallel to the field (true at saturation, and close
to it when the field is along an easy axis). For A ≠ B a transverse field (N_A − N_B) sinφ cosφ M is ignored; for a
square sample it is zero and the result is exact.

Cylinder 3 mm across and 6 mm long, field along its axis (C) or across it (A):

```bash
python3 demag.py loop.dat --shape cylinder --dims 3 3 6 -o loop_corrected.csv
python3 demag.py loop.dat --shape cylinder --dims 3 3 6 --axis A -o loop_corrected.csv
```

Thin film measured out of plane (5 x 5 mm, 50 nm thick):

```bash
python3 demag.py loop.dat --dims 5 5 50e-6 -o loop_corrected.csv
```

### Substrate background

The background is subtracted from the moment before the demagnetizing correction, so `--dims`/volume describe the
magnetic layer only. Choose one:

- `--bg-fit HMIN`: fit m = χH + M_s·sign(H) + c to the points with |H| ≥ HMIN (where the sample is saturated) and
  subtract χH. The fitted χ is printed. Only χH is removed; a vertical offset c stays in the data.
- `--bg-slope CHI`: subtract CHI·H, CHI in emu/Oe. Write negative values as `--bg-slope=-3e-8`.
- `--bg-file PATH --bg-scale S`: subtract a measured bare-substrate (or empty-holder) loop, interpolated onto the
  sample's fields and multiplied by S (e.g. the mass or area ratio of the two substrate pieces). It is read with the
  same column names as the sample, and both of its branches are averaged.

```bash
python3 demag.py film.dat --dims 5 5 50e-6 --bg-fit 20000 -o film_corrected.csv
```

Known N (sphere 1/3, ellipsoid, cylinder from a table), volume in cm^3:

```bash
python3 demag.py loop.dat --N 0.3333 --volume 0.0042 -o loop_corrected.csv
```

Input is a Quantum Design MPMS3 (or older MPMS) `.dat` file, or any CSV with a header row. In an MPMS3 file the
tool reads the data after `[Data]` and, by default, uses the moment column that actually holds data: `Moment (emu)`
for VSM measurements, `DC Moment Free Ctr (emu)` for DC scans (or name another with `--moment-col`, e.g.
`"Fixed Ctr"`). Rows without a value (failed points, comments) are skipped. The GUI shows the sample notes stored
in the file header (material, comment, size, shape, mass) to help with the dimensions. The correction assumes Oe and
emu; if the column names say otherwise, the tool warns. AC-susceptibility runs carry no DC moment and are refused.

Output CSV: `H_applied (Oe), moment (emu), background (emu), M (emu/cm3), H_int (Oe)`.

## Check

```bash
python3 test_demag.py
```

## License

MIT, see [LICENSE](LICENSE).
