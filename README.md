# squid-demag

Corrects SQUID magnetometry hysteresis loops for the demagnetizing field of the specimen:

    M     = moment / V                 (emu/cm^3)
    H_int = H_applied - 4*pi*N*M       (Oe)

`N` for a rectangular prism (crystal, plate, thin film) comes from Aharoni's closed form,
J. Appl. Phys. **83**, 3432 (1998). For any other shape, pass `N` yourself.

## GUI

```bash
python3 gui.py            # or: python3 gui.py loop.dat
```

Open a data file, pick the field and moment columns, enter the prism sides A, B, C (mm) and the side the
field is applied along. The plot shows M against the applied and the internal field and updates as you type.
Fill in N or the volume to override the prism values. Save writes the corrected CSV.

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

Thin film measured out of plane (5 x 5 mm, 50 nm thick):

```bash
python3 demag.py loop.dat --dims 5 5 50e-6 -o loop_corrected.csv
```

Known N (sphere 1/3, ellipsoid, cylinder from a table), volume in cm^3:

```bash
python3 demag.py loop.dat --N 0.3333 --volume 0.0042 -o loop_corrected.csv
```

Input is an MPMS `.dat` (data after `[Data]`) or any CSV with a header row. Columns are found by
substring (`--field-col`, default "Magnetic Field"; `--moment-col`, default "Moment"); the tool prints
which columns it picked, so check that it took the moment column you want (e.g. "DC Moment Free Ctr").

Output CSV: `H_applied (Oe), moment (emu), M (emu/cm3), H_int (Oe)`.

## Check

```bash
python3 test_demag.py
```
