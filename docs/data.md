# Data files

[README](../README.md) · [Guide](guide.md) · [Accuracy](accuracy.md) · [Data files](data.md) · [Methods](methods.md)

## Data files

| File                         | Size   | brotli | Contents                                      |
|------------------------------|--------|--------|-----------------------------------------------|
| `data/eph-1as.bin`           | 145 KB | 112 KB | series truncated to ~0.5″ (4,719 terms)        |
| `data/eph-01as.bin`          | 303 KB | 213 KB | series truncated to ~0.2″ (14,176 terms)       |
| `data/eph-full.bin`          | 649 KB | 420 KB | the full series (34,530 terms)                 |
| `data/eph-*-chiron.bin`      | +60 KB | +55 KB | each level above with Chiron (1600–2500)       |
| `data/eph-chiron.bin`        | 60 KB  | 56 KB  | Chiron alone, to add to any series file        |
| `data/eph-asteroids.bin`     | 425 KB | 397 KB | Ceres, Pallas, Juno, Vesta, Chiron, 1600–2500  |
| `data/eph-stars.bin`         | 26 KB  |        | 410 named stars and the Galactic Center        |
| `data/de440/de440-YEAR.bin`  | ~330 KB each | ~310 KB | DE440 table for YEAR..YEAR+100, ~0.002″ |

Chiron comes three ways with the same data: inside a series file
(`eph-1as-chiron.bin`), as its own file next to one (`eph-chiron.bin`),
or with the other minor bodies (`eph-asteroids.bin`).

The series files cover 1550–2650 (Pluto and the corrections) and work,
less accurately, outside it.  Load one series file, plus the asteroids if
you need them.  The series files include Pluto and the DE440 corrections;
truncation drops the VSOP87/ELP terms whose amplitude (times T^α at
1550/2650) is below a per-planet threshold chosen so the level holds over
1550–2650 (see `LEVEL_*` in the `Makefile`).

Several files can be loaded at once; a later file overrides what an
earlier one provides.  Each file is a small container (`src/data.h`):
typed sections used in place after loading, with a CRC-32, so a damaged
or mismatched file is rejected rather than giving wrong positions.

**DE440 tables.**  `data/de440/` holds one file per century (1550–1600,
1600–1700, …, 2600–2650; 3.7 MB in all): load the ones you need, side by
side.  Each stores **DE440 minus the library's own series**: the residual is small and smooth, so it packs into fixed-length
Chebyshev segments with one absolute quantization step per body and a
fixed bit width per coefficient.  The series under the Moon and the inner
planets are truncated to terms above 1e-7 rad (~0.02") and stored in the
same file (41 KB), so a table needs no other file, and the table records
a checksum of that base, so tables built against different series are
rejected; the outer planets and
the Earth-Moon barycentre are stored without a base series.  Each body is
fitted to 0.002" (seen from the Earth at its closest) and its speed to
0.01"/day; the Earth comes from the barycentre and the Moon.  For dates
no loaded table covers, a series file is used if one is loaded.
Regenerate with `python3 tools/pack_ephemeris.py` (needs numpy + skyfield
+ pyerfa, the DE440 kernel, and `make && make tools`; `--target` sets the
accuracy, `--chunk` the years per file).
The format is our own, built from public-domain JPL data, so no
third-party license applies.

## Compact encoding

The coefficient tables use a quantized layout to keep the files small
without measurable loss:

* **VSOP87D** — the secular (`C == 0`) terms are folded into an exact
  `double` polynomial (their amplitudes reach ~10⁴ rad/millennium and must
  not be float32); the small periodic terms are stored as
  `{float amplitude, uint32 phase, double frequency}` — 24 → 16 bytes/term.
* **ELP2000-82B** — the bounded phase (mod 2π) is stored in 30 bits, with
  the time power packed into the top 2 bits and only the rate kept as
  `double` — 48/24 → 40/16 bytes/term.
* **Pluto, asteroids, corrections** — float32 Chebyshev coefficients.

## Coefficient provenance and regeneration

The VSOP87D tables are downloaded from IMCCE (Observatoire de Paris) and
written to `data/vsop87d_data.cpp`; ELP2000-82B, Pluto, the DE440
corrections and the asteroid fits have their own generators in `tools/`
(the DE440 fits use the full `de440.bsp` kernel, 1550–2650, which is a
development download and never shipped).
These generated sources are not compiled into the library:
`tools/pack_series.cpp` packs them into the data files.

```sh
python3 tools/gen_vsop87.py        # downloads data/vsop87d/* and regenerates the tables
make tools && make data            # rebuilds data/eph-*.bin
```

The test reference file `tests/reference_data.h` is generated once from a
JPL DE440s kernel with Skyfield (reduced with ERFA/IAU-2006); the kernel is **not** part of the library
and is never shipped.  See `tools/make_reference.py` and
`tools/gen_reference.py` if you want to regenerate it.

## Layout

```
include/ephemeris/ephemeris.h   public C API
src/                            VSOP87, ELP2000 Moon, Pluto, asteroids,
                                DE440 corrections, data loader, DE440
                                table, frames, time, houses, searches,
                                phenomena, stars, eclipses
data/*_data.{h,cpp}             generated coefficient tables (sources of
                                the data files; only delta_t is compiled in)
data/*.bin                      data files
tools/                          coefficient, fit, data-file and reference
                                generators
tests/                          tests and their reference data (DE440,
                                ERFA, NASA, Indian Astronomical Ephemeris)
bindings/js/                    npm package @astrologer/ephemera
examples/positions.c            CLI that prints a chart
```
