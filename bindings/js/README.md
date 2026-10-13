<p align="center"><img alt="Ephemera" src="https://raw.githubusercontent.com/intellecat/ephemera/main/assets/logo/lockup-horizontal-light.svg" width="420"></p>

# @astrologer/ephemera

**A pocket-sized, MIT-licensed ephemeris for astrology and astronomy.**
JPL-grade planetary positions, houses, ayanamsas, eclipses, fixed stars and
1,000+ asteroids, in a 49 KB WebAssembly core with no dependencies. It
works in Node 18+, browsers and bundlers.

- **Tiny:** 49 KB core plus 145 KB of planetary data, built in by default.
- **Accurate:** within 0.6″ of NASA JPL's DE440 over 1550–2650, or 0.03″
  with the high-precision tables.
- **MIT licensed:** use it in closed-source and commercial apps.
- **Complete for astrology:** 14 house systems, 42 ayanamsas, nodes and
  Lilith, Uranian points, eclipses, returns and stations, sidereal
  returns, planetary nodes and apsides, and 411 fixed stars.
- **Independently verified** against JPL, NASA's eclipse canon and
  India's official ephemeris. See the
  [accuracy page](https://github.com/intellecat/ephemera/blob/main/docs/accuracy.md).

```js
import { createEphemera, Body, HouseSystem, Ayanamsa } from '@astrologer/ephemera';

const eph = await createEphemera();                 // 1" data built in
const jd = eph.julianDate(2026, 1, 1.5);           // UT
const moon = eph.position(Body.Moon, jd);          // degrees, AU, per day
const cusps = eph.houses(jd, -0.1278, 51.5074, HouseSystem.Placidus);
const lahiri = eph.ayanamsa(Ayanamsa.Lahiri, eph.utToTt(jd, 2026));
```

Angles are in degrees and times are Julian days, UT or TT as each method's
name says.  Besides positions, houses (`housesInfo` also reports the polar
fallback to Porphyry), ayanamsas, searches, eclipses and stars (`starCount`,
`starName`, ...), there are the time and frame helpers: `utcToJd`,
`calendarDate`, `ttToUt`, `siderealTime`, `equationOfTime`, `nutation`,
`eclipticToEquatorial` (declination for parallels), `heliocentric`,
`midpoint` and `splitLongitude`.  For sidereal returns and ingresses pass
`{ ayanamsa: Ayanamsa.Lahiri }` (or a custom `{ jdT0, valueAtT0 }`) to
`longitudeCrossing` or `heliocentricCrossing`.  `index.d.ts` documents them all.

## Modules and data

`createEphemera({ module, data })` picks a WebAssembly module and adds data files
from the package:

| `module` | Size (brotli) | Data built in |
|---|---|---|
| `'1as'` (default) | 162 KB | 1" series, 1550–2650 |
| `'1as-chiron'` | 217 KB | the same plus Chiron |
| `'core'` | 49 KB | none: pass `data` |
| `'lite'` | 42 KB | none; no eclipses or fixed stars |

The `data` names are `'1as'`, `'01as'` (~0.2"), `'chiron'`,
`'1as-chiron'`, `'01as-chiron'`, `'asteroids'` (Ceres–Vesta and Chiron)
and `'stars'`.  Raw bytes work too, e.g. a DE440 table (~0.002") from the
repository's `data/de440/`:

```js
const eph = await createEphemera({ module: 'core', data: ['01as', 'stars', tableBytes] });
```

### Minor planets

More asteroids, centaurs and TNOs (1900–2100, ~0.5″) come in two optional
files from the repository, not in this package:

- `eph-minor.bin`, 180 popular bodies, 5.7 MB;
- `eph-minor-1000.bin`, numbered 1–1000, 17 MB.

```js
const url = 'https://raw.githubusercontent.com/intellecat/ephemera/main/data/eph-minor.bin';
const minor = new Uint8Array(await (await fetch(url)).arrayBuffer());
const eph = await createEphemera({ data: [minor] });
const eros = eph.minorBody('Eros');          // or '433', '433 Eros'
const pos = eph.position(eros, jd);
```

In a browser without a bundler, serve the package directory and
`createEphemera()` fetches the files next to `index.mjs`.  `demo.html` is a small
example page.

## More

The C library, the accuracy figures, the methods and the references are
in the main repository: https://github.com/intellecat/ephemera.
License: MIT.
