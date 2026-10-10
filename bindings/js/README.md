# @astrologer/ephemera

A pocket-sized, MIT-licensed ephemeris for astrology and astronomy, as a
WebAssembly module with no dependencies and no imports.  It works in
Node (18+), browsers and bundlers.

* The Sun, Moon, planets, Pluto, Chiron, Ceres, Pallas, Juno and Vesta,
  plus the mean and true lunar nodes and mean and osculating Lilith.
* Planetary nodes and apsides, mean and osculating (`orbit()`).
* The Uranian points: Cupido, Hades, Zeus, Kronos, Apollon, Admetos,
  Vulkanus and Poseidon.
* Houses (Placidus, Koch, Whole Sign, Equal, Porphyry, Regiomontanus,
  Campanus, Alcabitius, Topocentric and more) and the angles: Vertex and
  co-ascendants.
* 42 ayanamsas: Lahiri, Fagan/Bradley, Krishnamurti, Raman and others,
  including the star-based ones.
* 410 named fixed stars.
* Solar and lunar eclipses, both global and local circumstances.
* Phases and magnitudes, ingresses and returns, and stations.
* Accuracy: about 0.5" with the built-in data over 1550–2650, matched
  against JPL DE440.

```js
import { load, Body, HouseSystem, Ayanamsa } from '@astrologer/ephemera';

const eph = await load();                          // 1" data built in
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
`midpoint` and `splitLongitude`.  `index.d.ts` documents them all.

## Modules and data

`load({ module, data })` picks a WebAssembly module and adds data files
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
const eph = await load({ module: 'core', data: ['01as', 'stars', tableBytes] });
```

In a browser without a bundler, serve the package directory and
`load()` fetches the files next to `index.mjs`.  `demo.html` is a small
example page.

## More

The C library, the accuracy figures, the methods and the references are
in the main repository: https://github.com/intellecat/ephemera.
License: MIT.
