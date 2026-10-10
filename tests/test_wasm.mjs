// Checks the wasm builds against the native library (make wasm-test).
// Usage: node tests/test_wasm.mjs <build dir> <data dir>
import { readFileSync } from 'fs';
import { execFileSync } from 'child_process';
import { createEphemeris, Body, EclipseType, HouseSystem } from '../bindings/js/ephemeris.mjs';

const [build, data] = process.argv.slice(2);
let fails = 0;
const check = (ok, what) => { if (!ok) { console.log('FAIL:', what); fails++; } };

// Native reference: apparent positions at J2000 TT with eph-1as.bin.
const native = execFileSync(`${build}/test_wasm_ref`, [`${data}/eph-1as.bin`])
  .toString().trim().split('\n').map((l) => l.split(' ').map(Number));

const compare = (eph, label) => {
  for (const [body, lon, lat, speed] of native) {
    const p = eph.positionTT(body, 2451545.0);
    check(p !== null, `${label}: body ${body} available`);
    if (!p) continue;
    const err = Math.max(Math.abs(p.longitude - lon), Math.abs(p.latitude - lat),
                         Math.abs(p.longitudeSpeed - speed));
    check(err < 1e-9, `${label}: body ${body} matches native (${err})`);
  }
};

// Core: nothing without data, then load the 1" file from bytes.
const core = await createEphemeris(readFileSync(`${build}/wasm/ephemeris.wasm`));
check(core.positionTT(Body.Sun, 2451545.0) === null, 'core without data');
check(Math.abs(core.julianDate(2000, 1, 1.5) - 2451545.0) < 1e-9, 'julian date');
check(core.houses(2451545.0, 10, 50) !== null, 'houses without data');
{
  const c = core.houses(2451545.0, 10, 50, HouseSystem.Placidus);
  check(c !== null && Math.abs(core.housePosition(c, c[4] + 1e-7) - 5) < 1e-6,
        'Placidus cusps and house position');
}
check(core.loadEmbedded() === false, 'core has no built-in data');
check(core.loadData(readFileSync(`${data}/eph-1as.bin`)), 'loadData');
compare(core, 'core + eph-1as.bin');
check(!core.loadData(new Uint8Array(100)), 'bad data rejected');

// Built-in 1" data: works without loading anything.
const builtIn = await createEphemeris(readFileSync(`${build}/wasm/ephemeris-1as.wasm`));
compare(builtIn, 'ephemeris-1as.wasm');
builtIn.unloadData();
check(builtIn.positionTT(Body.Sun, 2451545.0) === null, 'unload');
check(builtIn.loadEmbedded(), 'reload built-in data');
check(builtIn.positionTT(Body.Moon, 2451545.0) !== null, 'after reload');

// Built-in 1" data with and without Chiron.
check(builtIn.positionTT(Body.Chiron, 2451545.0) === null, 'ephemeris-1as.wasm has no Chiron');
const withChiron = await createEphemeris(readFileSync(`${build}/wasm/ephemeris-1as-chiron.wasm`));
check(withChiron.positionTT(Body.Chiron, 2451545.0) !== null, 'ephemeris-1as-chiron.wasm: Chiron');
compare(withChiron, 'ephemeris-1as-chiron.wasm');

// Eclipses: 2024 Apr 8 (total, saros 139, gamma 0.3431 in NASA's canon),
// 2025 Mar 14 (total lunar), and the 2024 Apr 8 eclipse seen from Mazatlan.
{
  const s = core.solarEclipse(2460400.0);
  check(s && s.type === EclipseType.Total && s.saros === 139 &&
        Math.abs(s.gamma - 0.3431) < 2e-4 && Math.abs(s.maximum - 2460409.262836) < 2e-5,
        `solar eclipse ${JSON.stringify(s)}`);
  const l = core.lunarEclipse(2460740.0);
  check(l && l.type === EclipseType.Total && l.saros === 123 &&
        Math.abs(l.umbralMagnitude - 1.1784) < 1e-3, `lunar eclipse ${JSON.stringify(l)}`);
  const where = core.solarEclipseWhere(s.maximum);
  check(where && Math.abs(where.latitude - s.latitude) < 1e-6, 'solarEclipseWhere');
  const how = core.solarEclipseHow(s.maximum, s.longitude, s.latitude);
  check(how && how.type === EclipseType.Total && Math.abs(how.magnitude - s.magnitude) < 1e-6,
        `solarEclipseHow ${JSON.stringify(how)}`);
  const loc = core.solarEclipseLocal(2460400.0, -106.42, 23.22);
  check(loc && loc.type === EclipseType.Total && loc.secondContact > 0 &&
        loc.altitude > 60, `local eclipse ${JSON.stringify(loc)}`);
}

// The lite core: everything but eclipses and stars.
{
  const lite = await createEphemeris(readFileSync(`${build}/wasm/ephemeris-lite.wasm`));
  check(lite.loadData(readFileSync(`${data}/eph-1as.bin`)), 'lite: loadData');
  compare(lite, 'ephemeris-lite.wasm + eph-1as.bin');
  check(lite.houses(2451545.0, 10, 50, HouseSystem.Placidus) !== null, 'lite: houses');
  check(lite.solarEclipse(2460400.0) === null && lite.lunarEclipse(2460400.0) === null,
        'lite: no eclipses');
  lite.loadData(readFileSync(`${data}/eph-stars.bin`));
  check(lite.starCount() === 0, 'lite: no stars listed');
  check(lite.starFind('Regulus') === -1, 'lite: no stars');
}

// The utility and listing wrappers, against known values.
{
  const near = (a, b, tol) => Math.abs(a - b) <= tol;
  check(near(core.ttToUt(core.utToTt(2460000.5, 2023), 2023), 2460000.5, 1e-9), 'ttToUt');
  const cd = core.calendarDate(2451545.25);
  check(cd.year === 2000 && cd.month === 1 && near(cd.day, 1.75, 1e-9), `calendarDate ${JSON.stringify(cd)}`);
  const u = core.utcToJd(2000, 1, 1, 12, 0, 0);
  check(u && near(u.jdTt, 2451545.0 + 64.184 / 86400, 1e-9), `utcToJd ${JSON.stringify(u)}`);
  check(core.utcToJd(2000, 13, 1) === null, 'utcToJd rejects month 13');
  check(Math.abs(core.tdbMinusTt(2451545.0)) < 0.002, 'tdbMinusTt');
  check(near(core.siderealTime(2451545.0), 280.46, 0.01), 'siderealTime');
  const eot = core.equationOfTime(core.julianDate(2026, 11, 3.5));
  check(eot !== null && near(eot, 16.4, 0.2), `equationOfTime ${eot}`);
  check(near(core.meanObliquity(2451545.0) * 3600, 84381.406, 1e-6), 'meanObliquity (IAU 2006)');
  const nut = core.nutation(2451545.0);
  check(near(nut.longitude * 3600, -13.93, 0.05) && near(nut.obliquity * 3600, -5.77, 0.05),
        `nutation ${JSON.stringify(nut)}`);
  const eq = core.eclipticToEquatorial(90, 0, 2451545.0);
  check(near(eq.rightAscension, 90, 1e-9) &&
        near(eq.declination, core.meanObliquity(2451545.0) + nut.obliquity, 1e-9),
        `eclipticToEquatorial ${JSON.stringify(eq)}`);
  const h = core.heliocentric(Body.Mars, 2451545.0);
  check(h && near(h.distance, 1.39, 0.01) && core.heliocentric(Body.Moon, 2451545.0) === null,
        `heliocentric ${JSON.stringify(h)}`);
  check(core.bodyName(Body.Moon) === 'Moon' && core.ayanamsaName(0) === 'Lahiri', 'names');
  check(near(core.midpoint(350, 10), 0, 1e-9) || near(core.midpoint(350, 10), 360, 1e-9), 'midpoint');
  const sp = core.splitLongitude(123.5);
  check(sp.sign === 4 && sp.degrees === 3 && sp.minutes === 30 && near(sp.seconds, 0, 1e-6),
        `splitLongitude ${JSON.stringify(sp)}`);
  const polar = core.housesInfo(2460000.5, 25, 70, HouseSystem.Placidus);
  const temperate = core.housesInfo(2460000.5, 25, 50, HouseSystem.Placidus);
  check(polar?.fallback === true && temperate?.fallback === false &&
        temperate.cusps.every((c, i) => c === core.houses(2460000.5, 25, 50, HouseSystem.Placidus)[i]),
        'housesInfo fallback flag');
  check(core.heliocentricCrossing(Body.Mars, 0, 2451545.0) > 2451545.0, 'heliocentricCrossing');
}

// The table and asteroids load into the core as well.
check(core.loadData(readFileSync(`${data}/eph-asteroids.bin`)), 'asteroids');
check(core.positionTT(Body.Ceres, 2451545.0) !== null, 'Ceres');
check(core.loadData(readFileSync(`${data}/de440/de440-2000.bin`)), 'DE440 table');

if (fails === 0) console.log('All wasm checks passed.');
process.exit(fails === 0 ? 0 : 1);
