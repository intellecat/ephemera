// Checks the wasm builds against the native library (make wasm-test).
// Usage: node tests/test_wasm.mjs <build dir> <data dir>
import { readFileSync } from 'fs';
import { execFileSync } from 'child_process';
import { createEphemeris, Body, HouseSystem } from '../examples/wasm/ephemeris.mjs';

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

// The table and asteroids load into the core as well.
check(core.loadData(readFileSync(`${data}/eph-asteroids.bin`)), 'asteroids');
check(core.positionTT(Body.Ceres, 2451545.0) !== null, 'Ceres');
check(core.loadData(readFileSync(`${data}/de440/de440-2000.bin`)), 'DE440 table');

if (fails === 0) console.log('All wasm checks passed.');
process.exit(fails === 0 ? 0 : 1);
