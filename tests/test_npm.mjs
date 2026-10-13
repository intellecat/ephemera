// Smoke test of the staged npm package (make npm-test).
// Usage: node tests/test_npm.mjs <package dir>
import { pathToFileURL } from 'url';
import { join, resolve } from 'path';

const dir = resolve(process.argv[2]);
const { createEphemera, load, Body, HouseSystem, Ayanamsa } =
  await import(pathToFileURL(join(dir, 'index.mjs')).href);
let fails = 0;
const check = (ok, what) => { if (!ok) { console.log('FAIL:', what); fails++; } };

const eph = await createEphemera();
check(load === createEphemera, 'deprecated load alias');
const jd = eph.julianDate(2026, 1, 1.5);
check(eph.position(Body.Moon, jd) !== null, 'default module: Moon');
check(eph.houses(jd, -0.1278, 51.5074, HouseSystem.Placidus)?.length === 12, 'houses');
check(Math.abs(eph.ayanamsa(Ayanamsa.Lahiri, 2451545.0) - 23.857092325) < 1e-6, 'Lahiri');

for (const module of ['1as-chiron', 'core', 'lite']) {
  const e = await createEphemera({ module, data: module === 'core' ? ['01as', 'chiron', 'stars'] : [] });
  check(e.position(Body.Sun, jd) !== null || module === 'lite', `${module}: Sun`);
}
const core = await createEphemera({ module: 'core', data: ['01as', 'chiron', 'stars', 'asteroids'] });
check(core.position(Body.Chiron, jd) !== null && core.position(Body.Ceres, jd) !== null,
      'core + data: Chiron, Ceres');
check(core.starFind('Regulus') >= 0, 'core + stars: Regulus');
{
  const n = core.starCount();
  const reg = core.starFind('Regulus');
  let roundTrip = true;
  for (let i = 0; i < n; i++) if (core.starFind(core.starName(i)) !== i) roundTrip = false;
  check(n === 411 && roundTrip && core.starName(reg) === 'Regulus' &&
        core.starDesignation(reg) === 'alf Leo' && Math.abs(core.starMagnitude(reg) - 1.4) < 0.1,
        `star listing (${n} stars)`);
}
check(core.solarEclipse(2460400.0)?.saros === 139, 'core: eclipse');

let threw = false;
try { await createEphemera({ module: 'nope' }); } catch { threw = true; }
check(threw, 'unknown module rejected');

if (fails === 0) console.log('All npm package checks passed.');
process.exit(fails === 0 ? 0 : 1);
