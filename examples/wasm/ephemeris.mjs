// ephemeris.mjs -- JavaScript wrapper for build/wasm/ephemeris*.wasm
// (make wasm).  Works in browsers and Node; the module needs no imports.
//
//   import { createEphemeris, Body } from './ephemeris.mjs';
//   const eph = await createEphemeris(await (await fetch('ephemeris.wasm')).arrayBuffer());
//   eph.loadData(new Uint8Array(await (await fetch('eph-1as.bin')).arrayBuffer()));
//   const jd = eph.julianDate(2026, 1, 1.5);
//   const moon = eph.position(Body.Moon, jd);   // degrees, AU, per day
//
// With ephemeris-1as.wasm the 1" data is built in and loadData is optional.
// Angles are in degrees here (the C API uses radians).

export const Body = Object.freeze({
  Mercury: 0, Venus: 1, Earth: 2, Mars: 3, Jupiter: 4, Saturn: 5, Uranus: 6,
  Neptune: 7, Pluto: 8, Sun: 9, Moon: 10, Ceres: 11, Pallas: 12, Juno: 13,
  Vesta: 14, Chiron: 15, MeanNode: 16, TrueNode: 17, MeanLilith: 18,
  OscuLilith: 19,
});

export const Ayanamsa = Object.freeze({
  Lahiri: 0,
  FaganBradley: 1,
  Raman: 2,
  Krishnamurti: 3,
  Yukteshwar: 4,
  DeLuce: 5,
  UshaShashi: 6,
  DjwhalKhul: 7,
  JnBhasin: 8,
  BabylonianKugler1: 9,
  BabylonianKugler2: 10,
  BabylonianKugler3: 11,
  BabylonianHuber: 12,
  BabylonianEtaPiscium: 13,
  BabylonianAldebaran: 14,
  Hipparchos: 15,
  Sassanian: 16,
  J2000: 17,
  J1900: 18,
  B1950: 19,
  Suryasiddhanta: 20,
  SuryasiddhantaMeanSun: 21,
  Aryabhata: 22,
  AryabhataMeanSun: 23,
  SsRevati: 24,
  SsCitra: 25,
  Skydram: 26,
  Aryabhata522: 27,
  BabylonianBritton: 28,
  GalacticEquatorFiorenza: 29,
  VettiusValens: 30,
  Lahiri1940: 31,
  LahiriVp285: 32,
  KrishnamurtiSenthilathiban: 33,
  LahiriIcrc: 34,
});

export const HouseSystem = Object.freeze({
  WholeSign: 0, Equal: 1, Porphyry: 2, Placidus: 3, Koch: 4, Regiomontanus: 5,
  Campanus: 6, Alcabitius: 7, Topocentric: 8, Vehlow: 9, EqualMC: 10,
  EqualAries: 11, Meridian: 12, Morinus: 13,
});

const DEG = 180 / Math.PI;

// source: an ArrayBuffer / typed array with the .wasm bytes, or a
// WebAssembly.Module.
export async function createEphemeris(source) {
  const result = await WebAssembly.instantiate(source, {});
  const instance = result.instance ?? result;
  const e = instance.exports;
  if (e._initialize) e._initialize(); // runs constructors: loads built-in data

  const scratch = e.malloc(64); // eph_position (48 bytes) / eph_observer / out doubles
  const f64 = (ptr, n) => new Float64Array(e.memory.buffer, ptr, n);

  const observer = (lonDeg, latDeg, elevationM = 0) => {
    const p = e.malloc(24);
    f64(p, 3).set([lonDeg / DEG, latDeg / DEG, elevationM]);
    return p;
  };

  const readPosition = (rc) => {
    if (rc !== 0) return null; // no data loaded for this body or date
    const [lon, lat, dist, vlon, vlat, vdist] = f64(scratch, 6);
    return {
      longitude: lon * DEG, latitude: lat * DEG, distance: dist,
      longitudeSpeed: vlon * DEG, latitudeSpeed: vlat * DEG, distanceSpeed: vdist,
    };
  };

  return {
    // Load a data file (Uint8Array); returns true on success.
    loadData(bytes) {
      const p = e.malloc(bytes.length);
      new Uint8Array(e.memory.buffer, p, bytes.length).set(bytes);
      const rc = e.eph_load_data(p, bytes.length); // copies
      e.free(p);
      return rc === 0;
    },
    unloadData: () => e.eph_unload_data(),
    loadEmbedded: () => e.eph_load_embedded() === 0,

    julianDate: (year, month, day) => e.eph_julian_date(year, month, day),
    deltaT: (year) => e.eph_delta_t(year),
    utToTt: (jdUt, year) => e.eph_ut_to_tt(jdUt, year),

    // Apparent geocentric position for a UT Julian date (null if no data).
    // year (fractional) selects Delta-T; derived from jdUt if omitted.
    position(body, jdUt, { apparent = true, year } = {}) {
      const y = year ?? 2000 + (jdUt - 2451545.0) / 365.25;
      const fn = apparent ? e.eph_position_apparent_ut : e.eph_position_ut;
      return readPosition(fn(body, jdUt, y, scratch));
    },
    positionTT(body, jdTt, { apparent = true } = {}) {
      const fn = apparent ? e.eph_position_apparent_tt : e.eph_position_tt;
      return readPosition(fn(body, jdTt, scratch));
    },
    topocentric(body, jdUt, lonDeg, latDeg, elevationM = 0, { year } = {}) {
      const y = year ?? 2000 + (jdUt - 2451545.0) / 365.25;
      const obs = observer(lonDeg, latDeg, elevationM);
      const rc = e.eph_position_topocentric_ut(body, jdUt, y, obs, scratch);
      e.free(obs);
      return readPosition(rc);
    },

    ascMc(jdUt, lonDeg, latDeg) {
      const obs = observer(lonDeg, latDeg);
      e.eph_asc_mc(jdUt, obs, scratch, scratch + 8);
      e.free(obs);
      const [asc, mc] = f64(scratch, 2);
      return { asc: asc * DEG, mc: mc * DEG };
    },
    houses(jdUt, lonDeg, latDeg, system = HouseSystem.Porphyry) {
      const obs = observer(lonDeg, latDeg);
      const cusps = e.malloc(96);
      const rc = e.eph_houses(jdUt, obs, system, cusps);
      // rc 1: the system is undefined there (polar); Porphyry cusps.
      const out = rc >= 0 ? Array.from(f64(cusps, 12), (c) => c * DEG) : null;
      e.free(cusps);
      e.free(obs);
      return out;
    },
    // Asc, MC, ARMC, Vertex, equatorial Ascendant, co-Ascendants, polar
    // Ascendant (degrees).
    houseAngles(jdUt, lonDeg, latDeg) {
      const obs = observer(lonDeg, latDeg);
      const out = e.malloc(64);
      e.eph_house_angles(jdUt, obs, out);
      const [ascendant, mc, armc, vertex, equatorialAscendant, coAscendantKoch,
             coAscendantMunkasey, polarAscendant] = Array.from(f64(out, 8), (x) => x * DEG);
      e.free(out);
      e.free(obs);
      return { ascendant, mc, armc, vertex, equatorialAscendant, coAscendantKoch,
               coAscendantMunkasey, polarAscendant };
    },
    // House (1.0 .. 12.999) of an ecliptic longitude, cusps from houses().
    housePosition(cusps, lonDeg) {
      const p = e.malloc(96);
      f64(p, 12).set(cusps.map((c) => c / DEG));
      const h = e.eph_house_position(p, lonDeg / DEG);
      e.free(p);
      return h;
    },
    // Event searches (UT in and out; null if none within maxDays).
    longitudeCrossing(body, lonDeg, jdUt, direction = 1, maxDays = 400) {
      const out = scratch;
      return e.eph_longitude_crossing_ut(body, lonDeg / DEG, jdUt, direction, maxDays, out) === 0
        ? f64(out, 1)[0] : null;
    },
    // TT in and out: { jd, node: 'ascending' | 'descending', longitude }.
    nodeCrossing(body, jdTt, direction = 1, maxDays = 40) {
      const rc = e.eph_node_crossing(body, jdTt, direction, maxDays, scratch, scratch + 8);
      if (rc < 1) return null;
      const [jd, lon] = f64(scratch, 2);
      return { jd, node: rc === 1 ? 'ascending' : 'descending', longitude: lon * DEG };
    },
    // TT in and out: { jd, turns: 'retrograde' | 'direct' }.
    station(body, jdTt, direction = 1, maxDays = 800) {
      const rc = e.eph_station(body, jdTt, direction, maxDays, scratch);
      return rc < 1 ? null : { jd: f64(scratch, 1)[0], turns: rc === 1 ? 'retrograde' : 'direct' };
    },
    ayanamsa: (system, jdTt) => e.eph_ayanamsa(system, jdTt) * DEG,
    ayanamsaCustom: (jdT0, valueDegT0, jdTt) =>
      e.eph_ayanamsa_custom(jdT0, valueDegT0 / DEG, jdTt) * DEG,
    ayanamsaLahiri: (jdTt) => e.eph_ayanamsa_lahiri(jdTt) * DEG,
  };
}
