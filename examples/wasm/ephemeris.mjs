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
  TrueCitra: 35, TrueRevati: 36, TruePushya: 37, TrueMula: 38,
  GalacticCenter: 39, GalacticCenterBrand: 40, GalacticCenterCochrane: 41,
});

export const HouseSystem = Object.freeze({
  WholeSign: 0, Equal: 1, Porphyry: 2, Placidus: 3, Koch: 4, Regiomontanus: 5,
  Campanus: 6, Alcabitius: 7, Topocentric: 8, Vehlow: 9, EqualMC: 10,
  EqualAries: 11, Meridian: 12, Morinus: 13,
});

export const EclipseType = Object.freeze({
  Total: 1, Annular: 2, Hybrid: 4, Partial: 8, Penumbral: 16, NonCentral: 32,
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
    // { phaseAngle, phase (0..1), elongation, diameter (degrees), magnitude }.
    phenomena(body, jdTt) {
      if (e.eph_phenomena_tt(body, jdTt, scratch) !== 0) return null;
      const [phaseAngle, phase, elongation, diameter, magnitude] = f64(scratch, 5);
      return { phaseAngle: phaseAngle * DEG, phase, elongation: elongation * DEG,
               diameter: diameter * DEG, magnitude };
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
    // Eclipses (TT Julian days in and out; 0 marks a phase that does not
    // occur; null if not found).  type: EclipseType flags.
    solarEclipse(jdTt, direction = 1) {
      const p = e.malloc(104);
      const rc = e.eph_solar_eclipse_next(jdTt, direction, p);
      const type = new Int32Array(e.memory.buffer, p, 1)[0];
      const saros = new Int32Array(e.memory.buffer, p + 48, 1)[0];
      const [maximum, gamma, magnitude, lat, lon] = f64(p + 8, 5);
      const [partialBegin, partialEnd, centralBegin, centralEnd, centerLineBegin,
             centerLineEnd] = f64(p + 56, 6);
      e.free(p);
      return rc !== 0 ? null : {
        type, maximum, gamma, magnitude, latitude: lat * DEG, longitude: lon * DEG, saros,
        partialBegin, partialEnd, centralBegin, centralEnd, centerLineBegin, centerLineEnd,
      };
    },
    // The next/previous solar eclipse seen from a place (check altitude:
    // the Sun may be below the horizon).
    solarEclipseLocal(jdTt, lonDeg, latDeg, elevationM = 0, direction = 1) {
      const obs = observer(lonDeg, latDeg, elevationM);
      const p = e.malloc(80);
      const rc = e.eph_solar_eclipse_local(jdTt, obs, direction, p);
      const type = new Int32Array(e.memory.buffer, p, 1)[0];
      const saros = new Int32Array(e.memory.buffer, p + 72, 1)[0];
      const [maximum, firstContact, secondContact, thirdContact, fourthContact, magnitude,
             obscuration, altitude] = f64(p + 8, 8);
      e.free(p);
      e.free(obs);
      return rc !== 0 ? null : {
        type, maximum, firstContact, secondContact, thirdContact, fourthContact,
        magnitude, obscuration, altitude: altitude * DEG, saros,
      };
    },
    lunarEclipse(jdTt, direction = 1) {
      const p = e.malloc(96);
      const rc = e.eph_lunar_eclipse_next(jdTt, direction, p);
      const type = new Int32Array(e.memory.buffer, p, 1)[0];
      const saros = new Int32Array(e.memory.buffer, p + 40, 1)[0];
      const [maximum, gamma, umbralMagnitude, penumbralMagnitude] = f64(p + 8, 4);
      const [penumbralBegin, penumbralEnd, partialBegin, partialEnd, totalBegin,
             totalEnd] = f64(p + 48, 6);
      e.free(p);
      return rc !== 0 ? null : {
        type, maximum, gamma, umbralMagnitude, penumbralMagnitude, saros,
        penumbralBegin, penumbralEnd, partialBegin, partialEnd, totalBegin, totalEnd,
      };
    },
    // Fixed stars (load eph-stars.bin): index by name, "alf Leo" or "HIP n".
    starFind(name) {
      const bytes = new TextEncoder().encode(name + '\0');
      const p = e.malloc(bytes.length);
      new Uint8Array(e.memory.buffer, p, bytes.length).set(bytes);
      const i = e.eph_star_find(p);
      e.free(p);
      return i;
    },
    starPosition(index, jdUt, { year } = {}) {
      const y = year ?? 2000 + (jdUt - 2451545.0) / 365.25;
      return readPosition(e.eph_star_position_apparent_ut(index, jdUt, y, scratch));
    },
    ayanamsa: (system, jdTt) => e.eph_ayanamsa(system, jdTt) * DEG,
    ayanamsaCustom: (jdT0, valueDegT0, jdTt) =>
      e.eph_ayanamsa_custom(jdT0, valueDegT0 / DEG, jdTt) * DEG,
    ayanamsaLahiri: (jdTt) => e.eph_ayanamsa_lahiri(jdTt) * DEG,
  };
}
