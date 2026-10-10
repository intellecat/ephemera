// Types for @astrologer/ephemera.  Angles in degrees, distances in AU,
// times as Julian days (UT or TT as named).

export declare const Body: {
  readonly Mercury: 0;
  readonly Venus: 1;
  readonly Earth: 2;
  readonly Mars: 3;
  readonly Jupiter: 4;
  readonly Saturn: 5;
  readonly Uranus: 6;
  readonly Neptune: 7;
  readonly Pluto: 8;
  readonly Sun: 9;
  readonly Moon: 10;
  readonly Ceres: 11;
  readonly Pallas: 12;
  readonly Juno: 13;
  readonly Vesta: 14;
  readonly Chiron: 15;
  readonly MeanNode: 16;
  readonly TrueNode: 17;
  readonly MeanLilith: 18;
  readonly OscuLilith: 19;
};
export type Body = (typeof Body)[keyof typeof Body];

export declare const Ayanamsa: {
  readonly Lahiri: 0;
  readonly FaganBradley: 1;
  readonly Raman: 2;
  readonly Krishnamurti: 3;
  readonly Yukteshwar: 4;
  readonly DeLuce: 5;
  readonly UshaShashi: 6;
  readonly DjwhalKhul: 7;
  readonly JnBhasin: 8;
  readonly BabylonianKugler1: 9;
  readonly BabylonianKugler2: 10;
  readonly BabylonianKugler3: 11;
  readonly BabylonianHuber: 12;
  readonly BabylonianEtaPiscium: 13;
  readonly BabylonianAldebaran: 14;
  readonly Hipparchos: 15;
  readonly Sassanian: 16;
  readonly J2000: 17;
  readonly J1900: 18;
  readonly B1950: 19;
  readonly Suryasiddhanta: 20;
  readonly SuryasiddhantaMeanSun: 21;
  readonly Aryabhata: 22;
  readonly AryabhataMeanSun: 23;
  readonly SsRevati: 24;
  readonly SsCitra: 25;
  readonly Skydram: 26;
  readonly Aryabhata522: 27;
  readonly BabylonianBritton: 28;
  readonly GalacticEquatorFiorenza: 29;
  readonly VettiusValens: 30;
  readonly Lahiri1940: 31;
  readonly LahiriVp285: 32;
  readonly KrishnamurtiSenthilathiban: 33;
  readonly LahiriIcrc: 34;
  readonly TrueCitra: 35;
  readonly TrueRevati: 36;
  readonly TruePushya: 37;
  readonly TrueMula: 38;
  readonly GalacticCenter: 39;
  readonly GalacticCenterBrand: 40;
  readonly GalacticCenterCochrane: 41;
};
export type Ayanamsa = (typeof Ayanamsa)[keyof typeof Ayanamsa];

export declare const HouseSystem: {
  readonly WholeSign: 0;
  readonly Equal: 1;
  readonly Porphyry: 2;
  readonly Placidus: 3;
  readonly Koch: 4;
  readonly Regiomontanus: 5;
  readonly Campanus: 6;
  readonly Alcabitius: 7;
  readonly Topocentric: 8;
  readonly Vehlow: 9;
  readonly EqualMC: 10;
  readonly EqualAries: 11;
  readonly Meridian: 12;
  readonly Morinus: 13;
};
export type HouseSystem = (typeof HouseSystem)[keyof typeof HouseSystem];

export declare const EclipseType: {
  readonly Total: 1;
  readonly Annular: 2;
  readonly Hybrid: 4;
  readonly Partial: 8;
  readonly Penumbral: 16;
  readonly NonCentral: 32;
};
export type EclipseType = (typeof EclipseType)[keyof typeof EclipseType];

export interface Position {
  longitude: number;
  latitude: number;
  distance: number;
  longitudeSpeed: number;
  latitudeSpeed: number;
  distanceSpeed: number;
}

export interface HouseAngles {
  ascendant: number;
  mc: number;
  armc: number;
  vertex: number;
  equatorialAscendant: number;
  coAscendantKoch: number;
  coAscendantMunkasey: number;
  polarAscendant: number;
}

export interface Phenomena {
  phaseAngle: number;
  /** Illuminated fraction, 0..1. */
  phase: number;
  elongation: number;
  diameter: number;
  magnitude: number;
}

/** Times are TT Julian days; 0 marks a phase that does not occur. */
export interface SolarEclipse {
  /** EclipseType flags. */
  type: number;
  maximum: number;
  gamma: number;
  magnitude: number;
  latitude: number;
  longitude: number;
  saros: number;
  partialBegin: number;
  partialEnd: number;
  centralBegin: number;
  centralEnd: number;
  centerLineBegin: number;
  centerLineEnd: number;
}

export interface LocalSolarEclipse {
  type: number;
  maximum: number;
  firstContact: number;
  secondContact: number;
  thirdContact: number;
  fourthContact: number;
  magnitude: number;
  obscuration: number;
  /** Of the Sun at maximum; below 0: not visible. */
  altitude: number;
  saros: number;
}

export interface LunarEclipse {
  type: number;
  maximum: number;
  gamma: number;
  umbralMagnitude: number;
  penumbralMagnitude: number;
  saros: number;
  penumbralBegin: number;
  penumbralEnd: number;
  partialBegin: number;
  partialEnd: number;
  totalBegin: number;
  totalEnd: number;
}

export interface Ephemeris {
  /** Load a data file's bytes; false if rejected. */
  loadData(bytes: Uint8Array): boolean;
  unloadData(): void;
  /** Reload the built-in data (false for modules without any). */
  loadEmbedded(): boolean;

  julianDate(year: number, month: number, day: number): number;
  deltaT(year: number): number;
  utToTt(jdUt: number, year: number): number;

  /** Geocentric position, apparent by default; null without data. */
  position(body: Body, jdUt: number, options?: { apparent?: boolean; year?: number }): Position | null;
  positionTT(body: Body, jdTt: number, options?: { apparent?: boolean }): Position | null;
  topocentric(body: Body, jdUt: number, lon: number, lat: number, elevationM?: number,
              options?: { year?: number }): Position | null;

  ascMc(jdUt: number, lon: number, lat: number): { asc: number; mc: number };
  /** Twelve cusps; null on error. */
  houses(jdUt: number, lon: number, lat: number, system?: HouseSystem): number[] | null;
  houseAngles(jdUt: number, lon: number, lat: number): HouseAngles;
  /** House 1.0 .. 12.999 of a longitude, given cusps from houses(). */
  housePosition(cusps: number[], lon: number): number;

  phenomena(body: Body, jdTt: number): Phenomena | null;

  /** UT in and out; null if not found within maxDays. */
  longitudeCrossing(body: Body, lon: number, jdUt: number, direction?: 1 | -1,
                    maxDays?: number): number | null;
  /** TT in and out. */
  nodeCrossing(body: Body, jdTt: number, direction?: 1 | -1, maxDays?: number):
    { jd: number; node: 'ascending' | 'descending'; longitude: number } | null;
  /** TT in and out. */
  station(body: Body, jdTt: number, direction?: 1 | -1, maxDays?: number):
    { jd: number; turns: 'retrograde' | 'direct' } | null;

  solarEclipse(jdTt: number, direction?: 1 | -1): SolarEclipse | null;
  solarEclipseLocal(jdTt: number, lon: number, lat: number, elevationM?: number,
                    direction?: 1 | -1): LocalSolarEclipse | null;
  lunarEclipse(jdTt: number, direction?: 1 | -1): LunarEclipse | null;

  /** Star index by name, "alf Leo" or "HIP 49669"; -1 if unknown. */
  starFind(name: string): number;
  starPosition(index: number, jdUt: number, options?: { year?: number }): Position | null;

  ayanamsa(system: Ayanamsa, jdTt: number): number;
  ayanamsaCustom(jdT0: number, valueAtT0: number, jdTt: number): number;
  ayanamsaLahiri(jdTt: number): number;
}

/** Instantiate a module from its .wasm bytes (or a compiled module). */
export declare function createEphemeris(
  source: BufferSource | WebAssembly.Module): Promise<Ephemeris>;

export type ModuleName = '1as' | '1as-chiron' | 'core' | 'lite';
/** Packaged data files (data/eph-<name>.bin). */
export type DataName =
  | '1as' | '1as-chiron' | '01as' | '01as-chiron' | 'chiron' | 'asteroids' | 'stars';

/** Load a packaged module and data files. */
export declare function load(options?: {
  module?: ModuleName;
  data?: (DataName | Uint8Array)[];
}): Promise<Ephemeris>;
