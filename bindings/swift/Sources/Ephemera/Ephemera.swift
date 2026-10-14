// Ephemera for Swift: a pocket-sized, MIT-licensed ephemeris for astrology
// and astronomy.
//
//     import Ephemera
//
//     let jd = Ephemera.julianDate(2026, 1, 1.5)               // 1 Jan 2026, 12:00 UT
//     let moon = Ephemera.position(.moon, jdUT: jd)!           // degrees, AU, per day
//     let cusps = Ephemera.houses(jdUT: jd, longitude: -0.1278, latitude: 51.5074,
//                                 system: .placidus)
//
// The 1" planetary series with Chiron (1550-2650) is loaded on first use.
// Angles are degrees, distances AU, times Julian days (UT or TT as named);
// calls that find nothing return nil.  The data loaded is shared by the
// whole process: load files before computing from several threads.

import CEphemera
import Foundation

public enum Ephemera {}

// MARK: - Constants

extension Ephemera {
    /// A body: the planets, the Sun and Moon, asteroids, lunar points,
    /// Uranian points, or a loaded minor planet (`minorBody(_:)`).
    public struct Body: RawRepresentable, Hashable, Sendable {
        public let rawValue: Int32
        public init(rawValue: Int32) { self.rawValue = rawValue }

        public static let mercury = Body(rawValue: 0)
        public static let venus = Body(rawValue: 1)
        public static let earth = Body(rawValue: 2)
        public static let mars = Body(rawValue: 3)
        public static let jupiter = Body(rawValue: 4)
        public static let saturn = Body(rawValue: 5)
        public static let uranus = Body(rawValue: 6)
        public static let neptune = Body(rawValue: 7)
        public static let pluto = Body(rawValue: 8)
        public static let sun = Body(rawValue: 9)
        public static let moon = Body(rawValue: 10)
        public static let ceres = Body(rawValue: 11)
        public static let pallas = Body(rawValue: 12)
        public static let juno = Body(rawValue: 13)
        public static let vesta = Body(rawValue: 14)
        public static let chiron = Body(rawValue: 15)
        public static let meanNode = Body(rawValue: 16)
        public static let trueNode = Body(rawValue: 17)
        public static let meanLilith = Body(rawValue: 18)
        public static let oscuLilith = Body(rawValue: 19)
        // The Hamburg school's Uranian points (Witte, Sieggrün).
        public static let cupido = Body(rawValue: 20)
        public static let hades = Body(rawValue: 21)
        public static let zeus = Body(rawValue: 22)
        public static let kronos = Body(rawValue: 23)
        public static let apollon = Body(rawValue: 24)
        public static let admetos = Body(rawValue: 25)
        public static let vulkanus = Body(rawValue: 26)
        public static let poseidon = Body(rawValue: 27)

        /// Body id of the first loaded minor planet.
        public static let minorBase: Int32 = 1000
    }

    public enum Ayanamsa: Int32, CaseIterable, Sendable {
        case lahiri = 0, faganBradley, raman, krishnamurti, yukteshwar, deLuce, ushaShashi,
             djwhalKhul, jnBhasin, babylonianKugler1, babylonianKugler2, babylonianKugler3,
             babylonianHuber, babylonianEtaPiscium, babylonianAldebaran, hipparchos, sassanian,
             j2000, j1900, b1950, suryasiddhanta, suryasiddhantaMeanSun, aryabhata,
             aryabhataMeanSun, ssRevati, ssCitra, skydram, aryabhata522, babylonianBritton,
             galacticEquatorFiorenza, vettiusValens, lahiri1940, lahiriVP285,
             krishnamurtiSenthilathiban, lahiriICRC,
             // Star-based (need the stars: loadPackageData(.stars)).
             trueCitra, trueRevati, truePushya, trueMula, galacticCenter, galacticCenterBrand,
             galacticCenterCochrane
    }

    public enum HouseSystem: Int32, CaseIterable, Sendable {
        case wholeSign = 0, equal, porphyry, placidus, koch, regiomontanus, campanus, alcabitius,
             topocentric, vehlow, equalMC, equalAries, meridian, morinus
    }

    public enum OrbitMethod: Int32, Sendable {
        /// Secular mean elements (Simon et al. 1994 / Meeus ch. 31): Mercury..Neptune.
        case mean = 0
        /// From the position and velocity: the planets but the Earth, Pluto, minor bodies.
        case osculating = 1
    }

    public struct EclipseType: OptionSet, Hashable, Sendable {
        public let rawValue: Int32
        public init(rawValue: Int32) { self.rawValue = rawValue }
        public static let total = EclipseType(rawValue: 1)
        public static let annular = EclipseType(rawValue: 2)
        public static let hybrid = EclipseType(rawValue: 4)
        public static let partial = EclipseType(rawValue: 8)
        public static let penumbral = EclipseType(rawValue: 16)
        public static let nonCentral = EclipseType(rawValue: 32)
    }

    /// An ayanamsa for the sidereal searches: a system, or a custom value
    /// (degrees) at a TT Julian day, carried by the precession.
    public enum AyanamsaOption: Sendable {
        case system(Ayanamsa)
        case custom(jdT0: Double, valueAtT0: Double)
    }
}

// MARK: - Results

extension Ephemera {
    public struct Position: Hashable, Sendable {
        public let longitude, latitude, distance: Double
        public let longitudeSpeed, latitudeSpeed, distanceSpeed: Double
        public var isRetrograde: Bool { longitudeSpeed < 0 }
    }

    public struct CalendarDate: Hashable, Sendable {
        public let year, month: Int
        /// With the fraction of the day.
        public let day: Double
    }

    public struct JulianDates: Hashable, Sendable {
        public let tt, ut: Double
    }

    public struct Equatorial: Hashable, Sendable {
        public let rightAscension, declination: Double
    }

    public struct Heliocentric: Hashable, Sendable {
        public let longitude, latitude, distance: Double
    }

    public struct HousesInfo: Hashable, Sendable {
        public let cusps: [Double]
        /// The system is undefined here (polar): Porphyry was used.
        public let fallback: Bool
    }

    public struct HouseAngles: Hashable, Sendable {
        public let ascendant, mc, armc, vertex, equatorialAscendant: Double
        public let coAscendantKoch, coAscendantMunkasey, polarAscendant: Double
    }

    public struct Phenomena: Hashable, Sendable {
        public let phaseAngle: Double
        /// Illuminated fraction, 0...1.
        public let phase: Double
        public let elongation, diameter, magnitude: Double
    }

    public struct NodeCrossing: Hashable, Sendable {
        public let jd: Double
        public let ascending: Bool
        public let longitude: Double
    }

    public struct Station: Hashable, Sendable {
        public let jd: Double
        /// True when the body turns retrograde, false when it turns direct.
        public let turnsRetrograde: Bool
    }

    public struct OrbitPoints: Hashable, Sendable {
        public let ascendingNode, descendingNode, perihelion, aphelion: Position
    }

    /// A heliocentric orbit, mean ecliptic and equinox of date.
    public struct Orbit: Hashable, Sendable {
        public let semiMajorAxis, eccentricity, inclination, node, perihelion: Double
        public let heliocentric, geocentric: OrbitPoints
    }

    /// Times are TT Julian days; 0 marks a phase that does not occur.
    public struct SolarEclipse: Hashable, Sendable {
        public let type: EclipseType
        public let maximum, gamma, magnitude, latitude, longitude: Double
        public let saros: Int
        public let partialBegin, partialEnd, centralBegin, centralEnd: Double
        public let centerLineBegin, centerLineEnd: Double
    }

    public struct LocalSolarEclipse: Hashable, Sendable {
        public let type: EclipseType
        public let maximum, firstContact, secondContact, thirdContact, fourthContact: Double
        public let magnitude, obscuration: Double
        /// Of the Sun at maximum; below 0: not visible.
        public let altitude: Double
        public let saros: Int
    }

    public struct LocalEclipseState: Hashable, Sendable {
        /// Empty: no eclipse here and now.
        public let type: EclipseType
        public let magnitude, obscuration, altitude: Double
    }

    public struct LunarEclipse: Hashable, Sendable {
        public let type: EclipseType
        public let maximum, gamma, umbralMagnitude, penumbralMagnitude: Double
        public let saros: Int
        public let penumbralBegin, penumbralEnd, partialBegin, partialEnd: Double
        public let totalBegin, totalEnd: Double
    }

    public struct GeoPoint: Hashable, Sendable {
        public let latitude, longitude: Double
    }

    public struct MinorPlanet: Hashable, Sendable {
        public let body: Body
        public let number: Int
        public let name: String
    }

    public struct SplitLongitude: Hashable, Sendable {
        /// 0 = Aries ... 11 = Pisces.
        public let sign, degrees, minutes: Int
        public let seconds: Double
    }
}

// MARK: - Helpers

private let deg = 180.0 / Double.pi
private let rad = Double.pi / 180.0

/// The built-in data, loaded once on first use.
private let builtInLoaded: Bool = loadBuiltIn()

private func loadBuiltIn() -> Bool {
    guard let url = Bundle.module.url(forResource: "eph-1as-chiron", withExtension: "bin"),
          let data = try? Data(contentsOf: url) else { return false }
    return data.withUnsafeBytes { eph_load_data($0.baseAddress, $0.count) == 0 }
}

@inline(__always) private func ready() { _ = builtInLoaded }

private func year(_ jd: Double) -> Double { 2000.0 + (jd - 2451545.0) / 365.25 }

private func body(_ b: Ephemera.Body) -> eph_body { eph_body(rawValue: UInt32(bitPattern: b.rawValue)) }

private func makePosition(_ p: eph_position) -> Ephemera.Position {
    Ephemera.Position(longitude: p.longitude * deg, latitude: p.latitude * deg, distance: p.distance,
                      longitudeSpeed: p.longitude_speed * deg, latitudeSpeed: p.latitude_speed * deg,
                      distanceSpeed: p.distance_speed)
}

private func observer(_ lon: Double, _ lat: Double, _ elevation: Double = 0) -> eph_observer {
    eph_observer(longitude: lon * rad, latitude: lat * rad, elevation: elevation)
}

private func string(_ p: UnsafePointer<CChar>?) -> String { p.map { String(cString: $0) } ?? "" }

private func readPosition(_ fill: (UnsafeMutablePointer<eph_position>) -> Int32) -> Ephemera.Position? {
    var out = eph_position()
    return fill(&out) == 0 ? makePosition(out) : nil
}

// MARK: - Data files

extension Ephemera {
    /// Load a data file's bytes (copied); false if rejected.  Later files
    /// override what earlier ones provide.
    @discardableResult
    public static func loadData(_ data: Data) -> Bool {
        ready() // the built-in data first, so this file overrides it
        return data.withUnsafeBytes { eph_load_data($0.baseAddress, $0.count) == 0 }
    }

    /// Load a data file from a path or file URL.
    @discardableResult
    public static func loadData(contentsOf url: URL) -> Bool {
        ready()
        return eph_load_file(url.path) == 0
    }

    /// Unload every data file (the built-in data too: see `loadEmbedded`).
    public static func unloadData() {
        ready()
        eph_unload_data()
    }

    /// Reload the built-in data (the 1" series with Chiron).
    @discardableResult
    public static func loadEmbedded() -> Bool {
        ready()
        return loadBuiltIn()
    }
}

// MARK: - Time

extension Ephemera {
    /// Julian day of a calendar date (Gregorian; Julian before 1582-10-15);
    /// `day` carries the fraction of the day (1.5 = noon).
    public static func julianDate(_ year: Int, _ month: Int, _ day: Double) -> Double {
        eph_julian_date(Int32(year), Int32(month), day)
    }

    public static func calendarDate(_ jd: Double) -> CalendarDate {
        var y: Int32 = 0, m: Int32 = 0, d = 0.0
        eph_calendar_date(jd, &y, &m, &d)
        return CalendarDate(year: Int(y), month: Int(m), day: d)
    }

    /// TT - UT, seconds.
    public static func deltaT(_ year: Double) -> Double { eph_delta_t(year) }

    public static func utToTT(_ jdUT: Double, year: Double? = nil) -> Double {
        eph_ut_to_tt(jdUT, year ?? Self.yearOf(jdUT))
    }

    public static func ttToUT(_ jdTT: Double, year: Double? = nil) -> Double {
        eph_tt_to_ut(jdTT, year ?? Self.yearOf(jdTT))
    }

    /// Civil UTC (leap seconds from 1972) to TT and UT Julian days; nil for an invalid date.
    public static func utcToJD(_ year: Int, _ month: Int, _ day: Int, _ hour: Int = 0,
                               _ minute: Int = 0, _ second: Double = 0) -> JulianDates? {
        var tt = 0.0, ut = 0.0
        guard eph_utc_to_jd(Int32(year), Int32(month), Int32(day), Int32(hour), Int32(minute),
                            second, &tt, &ut) == 0 else { return nil }
        return JulianDates(tt: tt, ut: ut)
    }

    /// TT and UT Julian days of an instant.
    public static func julianDates(_ date: Date) -> JulianDates {
        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = TimeZone(identifier: "UTC")!
        let c = calendar.dateComponents([.year, .month, .day, .hour, .minute, .second, .nanosecond], from: date)
        let second = Double(c.second ?? 0) + Double(c.nanosecond ?? 0) / 1e9
        return utcToJD(c.year!, c.month!, c.day!, c.hour!, c.minute!, second)!
    }

    /// TDB - TT, seconds.
    public static func tdbMinusTT(_ jdTT: Double) -> Double { eph_tdb_minus_tt(jdTT) }

    /// Greenwich apparent sidereal time (IAU 2006), degrees.
    public static func siderealTime(_ jdUT: Double) -> Double { eph_sidereal_time(jdUT) * deg }

    /// Apparent minus mean solar time, minutes; nil without the Sun's data.
    public static func equationOfTime(_ jdUT: Double) -> Double? {
        ready()
        var e = 0.0
        return eph_equation_of_time(jdUT, &e) == 0 ? e * 1440 : nil
    }

    static func yearOf(_ jd: Double) -> Double { year(jd) }
}

// MARK: - Frames

extension Ephemera {
    /// Mean obliquity of the ecliptic (IAU 2006), degrees.
    public static func meanObliquity(_ jdTT: Double) -> Double { eph_mean_obliquity(jdTT) * deg }

    /// Nutation in longitude and obliquity (IAU 2000B), degrees.
    public static func nutation(_ jdTT: Double) -> (longitude: Double, obliquity: Double) {
        var dpsi = 0.0, deps = 0.0
        eph_nutation(jdTT, &dpsi, &deps)
        return (dpsi * deg, deps * deg)
    }

    /// Apparent ecliptic to equatorial coordinates (true obliquity), degrees.
    public static func eclipticToEquatorial(longitude: Double, latitude: Double, jdTT: Double) -> Equatorial {
        var ra = 0.0, dec = 0.0
        eph_ecliptic_to_equatorial(longitude * rad, latitude * rad, jdTT, &ra, &dec)
        return Equatorial(rightAscension: ra * deg, declination: dec * deg)
    }
}

// MARK: - Positions

extension Ephemera {
    /// Geocentric ecliptic position of date for a UT Julian day, apparent by
    /// default; nil without data.  `year` (fractional) selects ΔT.
    public static func position(_ b: Body, jdUT: Double, apparent: Bool = true,
                                year: Double? = nil) -> Position? {
        ready()
        let y = year ?? Self.yearOf(jdUT)
        return readPosition { apparent ? eph_position_apparent_ut(body(b), jdUT, y, $0)
                                       : eph_position_ut(body(b), jdUT, y, $0) }
    }

    public static func position(_ b: Body, jdTT: Double, apparent: Bool = true) -> Position? {
        ready()
        return readPosition { apparent ? eph_position_apparent_tt(body(b), jdTT, $0)
                                       : eph_position_tt(body(b), jdTT, $0) }
    }

    /// Heliocentric ecliptic of date; nil for the Sun, the Moon or without data.
    public static func heliocentric(_ b: Body, jdTT: Double) -> Heliocentric? {
        ready()
        var lon = 0.0, lat = 0.0, dist = 0.0
        guard eph_heliocentric_tt(body(b), jdTT, &lon, &lat, &dist) == 0 else { return nil }
        return Heliocentric(longitude: lon * deg, latitude: lat * deg, distance: dist)
    }

    /// "Sun", "Moon", ... (minor planets: their names).
    public static func name(of b: Body) -> String {
        ready()
        if b.rawValue >= Body.minorBase { return string(eph_minor_name(b.rawValue - Body.minorBase)) }
        return string(eph_body_name(body(b)))
    }

    /// Apparent position seen from a place (east longitude and latitude in
    /// degrees, elevation in metres).
    public static func topocentric(_ b: Body, jdUT: Double, longitude: Double, latitude: Double,
                                   elevation: Double = 0, year: Double? = nil) -> Position? {
        ready()
        var obs = observer(longitude, latitude, elevation)
        let y = year ?? Self.yearOf(jdUT)
        return readPosition { eph_position_topocentric_ut(body(b), jdUT, y, &obs, $0) }
    }
}

// MARK: - Houses

extension Ephemera {
    public static func ascMC(jdUT: Double, longitude: Double, latitude: Double) -> (asc: Double, mc: Double) {
        var obs = observer(longitude, latitude)
        var asc = 0.0, mc = 0.0
        eph_asc_mc(jdUT, &obs, &asc, &mc)
        return (asc * deg, mc * deg)
    }

    /// The twelve cusps, and whether the system is undefined at this latitude
    /// (Placidus, Koch, Alcabitius inside the polar circles: Porphyry is used).
    public static func housesInfo(jdUT: Double, longitude: Double, latitude: Double,
                                  system: HouseSystem = .porphyry) -> HousesInfo? {
        var obs = observer(longitude, latitude)
        var cusps = [Double](repeating: 0, count: 12)
        let rc = eph_houses(jdUT, &obs, eph_house_system(rawValue: UInt32(system.rawValue)), &cusps)
        guard rc >= 0 else { return nil }
        return HousesInfo(cusps: cusps.map { $0 * deg }, fallback: rc == 1)
    }

    /// Twelve cusps (degrees); nil on bad arguments.
    public static func houses(jdUT: Double, longitude: Double, latitude: Double,
                              system: HouseSystem = .porphyry) -> [Double]? {
        housesInfo(jdUT: jdUT, longitude: longitude, latitude: latitude, system: system)?.cusps
    }

    /// Ascendant, MC, ARMC, Vertex, equatorial Ascendant, co-Ascendants and
    /// polar Ascendant, degrees.
    public static func houseAngles(jdUT: Double, longitude: Double, latitude: Double) -> HouseAngles {
        var obs = observer(longitude, latitude)
        var a = eph_angles()
        eph_house_angles(jdUT, &obs, &a)
        return HouseAngles(ascendant: a.ascendant * deg, mc: a.mc * deg, armc: a.armc * deg,
                           vertex: a.vertex * deg, equatorialAscendant: a.equatorial_ascendant * deg,
                           coAscendantKoch: a.co_ascendant_koch * deg,
                           coAscendantMunkasey: a.co_ascendant_munkasey * deg,
                           polarAscendant: a.polar_ascendant * deg)
    }

    /// House 1.0 ..< 13.0 of a longitude, given cusps from `houses`.
    public static func housePosition(cusps: [Double], longitude: Double) -> Double {
        let c = cusps.map { $0 * rad }
        return eph_house_position(c, longitude * rad)
    }
}

// MARK: - Phenomena

extension Ephemera {
    /// Phase angle, illuminated fraction, elongation, apparent diameter
    /// (degrees) and visual magnitude; nil without data.
    public static func phenomena(_ b: Body, jdTT: Double) -> Phenomena? {
        ready()
        var p = eph_phenomena()
        guard eph_phenomena_tt(body(b), jdTT, &p) == 0 else { return nil }
        return Phenomena(phaseAngle: p.phase_angle * deg, phase: p.phase, elongation: p.elongation * deg,
                         diameter: p.diameter * deg, magnitude: p.magnitude)
    }
}

// MARK: - Searches

extension Ephemera {
    /// When the apparent longitude reaches `longitude` (returns, ingresses);
    /// UT in and out, nil if not within `maxDays`.  With `ayanamsa` the
    /// longitude is sidereal, the ayanamsa taken at the crossing itself.
    public static func longitudeCrossing(_ b: Body, longitude: Double, jdUT: Double,
                                         direction: Int = 1, maxDays: Double = 400,
                                         ayanamsa: AyanamsaOption? = nil) -> Double? {
        ready()
        var out = 0.0
        let dir = Int32(direction)
        guard let ayanamsa else {
            return eph_longitude_crossing_ut(body(b), longitude * rad, jdUT, dir, maxDays, &out) == 0 ? out : nil
        }
        let jdTT = utToTT(jdUT)
        let rc: Int32
        switch ayanamsa {
        case .system(let s):
            rc = eph_longitude_crossing_sidereal(body(b), longitude * rad,
                                                 eph_ayanamsa_system(rawValue: UInt32(s.rawValue)),
                                                 jdTT, dir, maxDays, &out)
        case .custom(let jdT0, let value):
            rc = eph_longitude_crossing_sidereal_custom(body(b), longitude * rad, jdT0, value * rad,
                                                        jdTT, dir, maxDays, &out)
        }
        return rc == 0 ? ttToUT(out) : nil
    }

    /// When the heliocentric longitude reaches `longitude`; TT in and out.
    public static func heliocentricCrossing(_ b: Body, longitude: Double, jdTT: Double,
                                            direction: Int = 1, maxDays: Double = 40000,
                                            ayanamsa: AyanamsaOption? = nil) -> Double? {
        ready()
        var out = 0.0
        let dir = Int32(direction)
        let rc: Int32
        switch ayanamsa {
        case nil:
            rc = eph_heliocentric_crossing(body(b), longitude * rad, jdTT, dir, maxDays, &out)
        case .system(let s)?:
            rc = eph_heliocentric_crossing_sidereal(body(b), longitude * rad,
                                                    eph_ayanamsa_system(rawValue: UInt32(s.rawValue)),
                                                    jdTT, dir, maxDays, &out)
        case .custom(let jdT0, let value)?:
            rc = eph_heliocentric_crossing_sidereal_custom(body(b), longitude * rad, jdT0, value * rad,
                                                           jdTT, dir, maxDays, &out)
        }
        return rc == 0 ? out : nil
    }

    /// When the body crosses the ecliptic; TT in and out.
    public static func nodeCrossing(_ b: Body, jdTT: Double, direction: Int = 1,
                                    maxDays: Double = 40) -> NodeCrossing? {
        ready()
        var jd = 0.0, lon = 0.0
        let rc = eph_node_crossing(body(b), jdTT, Int32(direction), maxDays, &jd, &lon)
        guard rc >= 1 else { return nil }
        return NodeCrossing(jd: jd, ascending: rc == 1, longitude: lon * deg)
    }

    /// When the body turns retrograde or direct; TT in and out.
    public static func station(_ b: Body, jdTT: Double, direction: Int = 1,
                               maxDays: Double = 800) -> Station? {
        ready()
        var jd = 0.0
        let rc = eph_station(body(b), jdTT, Int32(direction), maxDays, &jd)
        guard rc >= 1 else { return nil }
        return Station(jd: jd, turnsRetrograde: rc == 1)
    }
}

// MARK: - Planetary nodes and apsides

extension Ephemera {
    /// The orbit's elements (degrees, AU) and its nodes, perihelion and
    /// aphelion, heliocentric and geocentric; nil if unsupported or without data.
    public static func orbit(_ b: Body, jdTT: Double, method: OrbitMethod = .mean) -> Orbit? {
        ready()
        var o = eph_orbit()
        guard eph_orbit_tt(body(b), jdTT, eph_orbit_method(rawValue: UInt32(method.rawValue)), &o) == 0
        else { return nil }
        let h = o.heliocentric, g = o.geocentric
        return Orbit(semiMajorAxis: o.semi_major_axis, eccentricity: o.eccentricity,
                     inclination: o.inclination * deg, node: o.node * deg, perihelion: o.perihelion * deg,
                     heliocentric: OrbitPoints(ascendingNode: makePosition(h.0), descendingNode: makePosition(h.1),
                                               perihelion: makePosition(h.2), aphelion: makePosition(h.3)),
                     geocentric: OrbitPoints(ascendingNode: makePosition(g.0), descendingNode: makePosition(g.1),
                                             perihelion: makePosition(g.2), aphelion: makePosition(g.3)))
    }
}

// MARK: - Eclipses

extension Ephemera {
    /// The next (or previous, direction -1) solar eclipse anywhere; TT.
    public static func solarEclipse(jdTT: Double, direction: Int = 1) -> SolarEclipse? {
        ready()
        var e = eph_solar_eclipse()
        guard eph_solar_eclipse_next(jdTT, Int32(direction), &e) == 0 else { return nil }
        return SolarEclipse(type: EclipseType(rawValue: e.type), maximum: e.maximum, gamma: e.gamma,
                            magnitude: e.magnitude, latitude: e.latitude * deg, longitude: e.longitude * deg,
                            saros: Int(e.saros), partialBegin: e.partial_begin, partialEnd: e.partial_end,
                            centralBegin: e.central_begin, centralEnd: e.central_end,
                            centerLineBegin: e.center_line_begin, centerLineEnd: e.center_line_end)
    }

    /// The next solar eclipse seen from a place (check `altitude`: the Sun
    /// may be below the horizon).
    public static func solarEclipseLocal(jdTT: Double, longitude: Double, latitude: Double,
                                         elevation: Double = 0, direction: Int = 1) -> LocalSolarEclipse? {
        ready()
        var obs = observer(longitude, latitude, elevation)
        var e = eph_local_solar_eclipse()
        guard eph_solar_eclipse_local(jdTT, &obs, Int32(direction), &e) == 0 else { return nil }
        return LocalSolarEclipse(type: EclipseType(rawValue: e.type), maximum: e.maximum,
                                 firstContact: e.first_contact, secondContact: e.second_contact,
                                 thirdContact: e.third_contact, fourthContact: e.fourth_contact,
                                 magnitude: e.magnitude, obscuration: e.obscuration,
                                 altitude: e.altitude * deg, saros: Int(e.saros))
    }

    /// Where the shadow axis meets the Earth; nil if it misses.
    public static func solarEclipseWhere(jdTT: Double) -> GeoPoint? {
        ready()
        var lat = 0.0, lon = 0.0
        guard eph_solar_eclipse_where(jdTT, &lat, &lon) == 0 else { return nil }
        return GeoPoint(latitude: lat * deg, longitude: lon * deg)
    }

    /// The solar eclipse seen from a place at an instant (empty type: none).
    public static func solarEclipseHow(jdTT: Double, longitude: Double, latitude: Double,
                                       elevation: Double = 0) -> LocalEclipseState? {
        ready()
        var obs = observer(longitude, latitude, elevation)
        var s = eph_local_eclipse_state()
        guard eph_solar_eclipse_how(jdTT, &obs, &s) == 0 else { return nil }
        return LocalEclipseState(type: EclipseType(rawValue: s.type), magnitude: s.magnitude,
                                 obscuration: s.obscuration, altitude: s.altitude * deg)
    }

    /// The next (or previous) lunar eclipse, penumbral ones included; TT.
    public static func lunarEclipse(jdTT: Double, direction: Int = 1) -> LunarEclipse? {
        ready()
        var e = eph_lunar_eclipse()
        guard eph_lunar_eclipse_next(jdTT, Int32(direction), &e) == 0 else { return nil }
        return LunarEclipse(type: EclipseType(rawValue: e.type), maximum: e.maximum, gamma: e.gamma,
                            umbralMagnitude: e.umbral_magnitude, penumbralMagnitude: e.penumbral_magnitude,
                            saros: Int(e.saros), penumbralBegin: e.penumbral_begin, penumbralEnd: e.penumbral_end,
                            partialBegin: e.partial_begin, partialEnd: e.partial_end,
                            totalBegin: e.total_begin, totalEnd: e.total_end)
    }
}

// MARK: - Minor planets

extension Ephemera {
    /// The body of a loaded minor planet: "433", "Eros" or "433 Eros"; nil
    /// if none is loaded.  Works with the positions, searches, phenomena,
    /// orbits and heliocentric positions.
    public static func minorBody(_ name: String) -> Body? {
        ready()
        let i = eph_minor_find(name)
        return i < 0 ? nil : Body(rawValue: Body.minorBase + i)
    }

    public static var minorCount: Int { ready(); return Int(eph_minor_count()) }

    public static func minorList() -> [MinorPlanet] {
        ready()
        return (0..<eph_minor_count()).map {
            MinorPlanet(body: Body(rawValue: Body.minorBase + $0), number: Int(eph_minor_number($0)),
                        name: string(eph_minor_name($0)))
        }
    }
}

// MARK: - Fixed stars

extension Ephemera {
    /// Star index by name ("Regulus"), designation ("alf Leo") or "HIP 49669";
    /// nil if unknown (or the stars are not loaded).
    public static func starFind(_ name: String) -> Int? {
        ready()
        let i = eph_star_find(name)
        return i < 0 ? nil : Int(i)
    }

    /// Stars loaded (0 without eph-stars.bin); indices 0 ..< count.
    public static var starCount: Int { ready(); return Int(eph_star_count()) }
    public static func starName(_ index: Int) -> String { string(eph_star_name(Int32(index))) }
    public static func starDesignation(_ index: Int) -> String { string(eph_star_designation(Int32(index))) }
    public static func starMagnitude(_ index: Int) -> Double { eph_star_magnitude(Int32(index)) }

    public static func starPosition(_ index: Int, jdUT: Double, year: Double? = nil) -> Position? {
        ready()
        let y = year ?? Self.yearOf(jdUT)
        return readPosition { eph_star_position_apparent_ut(Int32(index), jdUT, y, $0) }
    }
}

// MARK: - Ayanamsas

extension Ephemera {
    /// Ayanamsa in degrees: sidereal = tropical - ayanamsa.
    public static func ayanamsa(_ system: Ayanamsa, jdTT: Double) -> Double {
        ready()
        return eph_ayanamsa(eph_ayanamsa_system(rawValue: UInt32(system.rawValue)), jdTT) * deg
    }

    /// A user-defined ayanamsa (degrees at jdT0) carried by the precession.
    public static func ayanamsaCustom(jdT0: Double, valueAtT0: Double, jdTT: Double) -> Double {
        eph_ayanamsa_custom(jdT0, valueAtT0 * rad, jdTT) * deg
    }

    /// "Lahiri", "Fagan/Bradley", ...
    public static func name(of system: Ayanamsa) -> String {
        string(eph_ayanamsa_name(eph_ayanamsa_system(rawValue: UInt32(system.rawValue))))
    }
}

// MARK: - Utilities

extension Ephemera {
    /// Midpoint of the shorter arc between two longitudes, degrees.
    public static func midpoint(_ a: Double, _ b: Double) -> Double { eph_midpoint(a * rad, b * rad) * deg }

    /// Zodiac sign (0 = Aries ... 11) and the position within it.
    public static func splitLongitude(_ longitude: Double) -> SplitLongitude {
        var sign: Int32 = 0, d: Int32 = 0, m: Int32 = 0, s = 0.0
        eph_split_longitude(longitude * rad, &sign, &d, &m, &s)
        return SplitLongitude(sign: Int(sign), degrees: Int(d), minutes: Int(m), seconds: s)
    }
}
