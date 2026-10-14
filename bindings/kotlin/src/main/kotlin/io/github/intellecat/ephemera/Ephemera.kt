package io.github.intellecat.ephemera

import java.io.File
import java.util.Calendar
import java.util.Date
import java.util.TimeZone
import kotlin.math.PI

/**
 * Ephemera for Kotlin and Java (Android, JVM): a pocket-sized, MIT-licensed
 * ephemeris for astrology and astronomy.
 *
 * ```kotlin
 * val jd = Ephemera.julianDate(2026, 1, 1.5)                 // 1 Jan 2026, 12:00 UT
 * val moon = Ephemera.position(Body.MOON, jd)!!              // degrees, AU, per day
 * val cusps = Ephemera.houses(jd, -0.1278, 51.5074, HouseSystem.PLACIDUS)
 * ```
 *
 * The 1" planetary series with Chiron (1550-2650) is built in.  Angles are
 * degrees, distances AU, times Julian days (UT or TT as named); calls that
 * find nothing return null.  The data loaded is shared by the whole process:
 * load files before computing from several threads.
 */
object Ephemera {
    private const val DEG = 180.0 / PI
    private const val RAD = PI / 180.0

    private fun year(jd: Double) = 2000.0 + (jd - 2451545.0) / 365.25

    private fun toPosition(v: DoubleArray?, offset: Int = 0): Position? = v?.let {
        Position(it[offset] * DEG, it[offset + 1] * DEG, it[offset + 2],
                 it[offset + 3] * DEG, it[offset + 4] * DEG, it[offset + 5])
    }

    private fun Double.orNull(): Double? = if (isNaN()) null else this

    // ---- data files ----

    /** Load a data file's bytes (copied); false if rejected.  Later files override earlier ones. */
    @JvmStatic fun loadData(bytes: ByteArray): Boolean = Native.loadData(bytes) == 0

    /** Load a data file from a path. */
    @JvmStatic fun loadData(file: File): Boolean = Native.loadFile(file.path) == 0

    /** Unload every data file (the built-in data too: see [loadEmbedded]). */
    @JvmStatic fun unloadData() = Native.unloadData()

    /** Reload the built-in data (the 1" series with Chiron). */
    @JvmStatic fun loadEmbedded(): Boolean = Native.loadEmbedded() == 0

    // ---- time ----

    /** Julian day of a calendar date (Gregorian; Julian before 1582-10-15); `day` carries the fraction. */
    @JvmStatic fun julianDate(year: Int, month: Int, day: Double): Double = Native.julianDate(year, month, day)

    @JvmStatic fun calendarDate(jd: Double): CalendarDate =
        Native.calendarDate(jd).let { CalendarDate(it[0].toInt(), it[1].toInt(), it[2]) }

    /** TT - UT, seconds. */
    @JvmStatic fun deltaT(year: Double): Double = Native.deltaT(year)

    @JvmStatic @JvmOverloads fun utToTt(jdUt: Double, year: Double = year(jdUt)): Double = Native.utToTt(jdUt, year)
    @JvmStatic @JvmOverloads fun ttToUt(jdTt: Double, year: Double = year(jdTt)): Double = Native.ttToUt(jdTt, year)

    /** Civil UTC (leap seconds from 1972) to TT and UT Julian days; null for an invalid date. */
    @JvmStatic @JvmOverloads
    fun utcToJd(year: Int, month: Int, day: Int, hour: Int = 0, minute: Int = 0, second: Double = 0.0): JulianDates? =
        Native.utcToJd(year, month, day, hour, minute, second)?.let { JulianDates(it[0], it[1]) }

    /** TT and UT Julian days of an instant, milliseconds since 1970 (System.currentTimeMillis()). */
    @JvmStatic fun julianDates(epochMillis: Long): JulianDates {
        val c = Calendar.getInstance(TimeZone.getTimeZone("UTC")).apply { timeInMillis = epochMillis }
        return utcToJd(c.get(Calendar.YEAR), c.get(Calendar.MONTH) + 1, c.get(Calendar.DAY_OF_MONTH),
                       c.get(Calendar.HOUR_OF_DAY), c.get(Calendar.MINUTE),
                       c.get(Calendar.SECOND) + c.get(Calendar.MILLISECOND) / 1000.0)!!
    }

    /** TT and UT Julian days of an instant. */
    @JvmStatic fun julianDates(date: Date): JulianDates = julianDates(date.time)

    /** TDB - TT, seconds. */
    @JvmStatic fun tdbMinusTt(jdTt: Double): Double = Native.tdbMinusTt(jdTt)

    /** Greenwich apparent sidereal time (IAU 2006), degrees. */
    @JvmStatic fun siderealTime(jdUt: Double): Double = Native.siderealTime(jdUt) * DEG

    /** Apparent minus mean solar time, minutes; null without the Sun's data. */
    @JvmStatic fun equationOfTime(jdUt: Double): Double? = Native.equationOfTime(jdUt).orNull()?.let { it * 1440 }

    // ---- frames ----

    /** Mean obliquity of the ecliptic (IAU 2006), degrees. */
    @JvmStatic fun meanObliquity(jdTt: Double): Double = Native.meanObliquity(jdTt) * DEG

    /** Nutation in longitude and obliquity (IAU 2000B), degrees. */
    @JvmStatic fun nutation(jdTt: Double): Nutation = Native.nutation(jdTt).let { Nutation(it[0] * DEG, it[1] * DEG) }

    /** Apparent ecliptic to equatorial coordinates (true obliquity), degrees. */
    @JvmStatic fun eclipticToEquatorial(lon: Double, lat: Double, jdTt: Double): Equatorial =
        Native.eclipticToEquatorial(lon * RAD, lat * RAD, jdTt).let { Equatorial(it[0] * DEG, it[1] * DEG) }

    // ---- positions ----

    /** Geocentric ecliptic position of date for a UT Julian day, apparent by default; null without data. */
    @JvmStatic @JvmOverloads
    fun position(body: Body, jdUt: Double, apparent: Boolean = true, year: Double = year(jdUt)): Position? =
        toPosition(Native.position(body.id, jdUt, year, if (apparent) 3 else 2))

    @JvmStatic @JvmOverloads
    fun positionTt(body: Body, jdTt: Double, apparent: Boolean = true): Position? =
        toPosition(Native.position(body.id, jdTt, 0.0, if (apparent) 1 else 0))

    /** Heliocentric ecliptic of date; null for the Sun, the Moon or without data. */
    @JvmStatic fun heliocentric(body: Body, jdTt: Double): Heliocentric? =
        Native.heliocentric(body.id, jdTt)?.let { Heliocentric(it[0] * DEG, it[1] * DEG, it[2]) }

    /** "Sun", "Moon", ... (minor planets: their names). */
    @JvmStatic fun bodyName(body: Body): String =
        if (body.id >= Body.MINOR_BASE) Native.minorName(body.id - Body.MINOR_BASE) else Native.bodyName(body.id)

    /** Apparent position seen from a place (east longitude, latitude in degrees, elevation in metres). */
    @JvmStatic @JvmOverloads
    fun topocentric(body: Body, jdUt: Double, lon: Double, lat: Double, elevation: Double = 0.0,
                    year: Double = year(jdUt)): Position? =
        toPosition(Native.topocentric(body.id, jdUt, year, lon * RAD, lat * RAD, elevation))

    // ---- houses ----

    @JvmStatic fun ascMc(jdUt: Double, lon: Double, lat: Double): AscMc =
        Native.ascMc(jdUt, lon * RAD, lat * RAD).let { AscMc(it[0] * DEG, it[1] * DEG) }

    /** The twelve cusps, and whether the system is undefined at this latitude (Porphyry is then used). */
    @JvmStatic @JvmOverloads
    fun housesInfo(jdUt: Double, lon: Double, lat: Double, system: HouseSystem = HouseSystem.PORPHYRY): HousesInfo? =
        Native.houses(jdUt, lon * RAD, lat * RAD, system.ordinal)?.let { v ->
            HousesInfo(List(12) { v[it] * DEG }, v[12] == 1.0)
        }

    /** Twelve cusps (degrees); null on bad arguments. */
    @JvmStatic @JvmOverloads
    fun houses(jdUt: Double, lon: Double, lat: Double, system: HouseSystem = HouseSystem.PORPHYRY): List<Double>? =
        housesInfo(jdUt, lon, lat, system)?.cusps

    /** Ascendant, MC, ARMC, Vertex, equatorial Ascendant, co-Ascendants and polar Ascendant, degrees. */
    @JvmStatic fun houseAngles(jdUt: Double, lon: Double, lat: Double): HouseAngles =
        Native.houseAngles(jdUt, lon * RAD, lat * RAD).map { it * DEG }.let {
            HouseAngles(it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7])
        }

    /** House 1.0 ..< 13.0 of a longitude, given cusps from [houses]. */
    @JvmStatic fun housePosition(cusps: List<Double>, lon: Double): Double =
        Native.housePosition(DoubleArray(12) { cusps[it] * RAD }, lon * RAD)

    // ---- phenomena ----

    /** Phase angle, illuminated fraction, elongation, apparent diameter (degrees) and magnitude. */
    @JvmStatic fun phenomena(body: Body, jdTt: Double): Phenomena? = Native.phenomena(body.id, jdTt)?.let {
        Phenomena(it[0] * DEG, it[1], it[2] * DEG, it[3] * DEG, it[4])
    }

    // ---- searches ----

    /**
     * When the apparent longitude reaches [lon] (returns, ingresses); UT in and out, null if not
     * within [maxDays].  With [ayanamsa] the longitude is sidereal, the ayanamsa taken at the
     * crossing itself (sidereal returns and ingresses).
     */
    @JvmStatic @JvmOverloads
    fun longitudeCrossing(body: Body, lon: Double, jdUt: Double, direction: Int = 1, maxDays: Double = 400.0,
                          ayanamsa: AyanamsaOption? = null): Double? {
        if (ayanamsa == null) return Native.longitudeCrossingUt(body.id, lon * RAD, jdUt, direction, maxDays).orNull()
        val jdTt = utToTt(jdUt)
        val tt = when (ayanamsa) {
            is AyanamsaOption.System ->
                Native.longitudeCrossingSidereal(body.id, lon * RAD, ayanamsa.system.ordinal, jdTt, direction, maxDays)
            is AyanamsaOption.Custom ->
                Native.longitudeCrossingSiderealCustom(body.id, lon * RAD, ayanamsa.jdT0, ayanamsa.valueAtT0 * RAD,
                                                       jdTt, direction, maxDays)
        }
        return tt.orNull()?.let { ttToUt(it) }
    }

    /** When the heliocentric longitude reaches [lon]; TT in and out. */
    @JvmStatic @JvmOverloads
    fun heliocentricCrossing(body: Body, lon: Double, jdTt: Double, direction: Int = 1, maxDays: Double = 40000.0,
                             ayanamsa: AyanamsaOption? = null): Double? = when (ayanamsa) {
        null -> Native.heliocentricCrossing(body.id, lon * RAD, jdTt, direction, maxDays)
        is AyanamsaOption.System ->
            Native.heliocentricCrossingSidereal(body.id, lon * RAD, ayanamsa.system.ordinal, jdTt, direction, maxDays)
        is AyanamsaOption.Custom ->
            Native.heliocentricCrossingSiderealCustom(body.id, lon * RAD, ayanamsa.jdT0, ayanamsa.valueAtT0 * RAD,
                                                      jdTt, direction, maxDays)
    }.orNull()

    /** When the body crosses the ecliptic; TT in and out. */
    @JvmStatic @JvmOverloads
    fun nodeCrossing(body: Body, jdTt: Double, direction: Int = 1, maxDays: Double = 40.0): NodeCrossing? =
        Native.nodeCrossing(body.id, jdTt, direction, maxDays)?.let { NodeCrossing(it[0], it[2] == 1.0, it[1] * DEG) }

    /** When the body turns retrograde or direct; TT in and out. */
    @JvmStatic @JvmOverloads
    fun station(body: Body, jdTt: Double, direction: Int = 1, maxDays: Double = 800.0): Station? =
        Native.station(body.id, jdTt, direction, maxDays)?.let { Station(it[0], it[1] == 1.0) }

    // ---- planetary nodes and apsides ----

    /** The orbit's elements (degrees, AU) and its nodes and apsides; null if unsupported or without data. */
    @JvmStatic @JvmOverloads
    fun orbit(body: Body, jdTt: Double, method: OrbitMethod = OrbitMethod.MEAN): Orbit? =
        Native.orbit(body.id, jdTt, method.ordinal)?.let { v ->
            fun points(base: Int) = OrbitPoints(toPosition(v, base)!!, toPosition(v, base + 6)!!,
                                                toPosition(v, base + 12)!!, toPosition(v, base + 18)!!)
            Orbit(v[0], v[1], v[2] * DEG, v[3] * DEG, v[4] * DEG, points(5), points(29))
        }

    // ---- eclipses ----

    /** The next (or previous, direction -1) solar eclipse anywhere; TT. */
    @JvmStatic @JvmOverloads
    fun solarEclipse(jdTt: Double, direction: Int = 1): SolarEclipse? = Native.solarEclipse(jdTt, direction)?.let {
        SolarEclipse(it[0].toInt(), it[1], it[2], it[3], it[4] * DEG, it[5] * DEG, it[6].toInt(),
                     it[7], it[8], it[9], it[10], it[11], it[12])
    }

    /** The next solar eclipse seen from a place (check `altitude`: the Sun may be below the horizon). */
    @JvmStatic @JvmOverloads
    fun solarEclipseLocal(jdTt: Double, lon: Double, lat: Double, elevation: Double = 0.0,
                          direction: Int = 1): LocalSolarEclipse? =
        Native.solarEclipseLocal(jdTt, lon * RAD, lat * RAD, elevation, direction)?.let {
            LocalSolarEclipse(it[0].toInt(), it[1], it[2], it[3], it[4], it[5], it[6], it[7], it[8] * DEG, it[9].toInt())
        }

    /** Where the shadow axis meets the Earth; null if it misses. */
    @JvmStatic fun solarEclipseWhere(jdTt: Double): GeoPoint? =
        Native.solarEclipseWhere(jdTt)?.let { GeoPoint(it[0] * DEG, it[1] * DEG) }

    /** The solar eclipse seen from a place at an instant (type 0: none). */
    @JvmStatic @JvmOverloads
    fun solarEclipseHow(jdTt: Double, lon: Double, lat: Double, elevation: Double = 0.0): LocalEclipseState? =
        Native.solarEclipseHow(jdTt, lon * RAD, lat * RAD, elevation)?.let {
            LocalEclipseState(it[0].toInt(), it[1], it[2], it[3] * DEG)
        }

    /** The next (or previous) lunar eclipse, penumbral ones included; TT. */
    @JvmStatic @JvmOverloads
    fun lunarEclipse(jdTt: Double, direction: Int = 1): LunarEclipse? = Native.lunarEclipse(jdTt, direction)?.let {
        LunarEclipse(it[0].toInt(), it[1], it[2], it[3], it[4], it[5].toInt(), it[6], it[7], it[8], it[9],
                     it[10], it[11])
    }

    // ---- minor planets ----

    /** The body of a loaded minor planet: "433", "Eros" or "433 Eros"; null if none is loaded. */
    @JvmStatic fun minorBody(name: String): Body? =
        Native.minorFind(name).takeIf { it >= 0 }?.let { Body(Body.MINOR_BASE + it) }

    @JvmStatic fun minorCount(): Int = Native.minorCount()

    @JvmStatic fun minorList(): List<MinorPlanet> = List(Native.minorCount()) {
        MinorPlanet(Body(Body.MINOR_BASE + it), Native.minorNumber(it), Native.minorName(it))
    }

    // ---- fixed stars ----

    /** Star index by name ("Regulus"), designation ("alf Leo") or "HIP 49669"; null if unknown. */
    @JvmStatic fun starFind(name: String): Int? = Native.starFind(name).takeIf { it >= 0 }

    /** Stars loaded (0 without eph-stars.bin); indices 0 until the count. */
    @JvmStatic fun starCount(): Int = Native.starCount()
    @JvmStatic fun starName(index: Int): String = Native.starName(index)
    @JvmStatic fun starDesignation(index: Int): String = Native.starDesignation(index)
    @JvmStatic fun starMagnitude(index: Int): Double = Native.starMagnitude(index)

    @JvmStatic @JvmOverloads
    fun starPosition(index: Int, jdUt: Double, year: Double = year(jdUt)): Position? =
        toPosition(Native.starPosition(index, jdUt, year))

    // ---- ayanamsas ----

    /** Ayanamsa in degrees: sidereal = tropical - ayanamsa. */
    @JvmStatic fun ayanamsa(system: Ayanamsa, jdTt: Double): Double = Native.ayanamsa(system.ordinal, jdTt) * DEG

    /** A user-defined ayanamsa (degrees at jdT0) carried by the precession. */
    @JvmStatic fun ayanamsaCustom(jdT0: Double, valueAtT0: Double, jdTt: Double): Double =
        Native.ayanamsaCustom(jdT0, valueAtT0 * RAD, jdTt) * DEG

    /** "Lahiri", "Fagan/Bradley", ... */
    @JvmStatic fun ayanamsaName(system: Ayanamsa): String = Native.ayanamsaName(system.ordinal)

    // ---- utilities ----

    /** Midpoint of the shorter arc between two longitudes, degrees. */
    @JvmStatic fun midpoint(a: Double, b: Double): Double = Native.midpoint(a * RAD, b * RAD) * DEG

    /** Zodiac sign (0 = Aries .. 11) and the position within it. */
    @JvmStatic fun splitLongitude(lon: Double): SplitLongitude =
        Native.splitLongitude(lon * RAD).let { SplitLongitude(it[0].toInt(), it[1].toInt(), it[2].toInt(), it[3]) }
}
