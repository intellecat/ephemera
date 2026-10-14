package io.github.intellecat.ephemera

/**
 * A body: the planets, the Sun and Moon, asteroids, lunar points, Uranian
 * points, or a loaded minor planet ([Ephemera.minorBody]).
 */
data class Body(val id: Int) {
    companion object {
        @JvmField val MERCURY = Body(0)
        @JvmField val VENUS = Body(1)
        @JvmField val EARTH = Body(2)
        @JvmField val MARS = Body(3)
        @JvmField val JUPITER = Body(4)
        @JvmField val SATURN = Body(5)
        @JvmField val URANUS = Body(6)
        @JvmField val NEPTUNE = Body(7)
        @JvmField val PLUTO = Body(8)
        @JvmField val SUN = Body(9)
        @JvmField val MOON = Body(10)
        @JvmField val CERES = Body(11)
        @JvmField val PALLAS = Body(12)
        @JvmField val JUNO = Body(13)
        @JvmField val VESTA = Body(14)
        @JvmField val CHIRON = Body(15)
        @JvmField val MEAN_NODE = Body(16)
        @JvmField val TRUE_NODE = Body(17)
        @JvmField val MEAN_LILITH = Body(18)
        @JvmField val OSCU_LILITH = Body(19)
        // The Hamburg school's Uranian points (Witte, Sieggrün).
        @JvmField val CUPIDO = Body(20)
        @JvmField val HADES = Body(21)
        @JvmField val ZEUS = Body(22)
        @JvmField val KRONOS = Body(23)
        @JvmField val APOLLON = Body(24)
        @JvmField val ADMETOS = Body(25)
        @JvmField val VULKANUS = Body(26)
        @JvmField val POSEIDON = Body(27)

        /** Body id of the first loaded minor planet. */
        const val MINOR_BASE = 1000
    }
}

/** In the C library's order: [Ayanamsa.ordinal] is its value. */
enum class Ayanamsa {
    LAHIRI, FAGAN_BRADLEY, RAMAN, KRISHNAMURTI, YUKTESHWAR, DE_LUCE, USHA_SHASHI, DJWHAL_KHUL,
    JN_BHASIN, BABYLONIAN_KUGLER1, BABYLONIAN_KUGLER2, BABYLONIAN_KUGLER3, BABYLONIAN_HUBER,
    BABYLONIAN_ETA_PISCIUM, BABYLONIAN_ALDEBARAN, HIPPARCHOS, SASSANIAN, J2000, J1900, B1950,
    SURYASIDDHANTA, SURYASIDDHANTA_MEAN_SUN, ARYABHATA, ARYABHATA_MEAN_SUN, SS_REVATI, SS_CITRA,
    SKYDRAM, ARYABHATA_522, BABYLONIAN_BRITTON, GALACTIC_EQUATOR_FIORENZA, VETTIUS_VALENS,
    LAHIRI_1940, LAHIRI_VP285, KRISHNAMURTI_SENTHILATHIBAN, LAHIRI_ICRC,
    // Star-based (need the stars: eph-stars.bin).
    TRUE_CITRA, TRUE_REVATI, TRUE_PUSHYA, TRUE_MULA, GALACTIC_CENTER, GALACTIC_CENTER_BRAND,
    GALACTIC_CENTER_COCHRANE,
}

/** In the C library's order: [HouseSystem.ordinal] is its value. */
enum class HouseSystem {
    WHOLE_SIGN, EQUAL, PORPHYRY, PLACIDUS, KOCH, REGIOMONTANUS, CAMPANUS, ALCABITIUS, TOPOCENTRIC,
    VEHLOW, EQUAL_MC, EQUAL_ARIES, MERIDIAN, MORINUS,
}

enum class OrbitMethod {
    /** Secular mean elements (Simon et al. 1994 / Meeus ch. 31): Mercury..Neptune. */
    MEAN,
    /** From the position and velocity: the planets but the Earth, Pluto, minor bodies. */
    OSCULATING,
}

/** Eclipse type flags (the `type` of an eclipse result). */
object EclipseType {
    const val TOTAL = 1
    const val ANNULAR = 2
    const val HYBRID = 4
    const val PARTIAL = 8
    const val PENUMBRAL = 16
    const val NON_CENTRAL = 32
}

/** An ayanamsa for the sidereal searches. */
sealed class AyanamsaOption {
    data class System(val system: Ayanamsa) : AyanamsaOption()
    /** A custom ayanamsa: its value (degrees) at a TT Julian day, carried by the precession. */
    data class Custom(val jdT0: Double, val valueAtT0: Double) : AyanamsaOption()
}

data class Position(
    val longitude: Double, val latitude: Double, val distance: Double,
    val longitudeSpeed: Double, val latitudeSpeed: Double, val distanceSpeed: Double,
) {
    val isRetrograde: Boolean get() = longitudeSpeed < 0
}

/** `day` carries the fraction of the day. */
data class CalendarDate(val year: Int, val month: Int, val day: Double)
data class JulianDates(val tt: Double, val ut: Double)
data class Nutation(val longitude: Double, val obliquity: Double)
data class Equatorial(val rightAscension: Double, val declination: Double)
data class Heliocentric(val longitude: Double, val latitude: Double, val distance: Double)
data class AscMc(val asc: Double, val mc: Double)

/** [fallback]: the system is undefined here (polar), Porphyry was used. */
data class HousesInfo(val cusps: List<Double>, val fallback: Boolean)

data class HouseAngles(
    val ascendant: Double, val mc: Double, val armc: Double, val vertex: Double,
    val equatorialAscendant: Double, val coAscendantKoch: Double, val coAscendantMunkasey: Double,
    val polarAscendant: Double,
)

/** [phase]: the illuminated fraction, 0..1. */
data class Phenomena(
    val phaseAngle: Double, val phase: Double, val elongation: Double, val diameter: Double, val magnitude: Double,
)

data class NodeCrossing(val jd: Double, val ascending: Boolean, val longitude: Double)

/** [turnsRetrograde]: true when the body turns retrograde, false when it turns direct. */
data class Station(val jd: Double, val turnsRetrograde: Boolean)

data class OrbitPoints(
    val ascendingNode: Position, val descendingNode: Position, val perihelion: Position, val aphelion: Position,
)

/** A heliocentric orbit, mean ecliptic and equinox of date. */
data class Orbit(
    val semiMajorAxis: Double, val eccentricity: Double, val inclination: Double, val node: Double,
    val perihelion: Double, val heliocentric: OrbitPoints, val geocentric: OrbitPoints,
)

/** Times are TT Julian days; 0 marks a phase that does not occur.  [type]: [EclipseType] flags. */
data class SolarEclipse(
    val type: Int, val maximum: Double, val gamma: Double, val magnitude: Double,
    val latitude: Double, val longitude: Double, val saros: Int,
    val partialBegin: Double, val partialEnd: Double, val centralBegin: Double, val centralEnd: Double,
    val centerLineBegin: Double, val centerLineEnd: Double,
)

/** [altitude]: of the Sun at maximum; below 0, not visible. */
data class LocalSolarEclipse(
    val type: Int, val maximum: Double, val firstContact: Double, val secondContact: Double,
    val thirdContact: Double, val fourthContact: Double, val magnitude: Double, val obscuration: Double,
    val altitude: Double, val saros: Int,
)

/** [type] 0: no eclipse here and now. */
data class LocalEclipseState(val type: Int, val magnitude: Double, val obscuration: Double, val altitude: Double)

data class LunarEclipse(
    val type: Int, val maximum: Double, val gamma: Double, val umbralMagnitude: Double,
    val penumbralMagnitude: Double, val saros: Int, val penumbralBegin: Double, val penumbralEnd: Double,
    val partialBegin: Double, val partialEnd: Double, val totalBegin: Double, val totalEnd: Double,
)

data class GeoPoint(val latitude: Double, val longitude: Double)
data class MinorPlanet(val body: Body, val number: Int, val name: String)

/** [sign]: 0 = Aries .. 11 = Pisces. */
data class SplitLongitude(val sign: Int, val degrees: Int, val minutes: Int, val seconds: Double)
