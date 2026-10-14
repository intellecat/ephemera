package io.github.intellecat.ephemera

/** The JNI functions (ephemera_jni.c).  Radians; use [Ephemera] instead. */
internal object Native {
    init {
        System.loadLibrary("ephemera")
    }

    @JvmStatic external fun loadData(data: ByteArray): Int
    @JvmStatic external fun loadFile(path: String): Int
    @JvmStatic external fun loadEmbedded(): Int
    @JvmStatic external fun unloadData()

    @JvmStatic external fun julianDate(year: Int, month: Int, day: Double): Double
    @JvmStatic external fun calendarDate(jd: Double): DoubleArray
    @JvmStatic external fun deltaT(year: Double): Double
    @JvmStatic external fun utToTt(jd: Double, year: Double): Double
    @JvmStatic external fun ttToUt(jd: Double, year: Double): Double
    @JvmStatic external fun utcToJd(year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double): DoubleArray?
    @JvmStatic external fun tdbMinusTt(jd: Double): Double
    @JvmStatic external fun siderealTime(jd: Double): Double
    @JvmStatic external fun equationOfTime(jd: Double): Double
    @JvmStatic external fun meanObliquity(jd: Double): Double
    @JvmStatic external fun nutation(jd: Double): DoubleArray
    @JvmStatic external fun eclipticToEquatorial(lon: Double, lat: Double, jd: Double): DoubleArray

    @JvmStatic external fun position(body: Int, jd: Double, year: Double, mode: Int): DoubleArray?
    @JvmStatic external fun heliocentric(body: Int, jd: Double): DoubleArray?
    @JvmStatic external fun bodyName(body: Int): String
    @JvmStatic external fun topocentric(body: Int, jd: Double, year: Double, lon: Double, lat: Double, elevation: Double): DoubleArray?

    @JvmStatic external fun ascMc(jd: Double, lon: Double, lat: Double): DoubleArray
    @JvmStatic external fun houses(jd: Double, lon: Double, lat: Double, system: Int): DoubleArray?
    @JvmStatic external fun houseAngles(jd: Double, lon: Double, lat: Double): DoubleArray
    @JvmStatic external fun housePosition(cusps: DoubleArray, lon: Double): Double

    @JvmStatic external fun phenomena(body: Int, jd: Double): DoubleArray?

    @JvmStatic external fun longitudeCrossingUt(body: Int, lon: Double, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun longitudeCrossingSidereal(body: Int, lon: Double, system: Int, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun longitudeCrossingSiderealCustom(body: Int, lon: Double, jd0: Double, value0: Double, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun heliocentricCrossing(body: Int, lon: Double, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun heliocentricCrossingSidereal(body: Int, lon: Double, system: Int, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun heliocentricCrossingSiderealCustom(body: Int, lon: Double, jd0: Double, value0: Double, jd: Double, dir: Int, max: Double): Double
    @JvmStatic external fun nodeCrossing(body: Int, jd: Double, dir: Int, max: Double): DoubleArray?
    @JvmStatic external fun station(body: Int, jd: Double, dir: Int, max: Double): DoubleArray?

    @JvmStatic external fun orbit(body: Int, jd: Double, method: Int): DoubleArray?

    @JvmStatic external fun solarEclipse(jd: Double, dir: Int): DoubleArray?
    @JvmStatic external fun solarEclipseLocal(jd: Double, lon: Double, lat: Double, elevation: Double, dir: Int): DoubleArray?
    @JvmStatic external fun solarEclipseWhere(jd: Double): DoubleArray?
    @JvmStatic external fun solarEclipseHow(jd: Double, lon: Double, lat: Double, elevation: Double): DoubleArray?
    @JvmStatic external fun lunarEclipse(jd: Double, dir: Int): DoubleArray?

    @JvmStatic external fun minorCount(): Int
    @JvmStatic external fun minorFind(name: String): Int
    @JvmStatic external fun minorNumber(index: Int): Int
    @JvmStatic external fun minorName(index: Int): String

    @JvmStatic external fun starCount(): Int
    @JvmStatic external fun starFind(name: String): Int
    @JvmStatic external fun starName(index: Int): String
    @JvmStatic external fun starDesignation(index: Int): String
    @JvmStatic external fun starMagnitude(index: Int): Double
    @JvmStatic external fun starPosition(index: Int, jd: Double, year: Double): DoubleArray?

    @JvmStatic external fun ayanamsa(system: Int, jd: Double): Double
    @JvmStatic external fun ayanamsaCustom(jd0: Double, value0: Double, jd: Double): Double
    @JvmStatic external fun ayanamsaName(system: Int): String

    @JvmStatic external fun midpoint(a: Double, b: Double): Double
    @JvmStatic external fun splitLongitude(lon: Double): DoubleArray
}
