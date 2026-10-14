// Smoke test of the Kotlin binding on the desktop JVM (make kotlin-test):
// the same Kotlin sources as the Android library, with the JNI library
// built for this machine.  Reference values come from the C library's
// tests and the npm, Python and Swift packages'.
import io.github.intellecat.ephemera.*
import java.io.File
import kotlin.math.abs
import kotlin.system.exitProcess

var fails = 0
fun check(ok: Boolean, what: String) { if (!ok) { println("FAIL: $what"); fails++ } }
fun near(a: Double?, b: Double, tol: Double) = a != null && abs(a - b) <= tol

fun main(args: Array<String>) {
    val data = File(args.getOrElse(0) { "data" })
    val jd = Ephemera.julianDate(2026, 1, 1.5)
    check(jd == 2461042.0, "julianDate")
    check(Ephemera.calendarDate(jd) == CalendarDate(2026, 1, 1.5), "calendarDate")
    val j = Ephemera.utcToJd(2026, 1, 1, 12)!!
    check(near((j.tt - j.ut) * 86400, Ephemera.deltaT(2026.0), 1.0), "utcToJd")
    check(Ephemera.utcToJd(2026, 2, 30) == null, "invalid date")
    check(near(Ephemera.julianDates(((jd - 2440587.5) * 86400000).toLong()).ut, j.ut, 1e-6),
          "julianDates(epochMillis)")

    check(near(Ephemera.position(Body.MOON, jd)?.longitude, 74.2252, 1e-3), "Moon")
    for (b in listOf(Body.SUN, Body.MERCURY, Body.PLUTO, Body.CHIRON, Body.TRUE_NODE, Body.MEAN_LILITH, Body.CUPIDO))
        check(Ephemera.position(b, jd) != null, "position $b")
    check(Ephemera.position(Body.CERES, jd) == null, "no Ceres without data")
    check(near(Ephemera.heliocentric(Body.EARTH, Ephemera.utToTt(jd))?.distance, 0.983, 0.001), "heliocentric")
    check(Ephemera.bodyName(Body.MOON) == "Moon", "bodyName")

    val cusps = Ephemera.houses(jd, -0.1278, 51.5074, HouseSystem.PLACIDUS)!!
    check(cusps.size == 12 && near(cusps[0], 25.536293992570464, 1e-9) && near(cusps[9], 280.137101405305, 1e-9),
          "Placidus cusps")
    check(near(Ephemera.housePosition(cusps, cusps[0] + 1e-9), 1.0, 1e-6), "housePosition")
    check(near(Ephemera.ayanamsa(Ayanamsa.LAHIRI, 2451545.0), 23.857092325, 1e-6), "Lahiri")
    check(Ephemera.ayanamsaName(Ayanamsa.FAGAN_BRADLEY) == "Fagan/Bradley", "ayanamsaName")

    val sun = Ephemera.position(Body.SUN, jd)!!.longitude
    val ret = Ephemera.longitudeCrossing(Body.SUN, sun, jd + 10)
    check(ret != null && ret - jd > 365 && ret - jd < 366, "solar return")
    check(Ephemera.longitudeCrossing(Body.SUN, 0.0, jd, ayanamsa = AyanamsaOption.System(Ayanamsa.LAHIRI)) != null,
          "sidereal ingress")
    check(Ephemera.station(Body.MERCURY, Ephemera.utToTt(jd)) != null, "station")
    val eclipse = Ephemera.solarEclipse(2460400.0)
    check(eclipse?.saros == 139 && eclipse.type and EclipseType.TOTAL != 0, "solar eclipse")
    check(Ephemera.lunarEclipse(2460400.0) != null, "lunar eclipse")
    check(near(Ephemera.orbit(Body.MARS, Ephemera.utToTt(jd))?.semiMajorAxis, 1.5237, 1e-3), "orbit")

    check(Ephemera.starCount() == 0 && Ephemera.starFind("Regulus") == null, "no stars yet")
    check(Ephemera.loadData(File(data, "eph-stars.bin")), "load stars")
    check(Ephemera.loadData(File(data, "eph-asteroids.bin").readBytes()), "load asteroids")
    check(Ephemera.position(Body.CERES, jd) != null, "Ceres")
    val n = Ephemera.starCount()
    check(n == 411 && (0 until n).all { Ephemera.starFind(Ephemera.starName(it)) == it }, "star listing ($n)")
    val reg = Ephemera.starFind("Regulus")!!
    check(Ephemera.starDesignation(reg) == "alf Leo", "Regulus")
    check(near(Ephemera.starPosition(reg, jd)?.longitude, 150.1, 0.2), "Regulus position")
    Ephemera.unloadData()
    check(Ephemera.position(Body.MOON, jd) == null && Ephemera.loadEmbedded() &&
          Ephemera.position(Body.MOON, jd) != null, "unload and reload the built-in data")

    val s = Ephemera.splitLongitude(123.5)
    check(s.sign == 4 && s.degrees == 3 && s.minutes == 30, "splitLongitude")

    if (fails == 0) println("All Kotlin binding checks passed.")
    exitProcess(if (fails == 0) 0 else 1)
}
