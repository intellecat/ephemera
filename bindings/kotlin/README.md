# Ephemera for Kotlin and Java

The Android library (and the same Kotlin API on the desktop JVM): the C++
library compiled from `src/` with the 1″ planetary series and Chiron
(1550–2650) built in, JNI glue, and the Kotlin API in
`io.github.intellecat.ephemera`. Android 5.0 (API 21) and later;
arm64-v8a, armeabi-v7a, x86_64 and x86.

```kotlin
import io.github.intellecat.ephemera.*

val jd = Ephemera.julianDate(2026, 1, 1.5)                    // 1 Jan 2026, 12:00 UT
val moon = Ephemera.position(Body.MOON, jd)!!                 // degrees, AU, per day
val cusps = Ephemera.houses(jd, -0.1278, 51.5074, HouseSystem.PLACIDUS)
val lahiri = Ephemera.ayanamsa(Ayanamsa.LAHIRI, Ephemera.utToTt(jd))
val eclipse = Ephemera.solarEclipse(Ephemera.utToTt(jd))
val now = Ephemera.julianDates(System.currentTimeMillis())    // TT and UT of an instant
```

From Java the same calls are static: `Ephemera.position(Body.MOON, jd)`.

Angles are degrees, distances AU, times Julian days (UT or TT as named);
calls that find nothing return `null`. Other data files (fixed stars,
asteroids, DE440 tables, minor planets) load with `Ephemera.loadData(bytes)`
or `loadData(file)`, e.g. from an app's assets; see
[data files](../../docs/data.md). The data loaded is shared by the whole
process.

Build the AAR with `./gradlew assembleRelease` here (needs the Android SDK
and NDK); `make kotlin-test` at the repository root runs the checks on the
desktop JVM.
