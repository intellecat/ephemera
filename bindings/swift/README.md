# Ephemera for Swift

The Swift package for iOS, macOS, tvOS, watchOS, visionOS and Linux. Its
`Package.swift` is at the root of the repository: the C++ library is compiled
from `src/`, and the 1″ planetary series with Chiron (1550–2650) is bundled
as a resource and loaded on first use.

In Xcode: **File → Add Package Dependencies…** with
`https://github.com/intellecat/ephemera`. In `Package.swift`:

```swift
.package(url: "https://github.com/intellecat/ephemera", from: "0.4.1"),
// and in the target: .product(name: "Ephemera", package: "ephemera")
```

```swift
import Ephemera

let jd = Ephemera.julianDate(2026, 1, 1.5)                 // 1 Jan 2026, 12:00 UT
let moon = Ephemera.position(.moon, jdUT: jd)!             // degrees, AU, per day
let cusps = Ephemera.houses(jdUT: jd, longitude: -0.1278, latitude: 51.5074, system: .placidus)
let lahiri = Ephemera.ayanamsa(.lahiri, jdTT: Ephemera.utToTT(jd))
let eclipse = Ephemera.solarEclipse(jdTT: Ephemera.utToTT(jd))
let now = Ephemera.julianDates(Date())                     // TT and UT of an instant
```

Angles are degrees, distances AU, times Julian days (UT or TT as named);
calls that find nothing return `nil`. Other data files (fixed stars,
asteroids, DE440 tables, minor planets) load with `Ephemera.loadData(_:)` or
`loadData(contentsOf:)`; see [data files](../../docs/data.md). The data
loaded is shared by the whole process.

Tests: `swift test` (or `make swift-test`) at the repository root.
