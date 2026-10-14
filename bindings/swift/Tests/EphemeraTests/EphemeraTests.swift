// Smoke tests of the Swift package (swift test).  Reference values come from
// the C library's tests and the npm and Python packages'.
import Ephemera
import Foundation
import XCTest

final class EphemeraTests: XCTestCase {
    let jd = Ephemera.julianDate(2026, 1, 1.5)  // 1 Jan 2026, 12:00 UT

    override func setUp() {
        Ephemera.unloadData()
        XCTAssertTrue(Ephemera.loadEmbedded())
    }

    func testTime() {
        XCTAssertEqual(jd, 2461042.0)
        let date = Ephemera.calendarDate(jd)
        XCTAssertEqual([Double(date.year), Double(date.month), date.day], [2026, 1, 1.5])
        let j = Ephemera.utcToJD(2026, 1, 1, 12)!
        XCTAssertEqual((j.tt - j.ut) * 86400, Ephemera.deltaT(2026), accuracy: 1)
        XCTAssertNil(Ephemera.utcToJD(2026, 2, 30))
        let instant = Date(timeIntervalSince1970: (jd - 2440587.5) * 86400)
        XCTAssertEqual(Ephemera.julianDates(instant).ut, j.ut, accuracy: 1e-6)
    }

    func testPositions() throws {
        let moon = try XCTUnwrap(Ephemera.position(.moon, jdUT: jd))
        XCTAssertEqual(moon.longitude, 74.2252, accuracy: 1e-3)
        for b: Ephemera.Body in [.sun, .mercury, .pluto, .chiron, .trueNode, .meanLilith, .cupido] {
            XCTAssertNotNil(Ephemera.position(b, jdUT: jd), "\(b)")
        }
        XCTAssertNil(Ephemera.position(.ceres, jdUT: jd))
        XCTAssertEqual(Ephemera.heliocentric(.earth, jdTT: Ephemera.utToTT(jd))!.distance, 0.983, accuracy: 0.001)
        XCTAssertEqual(Ephemera.name(of: .moon), "Moon")
    }

    func testHousesAndAyanamsa() throws {
        let cusps = try XCTUnwrap(Ephemera.houses(jdUT: jd, longitude: -0.1278, latitude: 51.5074, system: .placidus))
        XCTAssertEqual(cusps.count, 12)
        XCTAssertEqual(cusps[0], 25.536293992570464, accuracy: 1e-9)
        XCTAssertEqual(cusps[9], 280.137101405305, accuracy: 1e-9)
        XCTAssertEqual(Ephemera.housePosition(cusps: cusps, longitude: cusps[0] + 1e-9), 1.0, accuracy: 1e-6)
        XCTAssertEqual(Ephemera.ayanamsa(.lahiri, jdTT: 2451545.0), 23.857092325, accuracy: 1e-6)
        XCTAssertEqual(Ephemera.name(of: Ephemera.Ayanamsa.faganBradley), "Fagan/Bradley")
    }

    func testSearchesAndEclipses() throws {
        let sun = try XCTUnwrap(Ephemera.position(.sun, jdUT: jd)).longitude
        let ret = try XCTUnwrap(Ephemera.longitudeCrossing(.sun, longitude: sun, jdUT: jd + 10))
        XCTAssertTrue(ret - jd > 365 && ret - jd < 366)
        XCTAssertNotNil(Ephemera.longitudeCrossing(.sun, longitude: 0, jdUT: jd, ayanamsa: .system(.lahiri)))
        XCTAssertNotNil(Ephemera.station(.mercury, jdTT: Ephemera.utToTT(jd)))
        let eclipse = try XCTUnwrap(Ephemera.solarEclipse(jdTT: 2460400.0))
        XCTAssertEqual(eclipse.saros, 139)
        XCTAssertTrue(eclipse.type.contains(.total))
        XCTAssertNotNil(Ephemera.lunarEclipse(jdTT: 2460400.0))
        XCTAssertEqual(Ephemera.orbit(.mars, jdTT: Ephemera.utToTT(jd))!.semiMajorAxis, 1.5237, accuracy: 1e-3)
    }

    func testDataFiles() throws {
        XCTAssertEqual(Ephemera.starCount, 0)
        XCTAssertNil(Ephemera.starFind("Regulus"))
        // The repository's data files, when the tests run from a checkout.
        let root = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .appendingPathComponent("../../../../data").standardized
        let stars = root.appendingPathComponent("eph-stars.bin")
        try XCTSkipUnless(FileManager.default.fileExists(atPath: stars.path), "no data/ checkout")
        XCTAssertTrue(Ephemera.loadData(contentsOf: stars))
        XCTAssertTrue(Ephemera.loadData(try Data(contentsOf: root.appendingPathComponent("eph-asteroids.bin"))))
        XCTAssertNotNil(Ephemera.position(.ceres, jdUT: jd))
        XCTAssertEqual(Ephemera.starCount, 411)
        let reg = try XCTUnwrap(Ephemera.starFind("Regulus"))
        XCTAssertEqual(Ephemera.starDesignation(reg), "alf Leo")
        XCTAssertEqual(Ephemera.starPosition(reg, jdUT: jd)!.longitude, 150.1, accuracy: 0.2)
    }

    func testUtilities() {
        let s = Ephemera.splitLongitude(123.5)
        XCTAssertEqual([s.sign, s.degrees, s.minutes], [4, 3, 30])
        XCTAssertTrue(Ephemera.siderealTime(jd).isFinite)
    }
}
