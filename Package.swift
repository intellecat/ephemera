// swift-tools-version:5.9
// Ephemera for Swift (iOS, macOS, tvOS, watchOS, visionOS, Linux).  The C++
// library is compiled from src/; the Swift API is in bindings/swift.
import PackageDescription

let package = Package(
    name: "Ephemera",
    platforms: [.iOS(.v13), .macOS(.v10_15), .tvOS(.v13), .watchOS(.v6), .visionOS(.v1)],
    products: [
        .library(name: "Ephemera", targets: ["Ephemera"]),
    ],
    targets: [
        .target(
            name: "CEphemera",
            path: ".",
            sources: ["src", "data/delta_t_data.cpp"],
            publicHeadersPath: "include",
            cxxSettings: [.headerSearchPath("src"), .headerSearchPath("data")]
        ),
        .target(
            name: "Ephemera",
            dependencies: ["CEphemera"],
            path: "bindings/swift/Sources/Ephemera",
            resources: [.copy("Resources/eph-1as-chiron.bin")]
        ),
        .testTarget(
            name: "EphemeraTests",
            dependencies: ["Ephemera"],
            path: "bindings/swift/Tests/EphemeraTests"
        ),
    ],
    cxxLanguageStandard: .cxx17
)
