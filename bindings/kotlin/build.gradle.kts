// Ephemera for Android: the C++ library (compiled from the repository's src/
// with the 1" series and Chiron built in) with JNI glue and the Kotlin API.
plugins {
    id("com.android.library") version "8.7.3"
    id("org.jetbrains.kotlin.android") version "2.1.0"
    `maven-publish`
}

// JitPack sets GROUP and VERSION (com.github.intellecat, the tag).
group = System.getenv("GROUP") ?: "io.github.intellecat"
version = System.getenv("VERSION") ?: "0.4.1"

android {
    namespace = "io.github.intellecat.ephemera"
    compileSdk = 35
    ndkVersion = "27.2.12479018"

    defaultConfig {
        minSdk = 21
        consumerProguardFiles("consumer-rules.pro")
        ndk {
            abiFilters += listOf("arm64-v8a", "armeabi-v7a", "x86_64", "x86")
        }
        externalNativeBuild {
            cmake {
                arguments += "-DANDROID_STL=c++_static"
            }
        }
    }

    externalNativeBuild {
        cmake {
            path = file("src/main/cpp/CMakeLists.txt")
            version = "3.22.1"
        }
    }

    compileOptions {
        sourceCompatibility = JavaVersion.VERSION_1_8
        targetCompatibility = JavaVersion.VERSION_1_8
    }

    publishing {
        singleVariant("release") {
            withSourcesJar()
        }
    }
}

kotlin {
    compilerOptions {
        jvmTarget.set(org.jetbrains.kotlin.gradle.dsl.JvmTarget.JVM_1_8)
    }
}

publishing {
    publications {
        register<MavenPublication>("release") {
            artifactId = "ephemera"
            afterEvaluate { from(components["release"]) }
            pom {
                name.set("Ephemera")
                description.set("Pocket-sized, MIT-licensed ephemeris for astrology and astronomy.")
                url.set("https://github.com/intellecat/ephemera")
                licenses {
                    license {
                        name.set("MIT")
                        url.set("https://opensource.org/licenses/MIT")
                    }
                }
            }
        }
    }
}
