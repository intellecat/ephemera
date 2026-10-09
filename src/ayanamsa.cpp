#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

namespace {

/* General precession in longitude since J2000 (IAU 2006), radians. */
double precession(double jd_tt) {
    const double T = (jd_tt - eph::J2000) / eph::DAYS_PER_CENTURY;
    return T * (5028.796195 +
                T * (1.1054348 +
                     T * (0.00007964 + T * (-0.000023857 - T * 0.0000000383)))) *
           eph::ARCSEC;
}

/* Each system's mean ayanamsa at J2000 (degrees), following its definition
 * (a value at a reference date) as realised by Swiss Ephemeris, e.g.
 * Lahiri: 23 deg 15' 00.658" (true) on 1956-03-21 (Indian Calendar Reform
 * Committee); Fagan/Bradley: SVP 24 deg 02' 31.36" in 1950; J1900: zero in
 * 1900. */
const struct {
    double j2000;
    const char* name;
} kSystems[] = {
    {23.857092325, "Lahiri"},
    {24.740299966, "Fagan/Bradley"},
    {22.410791012, "Raman"},
    {23.760240012, "Krishnamurti (KP)"},
    {22.478803000, "Sri Yukteshwar"},
    {27.815752774, "De Luce"},
    {20.057541000, "Usha/Shashi"},
    {28.359678600, "Djwhal Khul"},
    {22.762137000, "J.N. Bhasin"},
    {23.533639863, "Babylonian (Kugler 1)"},
    {24.933639863, "Babylonian (Kugler 2)"},
    {25.783639863, "Babylonian (Kugler 3)"},
    {24.733639863, "Babylonian (Huber)"},
    {24.522527935, "Babylonian (Eta Piscium)"},
    {24.758923883, "Babylonian (Aldebaran = 15 Tau)"},
    {20.247788094, "Hipparchos"},
    {19.992959383, "Sassanian"},
    {0.000000000, "J2000"},
    {1.396581000, "J1900"},
    {0.698370129, "B1950"},
    {20.895058836, "Suryasiddhanta"},
    {20.680424886, "Suryasiddhanta (mean Sun)"},
    {20.895059722, "Aryabhata"},
    {20.657427342, "Aryabhata (mean Sun)"},
    {20.103388376, "SS Revati"},
    {23.005763276, "SS Citra"},
    {30.017793877, "Skydram (Mardyks)"},
    {20.575846967, "Aryabhata 522"},
    {24.615752774, "Babylonian (Britton)"},
    {25.000019094, "Galactic Equator (Fiorenza)"},
    {22.795609299, "Vettius Valens"},
    {23.842323232, "Lahiri 1940"},
    {23.863481202, "Lahiri VP285"},
    {23.780364957, "Krishnamurti-Senthilathiban"},
    {23.856788988, "Lahiri ICRC"},
};
static_assert(sizeof kSystems / sizeof kSystems[0] == EPH_AYANAMSA_COUNT,
              "one entry per eph_ayanamsa_system");

}  // namespace

extern "C" {

double eph_ayanamsa(eph_ayanamsa_system system, double jd_tt) {
    if (system < 0 || system >= EPH_AYANAMSA_COUNT) return std::nan("");
    return kSystems[system].j2000 * eph::DEG + precession(jd_tt);
}

const char* eph_ayanamsa_name(eph_ayanamsa_system system) {
    if (system < 0 || system >= EPH_AYANAMSA_COUNT) return nullptr;
    return kSystems[system].name;
}

double eph_ayanamsa_custom(double jd_t0, double value_t0, double jd_tt) {
    return value_t0 + precession(jd_tt) - precession(jd_t0);
}

double eph_ayanamsa_lahiri(double jd_tt) {
    return eph_ayanamsa(EPH_AYANAMSA_LAHIRI, jd_tt);
}

}  // extern "C"
