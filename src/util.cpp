#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/* Time helpers and small angle utilities. */

namespace {

/* TAI - UTC (seconds) from each date (1 January / 1 July, 0h UTC) as JD. */
struct LeapSecond {
    double jd;
    int tai_utc;
};
const LeapSecond kLeapSeconds[] = {
    {2441317.5, 10}, {2441499.5, 11}, {2441683.5, 12}, {2442048.5, 13},
    {2442413.5, 14}, {2442778.5, 15}, {2443144.5, 16}, {2443509.5, 17},
    {2443874.5, 18}, {2444239.5, 19}, {2444786.5, 20}, {2445151.5, 21},
    {2445516.5, 22}, {2446247.5, 23}, {2447161.5, 24}, {2447892.5, 25},
    {2448257.5, 26}, {2448804.5, 27}, {2449169.5, 28}, {2449534.5, 29},
    {2450083.5, 30}, {2450630.5, 31}, {2451179.5, 32}, {2453736.5, 33},
    {2454832.5, 34}, {2456109.5, 35}, {2457204.5, 36}, {2457754.5, 37},
};

}  // namespace

extern "C" {

double eph_sidereal_time(double jd_ut) {
    return eph::normalize(eph::gast_degrees(jd_ut) * eph::DEG);
}

int eph_utc_to_jd(int year, int month, int day, int hour, int minute,
                  double second, double* jd_tt, double* jd_ut) {
    if (month < 1 || month > 12 || day < 1 || day > 31 || hour < 0 ||
        hour > 23 || minute < 0 || minute > 59 || second < 0.0 ||
        second >= 61.0)
        return -1;
    int y2, m2;
    double d2;
    const double jd0 = eph_julian_date(year, month, day);
    eph_calendar_date(jd0, &y2, &m2, &d2);
    if (y2 != year || m2 != month || (int)d2 != day) return -1; /* e.g. 31 April */
    /* A leap second (second >= 60) only exists at the end of a day before
     * a table entry. */
    const double jd = jd0 + (hour + (minute + second / 60.0) / 60.0) / 24.0;
    const double year_frac = year + (jd0 - eph_julian_date(year, 1, 1.0)) / 365.25;
    if (jd0 < kLeapSeconds[0].jd) {
        /* Before 1972: UTC was kept within ~0.1 s of UT. */
        if (jd_ut) *jd_ut = jd;
        if (jd_tt) *jd_tt = eph_ut_to_tt(jd, year_frac);
        return 0;
    }
    int leap = kLeapSeconds[0].tai_utc;
    for (const LeapSecond& l : kLeapSeconds)
        if (jd0 >= l.jd) leap = l.tai_utc;
    /* TT = TAI + 32.184 s exactly; UT1 from Delta-T (within ~0.1 s). */
    const double tt = jd + (leap + 32.184) / 86400.0;
    if (jd_tt) *jd_tt = tt;
    if (jd_ut) *jd_ut = eph_tt_to_ut(tt, year_frac);
    return 0;
}

int eph_equation_of_time(double jd_ut, double* e_days) {
    if (!e_days) return -1;
    /* Apparent minus mean solar time: the Sun's true hour angle at
     * Greenwich + 12 h, minus the UT of day. */
    int y, m;
    double d;
    eph_calendar_date(jd_ut, &y, &m, &d);
    const double year = y + (jd_ut - eph_julian_date(y, 1, 1.0)) / 365.25;
    const double jd_tt = eph_ut_to_tt(jd_ut, year);
    eph_position sun;
    if (eph_position_apparent_tt(EPH_SUN, jd_tt, &sun) != 0) return -1;
    double ra, dec;
    eph_ecliptic_to_equatorial(sun.longitude, sun.latitude, jd_tt, &ra, &dec);
    const double hour_angle = eph_sidereal_time(jd_ut) - ra;
    const double mean = (jd_ut + 0.5 - std::floor(jd_ut + 0.5)) * eph::TWO_PI;
    *e_days = (eph::normalize(hour_angle + eph::PI - mean + eph::PI) - eph::PI) /
              eph::TWO_PI;
    return 0;
}

double eph_midpoint(double a, double b) {
    /* Midpoint of the shorter arc from a to b. */
    const double d = eph::normalize(b - a + eph::PI) - eph::PI;
    return eph::normalize(a + d / 2.0);
}

void eph_split_longitude(double longitude, int* sign, int* degrees,
                         int* minutes, double* seconds) {
    /* Rounding to 1e-6" first keeps 29 deg 59' 60" from appearing. */
    double s = std::round(eph::normalize(longitude) / eph::ARCSEC * 1e6) / 1e6;
    if (s >= 1296000.0) s -= 1296000.0;
    const int whole = (int)std::floor(s);
    if (sign) *sign = whole / 108000;
    if (degrees) *degrees = (whole % 108000) / 3600;
    if (minutes) *minutes = (whole % 3600) / 60;
    if (seconds) *seconds = s - std::floor(s / 60.0) * 60.0;
}

}  // extern "C"
