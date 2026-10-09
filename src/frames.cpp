#include "internal.h"

namespace eph {

Vec3 spherical_to_cartesian(double lon, double lat, double r) {
    const double cl = std::cos(lat);
    return {r * cl * std::cos(lon), r * cl * std::sin(lon), r * std::sin(lat)};
}

void cartesian_to_spherical(const Vec3& v, double& lon, double& lat, double& r) {
    r = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    lon = normalize(std::atan2(v.y, v.x));
    lat = (r > 0.0) ? std::asin(v.z / r) : 0.0;
}

/*
 * Precession of ecliptic coordinates from J2000.0 to the mean equinox and
 * ecliptic of date, IAU 1976 (FK5) model.
 *
 * Route: ecliptic J2000 -> equatorial J2000 -> precess (accumulated
 * angles zeta, z, theta) -> equatorial of date -> ecliptic of date.
 * This reproduces IMCCE's own VSOP87B(J2000) to VSOP87D(of date) reduction
 * to better than 0.1 arcsec in both coordinates, whereas the classical
 * Meeus (21.5) low-order ecliptic formula develops large latitude errors
 * at high ecliptic latitude.
 *
 * Precession angles (Lieske 1979; ERFA prec76), t in Julian centuries:
 *   zeta  = (2306.2181 + 0.30188 t + 0.017998 t^2) t   arcsec
 *   z     = (2306.2181 + 1.09468 t + 0.018203 t^2) t   arcsec
 *   theta = (2004.3109 - 0.42665 t - 0.041833 t^2) t   arcsec
 * Matrix: V(date) = Rz(z) Ry(-theta) Rz(zeta) V(J2000).
 */
void precess_ecliptic_j2000_to_date(double jd_tt, double& lon, double& lat) {
    const double t = (jd_tt - J2000) / DAYS_PER_CENTURY;

    const double zeta = (2306.2181 + (0.30188 + 0.017998 * t) * t) * t * ARCSEC;
    const double z = (2306.2181 + (1.09468 + 0.018203 * t) * t) * t * ARCSEC;
    const double theta =
        (2004.3109 + (-0.42665 - 0.041833 * t) * t) * t * ARCSEC;

    /* Ecliptic J2000 -> equatorial J2000 (unit vector). */
    const double cb = std::cos(lat), sb = std::sin(lat);
    const double e0 = obliquity_iau1980(J2000);
    const double c0 = std::cos(e0), s0 = std::sin(e0);
    double x = cb * std::cos(lon);
    double y = cb * std::sin(lon) * c0 - sb * s0;
    double zz = cb * std::sin(lon) * s0 + sb * c0;

    /* Rz(zeta). */
    const double cz = std::cos(zeta), sz = std::sin(zeta);
    double t1 = x * cz - y * sz;
    double t2 = x * sz + y * cz;
    x = t1;
    y = t2;

    /* Ry(-theta). */
    const double ct = std::cos(theta), st = std::sin(theta);
    t1 = x * ct - zz * st;
    t2 = x * st + zz * ct;
    x = t1;
    zz = t2;

    /* Rz(z). */
    const double czz = std::cos(z), szz = std::sin(z);
    t1 = x * czz - y * szz;
    t2 = x * szz + y * czz;
    x = t1;
    y = t2;

    /* Equatorial of date -> ecliptic of date. */
    const double e = obliquity_iau1980(jd_tt);
    const double ce = std::cos(e), se = std::sin(e);
    const double yl = y * ce + zz * se;
    const double zl = -y * se + zz * ce;

    lon = normalize(std::atan2(yl, x));
    lat = std::asin(zl);
}

double obliquity_iau1980(double jd_tt) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    const double seconds =
        84381.448 + T * (-46.8150 + T * (-0.00059 + T * 0.001813));
    return seconds * ARCSEC;
}

double mean_obliquity(double jd_tt) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    const double seconds =
        84381.406 +
        T * (-46.836769 +
             T * (-0.0001831 +
                  T * (0.00200340 + T * (-0.000000576 - T * 0.0000000434))));
    return seconds * ARCSEC;
}

namespace {

struct NutTerm {
    int l, lp, f, d, om;    /* multipliers of the Delaunay arguments */
    double ps, pst, pc;     /* dpsi: (ps + pst*T) sin + pc cos      */
    double ec, ect, es;     /* deps: (ec + ect*T) cos + es sin      */
};

/* IAU 2000B luni-solar nutation, 77 terms, in 0.1 microarcsec: the
 * coefficients of McCarthy & Luzum (2003) as tabulated in ERFA nut00b.c. */
const NutTerm kNutation[] = {
    { 0,  0,  0,  0,  1, -172064161, -174666, 33386, 92052331, 9086, 15377},
    { 0,  0,  2, -2,  2, -13170906, -1675, -13696, 5730336, -3015, -4587},
    { 0,  0,  2,  0,  2, -2276413, -234, 2796, 978459, -485, 1374},
    { 0,  0,  0,  0,  2, 2074554, 207, -698, -897492, 470, -291},
    { 0,  1,  0,  0,  0, 1475877, -3633, 11817, 73871, -184, -1924},
    { 0,  1,  2, -2,  2, -516821, 1226, -524, 224386, -677, -174},
    { 1,  0,  0,  0,  0, 711159, 73, -872, -6750, 0, 358},
    { 0,  0,  2,  0,  1, -387298, -367, 380, 200728, 18, 318},
    { 1,  0,  2,  0,  2, -301461, -36, 816, 129025, -63, 367},
    { 0, -1,  2, -2,  2, 215829, -494, 111, -95929, 299, 132},
    { 0,  0,  2, -2,  1, 128227, 137, 181, -68982, -9, 39},
    {-1,  0,  2,  0,  2, 123457, 11, 19, -53311, 32, -4},
    {-1,  0,  0,  2,  0, 156994, 10, -168, -1235, 0, 82},
    { 1,  0,  0,  0,  1, 63110, 63, 27, -33228, 0, -9},
    {-1,  0,  0,  0,  1, -57976, -63, -189, 31429, 0, -75},
    {-1,  0,  2,  2,  2, -59641, -11, 149, 25543, -11, 66},
    { 1,  0,  2,  0,  1, -51613, -42, 129, 26366, 0, 78},
    {-2,  0,  2,  0,  1, 45893, 50, 31, -24236, -10, 20},
    { 0,  0,  0,  2,  0, 63384, 11, -150, -1220, 0, 29},
    { 0,  0,  2,  2,  2, -38571, -1, 158, 16452, -11, 68},
    { 0, -2,  2, -2,  2, 32481, 0, 0, -13870, 0, 0},
    {-2,  0,  0,  2,  0, -47722, 0, -18, 477, 0, -25},
    { 2,  0,  2,  0,  2, -31046, -1, 131, 13238, -11, 59},
    { 1,  0,  2, -2,  2, 28593, 0, -1, -12338, 10, -3},
    {-1,  0,  2,  0,  1, 20441, 21, 10, -10758, 0, -3},
    { 2,  0,  0,  0,  0, 29243, 0, -74, -609, 0, 13},
    { 0,  0,  2,  0,  0, 25887, 0, -66, -550, 0, 11},
    { 0,  1,  0,  0,  1, -14053, -25, 79, 8551, -2, -45},
    {-1,  0,  0,  2,  1, 15164, 10, 11, -8001, 0, -1},
    { 0,  2,  2, -2,  2, -15794, 72, -16, 6850, -42, -5},
    { 0,  0, -2,  2,  0, 21783, 0, 13, -167, 0, 13},
    { 1,  0,  0, -2,  1, -12873, -10, -37, 6953, 0, -14},
    { 0, -1,  0,  0,  1, -12654, 11, 63, 6415, 0, 26},
    {-1,  0,  2,  2,  1, -10204, 0, 25, 5222, 0, 15},
    { 0,  2,  0,  0,  0, 16707, -85, -10, 168, -1, 10},
    { 1,  0,  2,  2,  2, -7691, 0, 44, 3268, 0, 19},
    {-2,  0,  2,  0,  0, -11024, 0, -14, 104, 0, 2},
    { 0,  1,  2,  0,  2, 7566, -21, -11, -3250, 0, -5},
    { 0,  0,  2,  2,  1, -6637, -11, 25, 3353, 0, 14},
    { 0, -1,  2,  0,  2, -7141, 21, 8, 3070, 0, 4},
    { 0,  0,  0,  2,  1, -6302, -11, 2, 3272, 0, 4},
    { 1,  0,  2, -2,  1, 5800, 10, 2, -3045, 0, -1},
    { 2,  0,  2, -2,  2, 6443, 0, -7, -2768, 0, -4},
    {-2,  0,  0,  2,  1, -5774, -11, -15, 3041, 0, -5},
    { 2,  0,  2,  0,  1, -5350, 0, 21, 2695, 0, 12},
    { 0, -1,  2, -2,  1, -4752, -11, -3, 2719, 0, -3},
    { 0,  0,  0, -2,  1, -4940, -11, -21, 2720, 0, -9},
    {-1, -1,  0,  2,  0, 7350, 0, -8, -51, 0, 4},
    { 2,  0,  0, -2,  1, 4065, 0, 6, -2206, 0, 1},
    { 1,  0,  0,  2,  0, 6579, 0, -24, -199, 0, 2},
    { 0,  1,  2, -2,  1, 3579, 0, 5, -1900, 0, 1},
    { 1, -1,  0,  0,  0, 4725, 0, -6, -41, 0, 3},
    {-2,  0,  2,  0,  2, -3075, 0, -2, 1313, 0, -1},
    { 3,  0,  2,  0,  2, -2904, 0, 15, 1233, 0, 7},
    { 0, -1,  0,  2,  0, 4348, 0, -10, -81, 0, 2},
    { 1, -1,  2,  0,  2, -2878, 0, 8, 1232, 0, 4},
    { 0,  0,  0,  1,  0, -4230, 0, 5, -20, 0, -2},
    {-1, -1,  2,  2,  2, -2819, 0, 7, 1207, 0, 3},
    {-1,  0,  2,  0,  0, -4056, 0, 5, 40, 0, -2},
    { 0, -1,  2,  2,  2, -2647, 0, 11, 1129, 0, 5},
    {-2,  0,  0,  0,  1, -2294, 0, -10, 1266, 0, -4},
    { 1,  1,  2,  0,  2, 2481, 0, -7, -1062, 0, -3},
    { 2,  0,  0,  0,  1, 2179, 0, -2, -1129, 0, -2},
    {-1,  1,  0,  1,  0, 3276, 0, 1, -9, 0, 0},
    { 1,  1,  0,  0,  0, -3389, 0, 5, 35, 0, -2},
    { 1,  0,  2,  0,  0, 3339, 0, -13, -107, 0, 1},
    {-1,  0,  2, -2,  1, -1987, 0, -6, 1073, 0, -2},
    { 1,  0,  0,  0,  2, -1981, 0, 0, 854, 0, 0},
    {-1,  0,  0,  1,  0, 4026, 0, -353, -553, 0, -139},
    { 0,  0,  2,  1,  2, 1660, 0, -5, -710, 0, -2},
    {-1,  0,  2,  4,  2, -1521, 0, 9, 647, 0, 4},
    {-1,  1,  0,  1,  1, 1314, 0, 0, -700, 0, 0},
    { 0, -2,  2, -2,  1, -1283, 0, 0, 672, 0, 0},
    { 1,  0,  2,  2,  1, -1331, 0, 8, 663, 0, 4},
    {-2,  0,  2,  2,  2, 1383, 0, -2, -594, 0, -2},
    {-1,  0,  0,  0,  2, 1405, 0, 4, -610, 0, 2},
    { 1,  1,  2, -2,  2, 1290, 0, 0, -556, 0, 0},
};

}  // namespace

/* IAU 2000B nutation: within ~1 mas of IAU 2000A over 1995-2050 and a few
 * mas over 1850-2150, including the fixed offsets that stand in for the
 * planetary terms. */
void nutation(double jd_tt, double& dpsi, double& deps) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    const double turn = 1296000.0; /* arcsec */

    /* Delaunay arguments (Simon et al. 1994). */
    const double l = std::fmod(485868.249036 + 1717915923.2178 * T, turn) * ARCSEC;
    const double lp = std::fmod(1287104.79305 + 129596581.0481 * T, turn) * ARCSEC;
    const double f = std::fmod(335779.526232 + 1739527262.8478 * T, turn) * ARCSEC;
    const double d = std::fmod(1072260.70369 + 1602961601.2090 * T, turn) * ARCSEC;
    const double om = std::fmod(450160.398036 - 6962890.5431 * T, turn) * ARCSEC;

    double dp = 0.0, de = 0.0;
    for (int i = (int)(sizeof kNutation / sizeof kNutation[0]) - 1; i >= 0; --i) {
        const NutTerm& t = kNutation[i];
        const double arg = t.l * l + t.lp * lp + t.f * f + t.d * d + t.om * om;
        const double s = std::sin(arg), c = std::cos(arg);
        dp += (t.ps + t.pst * T) * s + t.pc * c;
        de += (t.ec + t.ect * T) * c + t.es * s;
    }
    dpsi = dp * 1e-7 * ARCSEC - 0.135e-3 * ARCSEC;
    deps = de * 1e-7 * ARCSEC + 0.388e-3 * ARCSEC;
}

void ecliptic_to_equatorial(double lon, double lat, double jd_tt, double* ra,
                            double* dec) {
    double dpsi, deps;
    nutation(jd_tt, dpsi, deps);
    const double eps = mean_obliquity(jd_tt) + deps;
    const double ce = std::cos(eps), se = std::sin(eps);
    const double sl = std::sin(lon), cl = std::cos(lon);
    const double sb = std::sin(lat), cb = std::cos(lat);
    const double x = cb * cl;
    const double y = ce * cb * sl - se * sb;
    const double z = se * cb * sl + ce * sb;
    *ra = normalize(std::atan2(y, x));
    *dec = std::asin(z);
}

void equatorial_to_ecliptic(double ra, double dec, double jd_tt, double* lon,
                            double* lat) {
    double dpsi, deps;
    nutation(jd_tt, dpsi, deps);
    const double eps = mean_obliquity(jd_tt) + deps;
    const double ce = std::cos(eps), se = std::sin(eps);
    const double x = std::cos(dec) * std::cos(ra);
    const double y = std::cos(dec) * std::sin(ra);
    const double z = std::sin(dec);
    const double yl = y * ce + z * se;
    const double zl = -y * se + z * ce;
    *lon = normalize(std::atan2(yl, x));
    *lat = std::asin(zl);
}

namespace {

struct Mat3 {
    double m[3][3];
};

Mat3 mul(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) r.m[i][j] += a.m[i][k] * b.m[k][j];
    return r;
}

/* Frame rotations, ERFA convention (positive angle rotates the axes). */
Mat3 rot1(double a) {
    const double c = std::cos(a), s = std::sin(a);
    return {{{1, 0, 0}, {0, c, s}, {0, -s, c}}};
}
Mat3 rot2(double a) {
    const double c = std::cos(a), s = std::sin(a);
    return {{{c, 0, -s}, {0, 1, 0}, {s, 0, c}}};
}
Mat3 rot3(double a) {
    const double c = std::cos(a), s = std::sin(a);
    return {{{c, s, 0}, {-s, c, 0}, {0, 0, 1}}};
}

Mat3 transpose(const Mat3& a) {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r.m[i][j] = a.m[j][i];
    return r;
}

}  // namespace

/*
 * IAU-1976 mean ecliptic of date -> IAU-2006 mean ecliptic of date.
 *
 * The 1976 frame is R1(eps80) P76 applied to the J2000 mean equator (which
 * the IAU-1976 reduction takes as the ICRS, i.e. no frame bias).  The 2006
 * ecliptic of date is R3(-psib) R1(phib) R3(gamb) applied to the GCRS, with
 * the Fukushima-Williams angles (ERFA pfw06), which include the frame bias.
 * The rotation is ~0.04" plus the 0.3"/century general precession change.
 */
Vec3 frame_1976_to_2006(double jd_tt, const Vec3& v) {
    const double t = (jd_tt - J2000) / DAYS_PER_CENTURY;

    const double zeta = (2306.2181 + (0.30188 + 0.017998 * t) * t) * t * ARCSEC;
    const double z = (2306.2181 + (1.09468 + 0.018203 * t) * t) * t * ARCSEC;
    const double theta =
        (2004.3109 + (-0.42665 - 0.041833 * t) * t) * t * ARCSEC;
    const Mat3 p76 = mul(rot3(-z), mul(rot2(theta), rot3(-zeta)));
    const Mat3 ecl76 = mul(rot1(obliquity_iau1980(jd_tt)), p76);

    const double gamb =
        (-0.052928 +
         t * (10.556378 +
              t * (0.4932044 +
                   t * (-0.00031238 + t * (-0.000002788 + t * 0.0000000260))))) *
        ARCSEC;
    const double phib =
        (84381.412819 +
         t * (-46.811016 +
              t * (0.0511268 +
                   t * (0.00053289 + t * (-0.000000440 - t * 0.0000000176))))) *
        ARCSEC;
    const double psib =
        (-0.041775 +
         t * (5038.481484 +
              t * (1.5584175 +
                   t * (-0.00018522 + t * (-0.000026452 - t * 0.0000000148))))) *
        ARCSEC;
    const Mat3 ecl06 = mul(rot3(-psib), mul(rot1(phib), rot3(gamb)));

    const Mat3 r = mul(ecl06, transpose(ecl76));
    return {r.m[0][0] * v.x + r.m[0][1] * v.y + r.m[0][2] * v.z,
            r.m[1][0] * v.x + r.m[1][1] * v.y + r.m[1][2] * v.z,
            r.m[2][0] * v.x + r.m[2][1] * v.y + r.m[2][2] * v.z};
}

/* IAU 2006 GMST: the Earth rotation angle plus the accumulated precession
 * in right ascension (Capitaine et al. 2003).  The polynomial is in TT;
 * using UT for it changes the result by ~1e-4". */
double gmst_degrees(double jd_ut) {
    const double d = jd_ut - J2000;
    const double era = 0.7790572732640 + 0.00273781191135448 * d +
                       (d - std::floor(d)) + std::floor(d);
    const double t = d / DAYS_PER_CENTURY;
    const double p = 0.014506 +
                     t * (4612.156534 +
                          t * (1.3915817 +
                               t * (-0.00000044 +
                                    t * (-0.000029956 - t * 0.0000000368))));
    return std::fmod(era, 1.0) * 360.0 + p / 3600.0;
}

double gast_degrees(double jd_ut) {
    double dpsi, deps;
    nutation(jd_ut, dpsi, deps);
    const double ee = dpsi * std::cos(mean_obliquity(jd_ut) + deps);
    return gmst_degrees(jd_ut) + ee / DEG;
}

}  // namespace eph
