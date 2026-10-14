/*
 * ephemera_jni.c -- JNI glue for the Kotlin/Java binding
 * (io.github.intellecat.ephemera.Native).  Angles stay in radians here; the
 * Kotlin layer converts to degrees.  Results with several values come back
 * as double arrays (null when the C call fails); searches return NaN when
 * nothing is found.
 */
#include <jni.h>
#include <math.h>
#include <stdlib.h>

#include "ephemeris/ephemeris.h"

#define FN(ret, name) JNIEXPORT ret JNICALL Java_io_github_intellecat_ephemera_Native_##name
#define ARGS JNIEnv* env, jclass cls

static jdoubleArray doubles(JNIEnv* env, const double* v, int n) {
    jdoubleArray a = (*env)->NewDoubleArray(env, n);
    if (a) (*env)->SetDoubleArrayRegion(env, a, 0, n, v);
    return a;
}

static jstring string(JNIEnv* env, const char* s) {
    return (*env)->NewStringUTF(env, s ? s : "");
}

static eph_observer observer(double lon, double lat, double elevation) {
    eph_observer o = {lon, lat, elevation};
    return o;
}

static jdoubleArray position(JNIEnv* env, int rc, const eph_position* p) {
    if (rc != 0) return NULL;
    double v[6] = {p->longitude, p->latitude, p->distance,
                   p->longitude_speed, p->latitude_speed, p->distance_speed};
    return doubles(env, v, 6);
}

/* ---- data files ---- */

FN(jint, loadData)(ARGS, jbyteArray data) {
    (void)cls;
    jsize n = (*env)->GetArrayLength(env, data);
    jbyte* bytes = (*env)->GetByteArrayElements(env, data, NULL);
    if (!bytes) return -1;
    int rc = eph_load_data(bytes, (size_t)n); /* copies */
    (*env)->ReleaseByteArrayElements(env, data, bytes, JNI_ABORT);
    return rc;
}

FN(jint, loadFile)(ARGS, jstring path) {
    (void)cls;
    const char* p = (*env)->GetStringUTFChars(env, path, NULL);
    if (!p) return -1;
    int rc = eph_load_file(p);
    (*env)->ReleaseStringUTFChars(env, path, p);
    return rc;
}

FN(jint, loadEmbedded)(ARGS) { (void)env; (void)cls; return eph_load_embedded(); }
FN(void, unloadData)(ARGS) { (void)env; (void)cls; eph_unload_data(); }

/* ---- time and frames ---- */

FN(jdouble, julianDate)(ARGS, jint y, jint m, jdouble d) { (void)env; (void)cls; return eph_julian_date(y, m, d); }
FN(jdouble, deltaT)(ARGS, jdouble year) { (void)env; (void)cls; return eph_delta_t(year); }
FN(jdouble, utToTt)(ARGS, jdouble jd, jdouble year) { (void)env; (void)cls; return eph_ut_to_tt(jd, year); }
FN(jdouble, ttToUt)(ARGS, jdouble jd, jdouble year) { (void)env; (void)cls; return eph_tt_to_ut(jd, year); }
FN(jdouble, tdbMinusTt)(ARGS, jdouble jd) { (void)env; (void)cls; return eph_tdb_minus_tt(jd); }
FN(jdouble, siderealTime)(ARGS, jdouble jd) { (void)env; (void)cls; return eph_sidereal_time(jd); }
FN(jdouble, meanObliquity)(ARGS, jdouble jd) { (void)env; (void)cls; return eph_mean_obliquity(jd); }

FN(jdoubleArray, calendarDate)(ARGS, jdouble jd) {
    (void)cls;
    int y, m;
    double d;
    eph_calendar_date(jd, &y, &m, &d);
    double v[3] = {y, m, d};
    return doubles(env, v, 3);
}

FN(jdoubleArray, utcToJd)(ARGS, jint y, jint mo, jint d, jint h, jint mi, jdouble s) {
    (void)cls;
    double v[2];
    if (eph_utc_to_jd(y, mo, d, h, mi, s, &v[0], &v[1]) != 0) return NULL;
    return doubles(env, v, 2);
}

FN(jdouble, equationOfTime)(ARGS, jdouble jd) {
    (void)env; (void)cls;
    double e;
    return eph_equation_of_time(jd, &e) == 0 ? e : NAN;
}

FN(jdoubleArray, nutation)(ARGS, jdouble jd) {
    (void)cls;
    double v[2];
    eph_nutation(jd, &v[0], &v[1]);
    return doubles(env, v, 2);
}

FN(jdoubleArray, eclipticToEquatorial)(ARGS, jdouble lon, jdouble lat, jdouble jd) {
    (void)cls;
    double v[2];
    eph_ecliptic_to_equatorial(lon, lat, jd, &v[0], &v[1]);
    return doubles(env, v, 2);
}

/* ---- positions ---- */

/* mode: 0 TT geometric, 1 TT apparent, 2 UT geometric, 3 UT apparent. */
FN(jdoubleArray, position)(ARGS, jint body, jdouble jd, jdouble year, jint mode) {
    (void)cls;
    eph_position p;
    int rc;
    switch (mode) {
    case 0: rc = eph_position_tt((eph_body)body, jd, &p); break;
    case 1: rc = eph_position_apparent_tt((eph_body)body, jd, &p); break;
    case 2: rc = eph_position_ut((eph_body)body, jd, year, &p); break;
    default: rc = eph_position_apparent_ut((eph_body)body, jd, year, &p); break;
    }
    return position(env, rc, &p);
}

FN(jdoubleArray, heliocentric)(ARGS, jint body, jdouble jd) {
    (void)cls;
    double v[3];
    if (eph_heliocentric_tt((eph_body)body, jd, &v[0], &v[1], &v[2]) != 0) return NULL;
    return doubles(env, v, 3);
}

FN(jstring, bodyName)(ARGS, jint body) { (void)cls; return string(env, eph_body_name((eph_body)body)); }

FN(jdoubleArray, topocentric)(ARGS, jint body, jdouble jd, jdouble year, jdouble lon, jdouble lat,
                              jdouble elevation) {
    (void)cls;
    eph_observer o = observer(lon, lat, elevation);
    eph_position p;
    return position(env, eph_position_topocentric_ut((eph_body)body, jd, year, &o, &p), &p);
}

/* ---- houses ---- */

FN(jdoubleArray, ascMc)(ARGS, jdouble jd, jdouble lon, jdouble lat) {
    (void)cls;
    eph_observer o = observer(lon, lat, 0);
    double v[2];
    eph_asc_mc(jd, &o, &v[0], &v[1]);
    return doubles(env, v, 2);
}

/* Twelve cusps and the return code (1: Porphyry fallback); null on error. */
FN(jdoubleArray, houses)(ARGS, jdouble jd, jdouble lon, jdouble lat, jint system) {
    (void)cls;
    eph_observer o = observer(lon, lat, 0);
    double v[13];
    int rc = eph_houses(jd, &o, (eph_house_system)system, v);
    if (rc < 0) return NULL;
    v[12] = rc;
    return doubles(env, v, 13);
}

FN(jdoubleArray, houseAngles)(ARGS, jdouble jd, jdouble lon, jdouble lat) {
    (void)cls;
    eph_observer o = observer(lon, lat, 0);
    eph_angles a;
    eph_house_angles(jd, &o, &a);
    double v[8] = {a.ascendant, a.mc, a.armc, a.vertex, a.equatorial_ascendant,
                   a.co_ascendant_koch, a.co_ascendant_munkasey, a.polar_ascendant};
    return doubles(env, v, 8);
}

FN(jdouble, housePosition)(ARGS, jdoubleArray cusps, jdouble lon) {
    (void)cls;
    double c[12] = {0};
    (*env)->GetDoubleArrayRegion(env, cusps, 0, 12, c);
    return eph_house_position(c, lon);
}

/* ---- phenomena ---- */

FN(jdoubleArray, phenomena)(ARGS, jint body, jdouble jd) {
    (void)cls;
    eph_phenomena p;
    if (eph_phenomena_tt((eph_body)body, jd, &p) != 0) return NULL;
    double v[5] = {p.phase_angle, p.phase, p.elongation, p.diameter, p.magnitude};
    return doubles(env, v, 5);
}

/* ---- searches (NaN: none) ---- */

#define RESULT(call) do { double out; return (call) == 0 ? out : NAN; } while (0)

FN(jdouble, longitudeCrossingUt)(ARGS, jint body, jdouble lon, jdouble jd, jint dir, jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_longitude_crossing_ut((eph_body)body, lon, jd, dir, max, &out));
}

FN(jdouble, longitudeCrossingSidereal)(ARGS, jint body, jdouble lon, jint system, jdouble jd, jint dir,
                                       jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_longitude_crossing_sidereal((eph_body)body, lon, (eph_ayanamsa_system)system, jd, dir,
                                           max, &out));
}

FN(jdouble, longitudeCrossingSiderealCustom)(ARGS, jint body, jdouble lon, jdouble jd0, jdouble value0,
                                             jdouble jd, jint dir, jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_longitude_crossing_sidereal_custom((eph_body)body, lon, jd0, value0, jd, dir, max, &out));
}

FN(jdouble, heliocentricCrossing)(ARGS, jint body, jdouble lon, jdouble jd, jint dir, jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_heliocentric_crossing((eph_body)body, lon, jd, dir, max, &out));
}

FN(jdouble, heliocentricCrossingSidereal)(ARGS, jint body, jdouble lon, jint system, jdouble jd, jint dir,
                                          jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_heliocentric_crossing_sidereal((eph_body)body, lon, (eph_ayanamsa_system)system, jd, dir,
                                              max, &out));
}

FN(jdouble, heliocentricCrossingSiderealCustom)(ARGS, jint body, jdouble lon, jdouble jd0,
                                                jdouble value0, jdouble jd, jint dir, jdouble max) {
    (void)env; (void)cls;
    RESULT(eph_heliocentric_crossing_sidereal_custom((eph_body)body, lon, jd0, value0, jd, dir, max, &out));
}

/* [jd, longitude, 1 ascending / 2 descending]; null if none. */
FN(jdoubleArray, nodeCrossing)(ARGS, jint body, jdouble jd, jint dir, jdouble max) {
    (void)cls;
    double v[3];
    int rc = eph_node_crossing((eph_body)body, jd, dir, max, &v[0], &v[1]);
    if (rc < 1) return NULL;
    v[2] = rc;
    return doubles(env, v, 3);
}

/* [jd, 1 retrograde / 2 direct]; null if none. */
FN(jdoubleArray, station)(ARGS, jint body, jdouble jd, jint dir, jdouble max) {
    (void)cls;
    double v[2];
    int rc = eph_station((eph_body)body, jd, dir, max, &v[0]);
    if (rc < 1) return NULL;
    v[1] = rc;
    return doubles(env, v, 2);
}

/* ---- orbits: 5 elements, then 8 positions of 6 (heliocentric, geocentric) ---- */

FN(jdoubleArray, orbit)(ARGS, jint body, jdouble jd, jint method) {
    (void)cls;
    eph_orbit o;
    if (eph_orbit_tt((eph_body)body, jd, (eph_orbit_method)method, &o) != 0) return NULL;
    double v[53] = {o.semi_major_axis, o.eccentricity, o.inclination, o.node, o.perihelion};
    for (int i = 0; i < 8; i++) {
        const eph_position* p = i < 4 ? &o.heliocentric[i] : &o.geocentric[i - 4];
        double* w = v + 5 + 6 * i;
        w[0] = p->longitude; w[1] = p->latitude; w[2] = p->distance;
        w[3] = p->longitude_speed; w[4] = p->latitude_speed; w[5] = p->distance_speed;
    }
    return doubles(env, v, 53);
}

/* ---- eclipses ---- */

FN(jdoubleArray, solarEclipse)(ARGS, jdouble jd, jint dir) {
    (void)cls;
    eph_solar_eclipse e;
    if (eph_solar_eclipse_next(jd, dir, &e) != 0) return NULL;
    double v[13] = {e.type, e.maximum, e.gamma, e.magnitude, e.latitude, e.longitude, e.saros,
                    e.partial_begin, e.partial_end, e.central_begin, e.central_end,
                    e.center_line_begin, e.center_line_end};
    return doubles(env, v, 13);
}

FN(jdoubleArray, solarEclipseLocal)(ARGS, jdouble jd, jdouble lon, jdouble lat, jdouble elevation, jint dir) {
    (void)cls;
    eph_observer o = observer(lon, lat, elevation);
    eph_local_solar_eclipse e;
    if (eph_solar_eclipse_local(jd, &o, dir, &e) != 0) return NULL;
    double v[10] = {e.type, e.maximum, e.first_contact, e.second_contact, e.third_contact,
                    e.fourth_contact, e.magnitude, e.obscuration, e.altitude, e.saros};
    return doubles(env, v, 10);
}

FN(jdoubleArray, solarEclipseWhere)(ARGS, jdouble jd) {
    (void)cls;
    double v[2];
    if (eph_solar_eclipse_where(jd, &v[0], &v[1]) != 0) return NULL;
    return doubles(env, v, 2);
}

FN(jdoubleArray, solarEclipseHow)(ARGS, jdouble jd, jdouble lon, jdouble lat, jdouble elevation) {
    (void)cls;
    eph_observer o = observer(lon, lat, elevation);
    eph_local_eclipse_state s;
    if (eph_solar_eclipse_how(jd, &o, &s) != 0) return NULL;
    double v[4] = {s.type, s.magnitude, s.obscuration, s.altitude};
    return doubles(env, v, 4);
}

FN(jdoubleArray, lunarEclipse)(ARGS, jdouble jd, jint dir) {
    (void)cls;
    eph_lunar_eclipse e;
    if (eph_lunar_eclipse_next(jd, dir, &e) != 0) return NULL;
    double v[12] = {e.type, e.maximum, e.gamma, e.umbral_magnitude, e.penumbral_magnitude, e.saros,
                    e.penumbral_begin, e.penumbral_end, e.partial_begin, e.partial_end,
                    e.total_begin, e.total_end};
    return doubles(env, v, 12);
}

/* ---- minor planets and stars ---- */

static int find(JNIEnv* env, jstring name, int (*fn)(const char*)) {
    const char* s = (*env)->GetStringUTFChars(env, name, NULL);
    if (!s) return -1;
    int i = fn(s);
    (*env)->ReleaseStringUTFChars(env, name, s);
    return i;
}

FN(jint, minorCount)(ARGS) { (void)env; (void)cls; return eph_minor_count(); }
FN(jint, minorFind)(ARGS, jstring name) { (void)cls; return find(env, name, eph_minor_find); }
FN(jint, minorNumber)(ARGS, jint i) { (void)env; (void)cls; return eph_minor_number(i); }
FN(jstring, minorName)(ARGS, jint i) { (void)cls; return string(env, eph_minor_name(i)); }

FN(jint, starCount)(ARGS) { (void)env; (void)cls; return eph_star_count(); }
FN(jint, starFind)(ARGS, jstring name) { (void)cls; return find(env, name, eph_star_find); }
FN(jstring, starName)(ARGS, jint i) { (void)cls; return string(env, eph_star_name(i)); }
FN(jstring, starDesignation)(ARGS, jint i) { (void)cls; return string(env, eph_star_designation(i)); }
FN(jdouble, starMagnitude)(ARGS, jint i) { (void)env; (void)cls; return eph_star_magnitude(i); }

FN(jdoubleArray, starPosition)(ARGS, jint i, jdouble jd, jdouble year) {
    (void)cls;
    eph_position p;
    return position(env, eph_star_position_apparent_ut(i, jd, year, &p), &p);
}

/* ---- ayanamsas and utilities ---- */

FN(jdouble, ayanamsa)(ARGS, jint system, jdouble jd) {
    (void)env; (void)cls;
    return eph_ayanamsa((eph_ayanamsa_system)system, jd);
}
FN(jdouble, ayanamsaCustom)(ARGS, jdouble jd0, jdouble value0, jdouble jd) {
    (void)env; (void)cls;
    return eph_ayanamsa_custom(jd0, value0, jd);
}
FN(jstring, ayanamsaName)(ARGS, jint system) {
    (void)cls;
    return string(env, eph_ayanamsa_name((eph_ayanamsa_system)system));
}

FN(jdouble, midpoint)(ARGS, jdouble a, jdouble b) { (void)env; (void)cls; return eph_midpoint(a, b); }

FN(jdoubleArray, splitLongitude)(ARGS, jdouble lon) {
    (void)cls;
    int sign, d, m;
    double s;
    eph_split_longitude(lon, &sign, &d, &m, &s);
    double v[4] = {sign, d, m, s};
    return doubles(env, v, 4);
}
