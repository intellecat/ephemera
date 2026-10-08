#include "internal.h"

namespace eph {

/*
 * Geocentric lunar position from the full ELP2000-82B theory.
 *
 * ELP returns rectangular coordinates in the J2000 ecliptic; these are
 * converted to spherical form and precessed to the mean ecliptic and
 * equinox of date (IAU 1976), matching the frame used by VSOP87D.
 */
void moon_position(double jd_tt, double& lon, double& lat, double& dist_au) {
    double x, y, z;
    elp_moon_j2000(jd_tt, x, y, z);

    const double r_km = std::sqrt(x * x + y * y + z * z);
    double l = std::atan2(y, x);
    double b = std::asin(z / r_km);
    precess_ecliptic_j2000_to_date(jd_tt, l, b);

    double dist = r_km / AU_KM;
    apply_moon_correction(jd_tt, l, b, dist);

    lon = normalize(l);
    lat = b;
    dist_au = dist;
}

}  // namespace eph
