#include "internal.h"

#include "data.h"

namespace eph {

/*
 * Pluto heliocentric position in the J2000 ecliptic frame, from Chebyshev
 * intervals fitted to JPL DE440 (tools/fit_pluto.py; 1550-2650, ~0.03").
 * Returns false outside the fitted range or when no data is loaded.
 */
bool pluto_heliocentric_j2000(double jd_tt, double& lon, double& lat, double& r) {
    const Data& d = data();
    Vec3 v;
    if (!d.pluto_ok || !chebyshev_body(d.pluto.t, d.pluto.c, d.pluto.nseg,
                                       d.pluto.ncoef, jd_tt, v))
        return false;
    cartesian_to_spherical(v, lon, lat, r);
    return true;
}

}  // namespace eph
