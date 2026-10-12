#!/usr/bin/env python3
"""Build a minor-planet data file (development only): data/eph-minor*.bin.

For each body, from public JPL data:

  1. Osculating elements from the JPL Small-Body Database (full precision).
  2. A numerical integration over the span (1900-2100 by default) with the
     Sun, the planets, the Moon and Pluto from DE440, the 16 massive
     asteroids JPL's own ephemerides use (integrated with each other first,
     masses from DE440), and the Sun's relativistic term.
  3. A check against JPL Horizons at five dates.  A body that misses the
     tolerance (near-Earth objects with close approaches) is integrated
     again in 10-year stretches, each started from a Horizons state, and
     checked again.
  4. Chebyshev pieces of the heliocentric ecliptic-J2000 position, split
     until each meets the tolerance: 0.5" as seen from the closest the
     Earth can get (q - 1 AU, at least 0.3 AU).  The constant terms are
     double, the others float32.

The section format is in src/data.h (kSectionMinor).  Needs numpy, scipy,
jplephem (Skyfield) and network access; DE440 is read from $DE440.

    DE440=de440.bsp python3 tools/make_minor_planets.py data/eph-minor.bin \\
        --numbers 1-1000 --extra astrology
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import re
import struct
import sys
import time
import urllib.parse
import urllib.request
import zlib

import numpy as np
from numpy.polynomial import chebyshev as Ch
from jplephem.spk import SPK
from scipy.integrate import solve_ivp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AU, DAY = 149597870.7, 86400.0
C_LIGHT = 299792.458 * DAY / AU
EPS0 = math.radians(84381.448 / 3600)  # ecliptic J2000 of JPL (IAU 1976)
ECL_TO_EQ = np.array([[1, 0, 0], [0, math.cos(EPS0), -math.sin(EPS0)],
                      [0, math.sin(EPS0), math.cos(EPS0)]])
SECTION_MINOR = 9
# The massive asteroids of JPL's ephemerides (sb441-n16), with masses in DE440.
BIG16 = [1, 2, 3, 4, 7, 10, 15, 16, 31, 52, 65, 87, 88, 107, 511, 704]

# Bodies beyond the first numbered ones that astrologers use: near-Earth
# objects, centaurs, trans-Neptunian objects and named favourites.
ASTROLOGY = [
    # near-Earth and Earth-crossing
    433, 1036, 1181, 1221, 1566, 1620, 1627, 1685, 1862, 1864, 1866, 1943, 1981,
    2061, 2063, 2100, 2101, 2212, 2340, 3122, 3200, 3362, 3753, 4179, 4660, 4769,
    5335, 6489, 25143, 69230, 99942, 101955, 162173,
    # named favourites
    1034, 1143, 1172, 1173, 1387, 1388, 1437, 1583, 1809, 2207, 3317, 3811, 4341,
    4450, 5011, 5261, 5731, 7088,
    # centaurs
    2060, 5145, 7066, 8405, 10199, 10370, 31824, 32532, 49036, 52872, 54598, 55576,
    60558, 83982,
    # trans-Neptunian objects
    15760, 19521, 20000, 28978, 38628, 42355, 47171, 50000, 58534, 65489, 66652,
    90377, 90482, 120347, 136108, 136199, 136472, 174567, 225088,
]

GM = {}
kernel = None
segments = {}


def load_de440(path):
    global kernel
    kernel = SPK.open(path)
    for s in kernel.segments:
        segments[(s.center, s.target)] = s
    text = fetch("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de440.tpc")
    for a, b in re.findall(r"BODY(\d+)_GM\s*=\s*\(\s*([0-9.D+E-]+)", text):
        GM[int(a)] = float(b.replace("D", "E")) * DAY**2 / AU**3  # AU^3/day^2


def fetch(url, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:  # rate limits, timeouts: back off
            if k == tries - 1:
                raise
            sys.stderr.write("  retry (%s)\n" % exc)
            time.sleep(10 * (k + 1))


def sun_state(jd):
    p, v = segments[(0, 10)].compute_and_differentiate(jd)
    return p / AU, v / AU


def planets(jd):
    """Barycentric positions (AU) and GMs of the Sun, planets, Moon, Pluto."""
    emb = segments[(0, 3)].compute(jd) / AU
    out = [(segments[(0, n)].compute(jd) / AU, GM[n]) for n in (10, 1, 2, 4, 5, 6, 7, 8, 9)]
    out.append((emb + segments[(3, 399)].compute(jd) / AU, GM[399]))
    out.append((emb + segments[(3, 301)].compute(jd) / AU, GM[301]))
    return out


def acceleration(t, r, v, extra):
    """Barycentric acceleration of particles at r (n x 3); `extra` gives
    more point masses [(position, GM)] at t."""
    acc = np.zeros_like(r)
    for p, gm in planets(t) + extra:
        d = p - r
        dist = np.linalg.norm(d, axis=1)
        dist[dist < 1e-9] = np.inf  # a body does not attract itself
        acc += gm * d / (dist**3)[:, None]
    ps, vs = sun_state(t)
    rh, vh = r - ps, v - vs
    rr, mu = np.linalg.norm(rh, axis=1), GM[10]
    acc += mu / (C_LIGHT**2 * rr**3)[:, None] * (
        (4 * mu / rr - np.sum(vh * vh, axis=1))[:, None] * rh
        + 4 * np.sum(rh * vh, axis=1)[:, None] * vh)
    return acc


def kepler_state(a, e, i, om, w, ma, mu):
    """Heliocentric ecliptic-J2000 state from osculating elements."""
    i, om, w, M = map(math.radians, (i, om, w, ma))
    if e < 1:
        E = M
        for _ in range(100):
            E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
        n = math.sqrt(mu / a**3)
        f = n / (1 - e * math.cos(E))
        vx, vy = -a * math.sin(E) * f, a * math.sqrt(1 - e * e) * math.cos(E) * f
    else:
        raise ValueError("unbound orbit")
    co, so = math.cos(om), math.sin(om)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(w), math.sin(w)
    P = np.array([co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si])
    Q = np.array([-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si])
    return x * P + y * Q, vx * P + vy * Q


def integrate(states, t_from, spans, extra_fn):
    """Integrate massless particles (barycentric equatorial states) from
    t_from to each end in `spans`; dense solutions."""
    n = len(states)
    y0 = np.concatenate([np.ravel([s[0] for s in states]), np.ravel([s[1] for s in states])])

    def rhs(t, y):
        r, v = y[:3 * n].reshape(n, 3), y[3 * n:].reshape(n, 3)
        return np.concatenate([v.ravel(), acceleration(t, r, v, extra_fn(t)).ravel()])

    out = []
    for t_end in spans:
        if abs(t_end - t_from) < 1e-9:
            out.append(None)
            continue
        out.append(solve_ivp(rhs, (t_from, t_end), y0, method="DOP853", rtol=1e-12,
                             atol=1e-15, dense_output=True))
    return out


class Track:
    """A body's heliocentric ecliptic-J2000 position from dense
    integrations: [(t0, t1, solution, index)]."""

    def __init__(self):
        self.parts = []

    def add(self, t0, t1, sol, k, n=None):
        self.parts.append((min(t0, t1), max(t0, t1), sol, k))

    def __call__(self, t):
        for a, b, sol, k in self.parts:
            if a - 1e-9 <= t <= b + 1e-9:
                y = sol.sol(t)
                return ECL_TO_EQ.T @ (y[3 * k:3 * k + 3] - sun_state(t)[0])
        raise ValueError("no track at %.1f" % t)


class BlendedTrack:
    """Integrations started from states at anchors ta[j], each covering
    ta[j] +- step, blended with smooth weights (a partition of unity of
    cos^2 ramps between neighbouring anchors) so the result has no jumps
    where one integration hands over to the next."""

    def __init__(self, anchors, step):
        self.anchors, self.step, self.tracks = anchors, step, []

    def __call__(self, t):
        j = int(np.clip(np.searchsorted(self.anchors, t) - 1, 0, len(self.anchors) - 2))
        a, b = self.anchors[j], self.anchors[j + 1]
        x = min(max((t - a) / (b - a), 0.0), 1.0)
        w = math.cos(x * math.pi / 2) ** 2
        if w > 1 - 1e-15:
            return self.tracks[j](t)
        if w < 1e-15:
            return self.tracks[j + 1](t)
        return w * self.tracks[j](t) + (1 - w) * self.tracks[j + 1](t)


def sbdb(numbers):
    """Full-precision elements, H, G, diameter for numbered bodies."""
    out = {}
    for i in range(0, len(numbers), 40):  # longer queries get HTTP 502
        chunk = numbers[i:i + 40]
        q = urllib.parse.urlencode({
            "fields": "pdes,name,epoch,e,a,q,i,om,w,ma,H,G,diameter", "full-prec": "true",
            "sb-cdata": json.dumps({"OR": ["pdes|EQ|%d" % n for n in chunk]})})
        d = json.loads(fetch("https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + q))
        for row in d["data"]:
            rec = dict(zip(d["fields"], row))
            out[int(rec["pdes"])] = rec
    return out


def horizons_vectors(num, jds, epoch):
    """Heliocentric ecliptic-J2000 positions and velocities (AU, AU/day)
    from JPL Horizons.  Horizons integrates a small body from its orbit
    solution's epoch through the requested dates in time order, so a request
    mixing dates before and after the epoch can cross a close approach the
    other way round; for chaotic bodies (Apophis after 2029) that changes
    the positions by up to 0.1 AU.  The dates before and after the epoch are
    therefore requested separately, each outward from it (25 per request)."""
    jds = list(jds)
    out_p, out_v = np.zeros((len(jds), 3)), np.zeros((len(jds), 3))
    before = sorted((j for j in range(len(jds)) if jds[j] < epoch), key=lambda j: -jds[j])
    after = sorted((j for j in range(len(jds)) if jds[j] >= epoch), key=lambda j: jds[j])
    for side in (before, after):
        for i in range(0, len(side), 25):
            idx = side[i:i + 25]
            p, v = horizons_request(num, sorted(jds[j] for j in idx))
            order = sorted(idx, key=lambda j: jds[j])  # Horizons answers in time order
            out_p[order], out_v[order] = p, v
    return out_p, out_v


def horizons_request(num, jds):
    q = urllib.parse.urlencode({
        "format": "text", "COMMAND": "'%d;'" % num, "CENTER": "'500@10'",
        "EPHEM_TYPE": "VECTORS", "REF_PLANE": "ECLIPTIC", "REF_SYSTEM": "J2000",
        "VEC_TABLE": "2", "OUT_UNITS": "AU-D", "CSV_FORMAT": "YES",
        "TLIST": " ".join("'%.6f'" % j for j in jds)})
    t = fetch("https://ssd.jpl.nasa.gov/api/horizons.api?" + q)
    if "$$SOE" not in t:  # e.g. dates outside the body's ephemeris
        lines = [ln.strip() for ln in t.splitlines() if ln.strip()]
        raise RuntimeError("Horizons: " + " | ".join(lines[-3:])[:300])
    rows = t[t.index("$$SOE") + 5:t.index("$$EOE")].strip().splitlines()
    vals = np.array([[float(x) for x in r.split(",")[2:8]] for r in rows])
    return vals[:, :3], vals[:, 3:]


def tolerance(q):
    """0.5" at the closest the Earth can get to the body (AU)."""
    return 0.5 / 206264.8 * max(q - 1.0, 0.3)


def fit_body(track, t0, t1, tol, degrees=(8, 10, 12, 14, 16, 20, 24)):
    """Chebyshev pieces meeting tol; the smallest over the degrees tried:
    (degree, boundaries, constant terms, float32 tails)."""
    best = None
    for D in degrees:
        nodes = np.cos(np.pi * (np.arange(D + 1) + 0.5) / (D + 1))
        check = np.linspace(-1, 1, 3 * D + 1)
        bounds, c0s, tails = [], [], []

        def fit(a, b, depth=0):
            ts = a + (nodes + 1) * (b - a) / 2
            P = np.array([track(t) for t in ts])
            coef = [Ch.chebfit(nodes, P[:, c], D) for c in range(3)]
            for c in range(3):
                coef[c][1:] = coef[c][1:].astype(np.float32)
            Q = np.array([track(a + (u + 1) * (b - a) / 2) for u in check])
            err = np.max(np.linalg.norm(
                np.array([Ch.chebval(check, coef[c]) for c in range(3)]).T - Q, axis=1))
            if err > tol and depth < 30 and b - a > 0.05:
                m = (a + b) / 2
                fit(a, m, depth + 1)
                fit(m, b, depth + 1)
                return
            bounds.append(a)
            c0s.append([coef[c][0] for c in range(3)])
            tails.append(np.array([coef[c][1:] for c in range(3)], dtype=np.float32))

        # Start from pieces of ~1/4 orbit (at least 50 days, at most 3200).
        L0 = 3200.0
        n0 = max(1, int(math.ceil((t1 - t0) / L0)))
        for k in range(n0):
            fit(t0 + k * (t1 - t0) / n0, t0 + (k + 1) * (t1 - t0) / n0)
        size = len(bounds) * (8 + 24 + 12 * D)
        if best is None or size < best[0]:
            best = (size, D, bounds + [t1], c0s, tails)
    return best[1:]


def pack(bodies):
    """The minor-planet section (format in src/data.h)."""
    pool, offsets = b"", []
    for b in bodies:
        offsets.append(len(pool))
        pool += b["name"].encode() + b"\0"
    pool += b"\0" * (-len(pool) % 8)
    head = 8 + 48 * len(bodies) + len(pool)
    blobs, recs, off = b"", b"", head
    for b, name_off in zip(bodies, offsets):
        D, bounds, c0s, tails = b["fit"]
        nseg = len(c0s)
        blob = (np.array(bounds, dtype="<f8").tobytes() + np.array(c0s, dtype="<f8").tobytes()
                + np.array(tails, dtype="<f4").tobytes())
        blob += b"\0" * (-len(blob) % 8)
        recs += struct.pack("<IIffffIIQQ", b["number"], name_off, b["H"], b["G"],
                            b["radius"], b["q"], nseg, D + 1, off + len(blobs), 0)
        blobs += blob
    return struct.pack("<II", len(bodies), len(pool)) + recs + pool + blobs


def write_file(path, payload, info):
    sections = [(7, 0, info.encode()), (SECTION_MINOR, 0, payload)]
    off = 16 + 24 * len(sections)
    index, body = b"", b""
    for t, sid, data in sections:
        aligned = (off + len(body) + 7) // 8 * 8
        body += b"\0" * (aligned - off - len(body))
        index += struct.pack("<IIQQ", t, sid, aligned, len(data))
        body += data
    rest = index + body
    with open(path, "wb") as f:
        f.write(b"EPHD" + struct.pack("<III", 2, len(sections), zlib.crc32(rest)) + rest)


def parse_numbers(spec):
    out = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += range(int(a), int(b) + 1)
        elif part:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--numbers", default="")
    ap.add_argument("--extra", default="", help="'astrology' or a list of numbers")
    ap.add_argument("--start", type=float, default=2415020.5)  # 1900-01-01
    ap.add_argument("--end", type=float, default=2488070.5)    # 2100-01-02
    ap.add_argument("--report", default=None, help="write the check results here")
    ap.add_argument("--cache", default=os.path.join(ROOT, "build", "minor-cache"),
                    help="finished bodies (a run can be resumed)")
    args = ap.parse_args()
    load_de440(os.environ.get("DE440", os.path.join(ROOT, "de440.bsp")))
    numbers = parse_numbers(args.numbers)
    numbers += ASTROLOGY if args.extra == "astrology" else parse_numbers(args.extra)
    numbers = sorted(set(numbers))
    t0, t1 = args.start, args.end
    log = sys.stderr.write

    els = sbdb(sorted(set(numbers + BIG16)))
    missing = [n for n in numbers if n not in els]
    if missing:
        log("not in SBDB: %s\n" % missing)
    numbers = [n for n in numbers if n in els]

    def state(n):
        e = els[n]
        r, v = kepler_state(float(e["a"]), float(e["e"]), float(e["i"]), float(e["om"]),
                            float(e["w"]), float(e["ma"]), GM[10])
        ps, vs = sun_state(float(e["epoch"]))
        return ECL_TO_EQ @ r + ps, ECL_TO_EQ @ v + vs

    # 1. The massive asteroids, attracting each other.
    epoch = float(els[1]["epoch"])
    assert all(float(els[n]["epoch"]) == epoch for n in BIG16), "BIG16 epochs differ"
    big_states = [state(n) for n in BIG16]
    big_gm = [GM[2000000 + n] for n in BIG16]
    nb = len(BIG16)

    def rhs_big(t, y):
        r, v = y[:3 * nb].reshape(nb, 3), y[3 * nb:].reshape(nb, 3)
        extra = [(r[k], big_gm[k]) for k in range(nb)]
        return np.concatenate([v.ravel(), acceleration(t, r, v, extra).ravel()])

    y0 = np.concatenate([np.ravel([s[0] for s in big_states]),
                         np.ravel([s[1] for s in big_states])])
    t_start = time.time()
    # Wide enough for the re-anchored integrations (+- 10 years beyond).
    big_lo, big_hi = t0 - 30 - 3653, t1 + 30 + 3653
    big_sols = [solve_ivp(rhs_big, (epoch, te), y0, method="DOP853", rtol=1e-12,
                          atol=1e-15, dense_output=True) for te in (big_lo, big_hi)]
    log("massive asteroids integrated (%.0f s)\n" % (time.time() - t_start))

    def big_at(t):
        if not big_lo <= t <= big_hi:  # never extrapolate the dense output
            raise ValueError("massive asteroids not integrated at %.1f" % t)
        sol = big_sols[0] if t < epoch else big_sols[1]
        y = sol.sol(t)
        return [(y[3 * k:3 * k + 3], big_gm[k]) for k in range(nb)]

    check_jds = [t0 + (t1 - t0) * f for f in (0.013, 0.262, 0.513, 0.737, 0.987)]
    os.makedirs(args.cache, exist_ok=True)

    def finish(n, track):
        """Check against Horizons (dates between any re-anchoring points;
        re-anchor if needed), fit, and cache the body."""
        e = els[n]
        q = float(e["q"])
        tol = tolerance(q)
        steps = (3652.5, 730.5, 182.625)  # re-anchoring intervals, coarse to fine
        for attempt in range(len(steps) + 1):
            ref, _ = horizons_vectors(n, check_jds, float(e["epoch"]))
            err = max(np.linalg.norm(track(t) - r) for t, r in zip(check_jds, ref))
            if err <= tol / 2 or attempt == len(steps):
                break
            # Re-anchor: Horizons states every `step` days, each integrated
            # over +-step and blended with its neighbours.
            step = steps[attempt]
            # Anchors inside the span (some bodies' Horizons ephemerides
            # start on 1900-01-02); the end ones are integrated beyond it.
            k = int(math.ceil((t1 - t0 - 4) / step)) + 1
            anchors = list(np.linspace(t0 + 2, t1 - 2, k))
            pos, vel = horizons_vectors(n, anchors, float(e["epoch"]))
            track = BlendedTrack(anchors, step)
            for j, ta in enumerate(anchors):
                ps, vs = sun_state(ta)
                st = (ECL_TO_EQ @ pos[j] + ps, ECL_TO_EQ @ vel[j] + vs)
                # Blending only needs each integration up to its neighbours.
                lo = anchors[j - 1] if j > 0 else t0 - 30
                hi = anchors[j + 1] if j < k - 1 else t1 + 30
                part = Track()
                for sol, te in zip(integrate([st], ta, (lo, hi), big_at), (lo, hi)):
                    if sol is not None:
                        part.add(ta, te, sol, 0)
                track.tracks.append(part)
            log("  %d %s: %.1e AU > %.1e, re-anchored every %.0f d\n"
                % (n, e["name"] or "", err, tol, step))
        if err > tol:
            log("  WARNING %d: %.1e AU after re-anchoring (tol %.1e)\n" % (n, err, tol))
        fit = fit_body(track, t0, t1, tol)
        diam = float(e["diameter"]) if e.get("diameter") else 0.0
        body = {"number": n, "name": e["name"] or e["pdes"], "q": q,
                "H": float(e["H"]) if e.get("H") else float("nan"),
                "G": float(e["G"]) if e.get("G") else 0.15,
                "radius": diam / 2, "fit": fit, "check": (err, tol)}
        with open(os.path.join(args.cache, "%d.pkl" % n), "wb") as f:
            pickle.dump(body, f)
        log("%7d %-14s %4d pieces, degree %2d, check %.1e AU (tol %.1e)\n"
            % (n, body["name"], len(fit[2]), fit[0], err, tol))

    skipped = []

    def finish_or_skip(n, track):
        try:
            finish(n, track)
        except Exception as exc:  # one body must not stop the run
            skipped.append(n)
            log("  SKIPPED %d %s: %s\n" % (n, els[n]["name"] or "", exc))

    def cached(n):
        return os.path.exists(os.path.join(args.cache, "%d.pkl" % n))

    for k, n in enumerate(BIG16):
        if n in numbers and not cached(n):
            track = Track()
            track.add(epoch, big_lo, big_sols[0], k)
            track.add(epoch, big_hi, big_sols[1], k)
            finish_or_skip(n, track)

    # 2. The others by epoch group (massless), near-Earth objects apart (they
    #    would force small steps on everyone), in chunks small enough for the
    #    dense output to fit in memory; each chunk is checked and fitted
    #    before the next (and cached, so a run can be resumed).
    rest = [n for n in numbers if n not in BIG16 and not cached(n)]
    groups = {}
    for n in rest:
        neo = float(els[n]["q"]) < 1.3
        groups.setdefault((float(els[n]["epoch"]), neo), []).append(n)
    for (ep, neo), members in sorted(groups.items()):
        for i in range(0, len(members), 40):
            chunk = members[i:i + 40]
            t_start = time.time()
            sols = integrate([state(n) for n in chunk], ep, (t0 - 30, t1 + 30), big_at)
            log("integrated %d bodies, epoch %.1f%s (%.0f s)\n"
                % (len(chunk), ep, " (NEO)" if neo else "", time.time() - t_start))
            for k, n in enumerate(chunk):
                track = Track()
                for sol, te in zip(sols, (t0 - 30, t1 + 30)):
                    if sol is not None:
                        track.add(ep, te, sol, k)
                finish_or_skip(n, track)
            del sols

    # 3. Pack.
    bodies = []
    for n in numbers:
        if not cached(n):
            continue
        with open(os.path.join(args.cache, "%d.pkl" % n), "rb") as f:
            bodies.append(pickle.load(f))
    if skipped:
        log("skipped %d bodies: %s\n" % (len(skipped), skipped))
    report = [(b["number"], b["name"], b["q"]) + b["check"] for b in bodies]
    payload = pack(bodies)
    info = ("Minor planets: %d bodies, JD %.1f-%.1f.  Integrated from JPL SBDB elements with "
            "DE440 and checked against JPL Horizons (tools/make_minor_planets.py)."
            % (len(bodies), t0, t1))
    write_file(args.out, payload, info)
    log("wrote %s: %d bodies, %.0f KB\n" % (args.out, len(bodies), os.path.getsize(args.out) / 1e3))
    if args.report:
        with open(args.report, "w") as f:
            for n, name, q, err, tol in report:
                f.write("%d\t%s\t%.4f\t%.3e\t%.3e\n" % (n, name, q, err, tol))


if __name__ == "__main__":
    main()
