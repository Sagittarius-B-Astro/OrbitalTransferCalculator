"""Reproduces every number quoted in the README / chat. Run from the repo root:

    python3 scripts/reproduce_results.py            # everything (~1 min)
    python3 scripts/reproduce_results.py bench      # only section 1
    python3 scripts/reproduce_results.py validate   # only section 2
    python3 scripts/reproduce_results.py unimodal   # only section 3
    python3 scripts/reproduce_results.py nodal      # only section 4
"""
import sys, time, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from planechange.lambert_izzo import lambert
from planechange.pparam import velocities_from_p, time_of_flight, p_bounds, transfer_angle
from planechange.kepler import propagate
from planechange.frames import Orbit, state_at
from planechange.min_dv import min_delta_v, delta_v, _p_from_s
from planechange.optimizers import nelder_mead

MU = 3.986004418e5


# ---------------------------------------------------------------- 1. benchmark
def bench(N=2000, repeats=5):
    r1 = np.array([7000.0, 0, 0]); r2 = np.array([-3000.0, 12000, 4000])
    n_sol = [len(lambert(r1, r2, 3000.0 + i, MU)) for i in range(0, N, 200)]
    print(f"solutions returned per Lambert call (sampled): {sorted(set(n_sol))}")
    def time_it(fn):
        best = 1e9
        for _ in range(repeats):
            t = time.perf_counter()
            for i in range(N):
                fn(i)
            best = min(best, (time.perf_counter() - t) / N)
        return best
    tl = time_it(lambda i: lambert(r1, r2, 3000.0 + i, MU))
    tp = time_it(lambda i: velocities_from_p(r1, r2, 9000.0 + i, MU))
    print(f"Izzo lambert():          {tl*1e6:7.1f} us/call (best of {repeats} x {N})")
    print(f"velocities_from_p():     {tp*1e6:7.1f} us/call")
    print(f"ratio:                   {tl/tp:7.1f}x")


# ---------------------------------------------------------------- 2. validation
def hohmann_dv(r1, r2):
    at = (r1 + r2) / 2
    v = lambda r, a: np.sqrt(MU * (2 / r - 1 / a))
    return abs(v(r1, at) - v(r1, r1)) + abs(v(r2, r2) - v(r2, at))


def validate():
    print("-- (a) Blanco 2025 worked example (mu = 3.986e5, r1 = 9000, r2 = 15000, 120 deg)")
    mu = 3.986e5
    r1 = np.array([9000.0, 0, 0]); th = np.radians(120); r2 = 15000 * np.array([np.cos(th), np.sin(th), 0])
    lo, hi, par = p_bounds(9000, 15000, th)
    print(f"   p_par,1/r1 = {lo/9000:.4f} (paper 0.6317)   p_par,2/r1 = {par/9000:.4f} (paper 1.8173)")
    vc = np.sqrt(mu / 9000); p = 1.3128 * 9000
    v1, v2, d = velocities_from_p(r1, r2, p, mu)
    print(f"   dv/vc = {np.linalg.norm(v1-[0,vc,0])/vc:.4f} (paper 0.1563)   "
          f"dt/Tc = {time_of_flight(9000,15000,d,p,v1,r1,mu)/(2*np.pi*9000/vc):.4f} (paper 0.4607)")

    print("-- (b) p-family vs Izzo (max |dv1|, |dv2| difference over p values)")
    for deg, ps in ((120, [7000, 9000, 11816, 16000, 25000]), (240, [2000, 5000, 9000, 16000])):
        t = np.radians(deg); r2 = 15000 * np.array([np.cos(t), np.sin(t), 0]); worst = 0
        for p in ps:
            v1, v2, d = velocities_from_p(r1, r2, p, mu)
            tof = time_of_flight(9000, 15000, d, p, v1, r1, mu)
            s = min(lambert(r1, r2, tof, mu), key=lambda s: np.linalg.norm(s["v1"] - v1))
            worst = max(worst, np.linalg.norm(s["v1"] - v1), np.linalg.norm(s["v2"] - v2))
        print(f"   dtheta = {deg} deg: worst velocity mismatch {worst:.2e} km/s")

    print("-- (c) Izzo solutions propagated with the independent Kepler propagator")
    worst = 0; count = 0
    r1 = np.array([9000.0, 0, 0])
    for deg in (60, 120, 175, 240, 300):
        t = np.radians(deg); r2 = 15000 * np.array([np.cos(t), np.sin(t), 0])
        for tof in (1000.0, 3915.0, 8000.0, 30000.0, 60000.0):
            for s in lambert(r1, r2, tof, MU):
                r, _ = propagate(r1, s["v1"], tof, MU)
                worst = max(worst, np.linalg.norm(r - r2)); count += 1
    print(f"   {count} solutions (M = 0,1,2,...), worst final-position error {worst:.2e} km")

    print("-- (d) free-time search vs closed forms (n_grid = 16)")
    r = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000), MU, n_grid=16)
    print(f"   coplanar 7000->14000: search {r['dv']:.5f}  Hohmann {hohmann_dv(7000,14000):.5f} km/s  a={r['a']:.1f} e={r['e']:.4f}")
    di = np.radians(20)
    r = min_delta_v(Orbit(7000, 7000), Orbit(7000, 7000, di), MU, n_grid=16)
    print(f"   pure 20 deg plane change: search {r['dv']:.5f}  2v sin(di/2) {2*np.sqrt(MU/7000)*np.sin(di/2):.5f} km/s")
    r = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000, di), MU, n_grid=16)
    at = 10500
    vp, va = (np.sqrt(MU * (2 / x - 1 / at)) for x in (7000, 14000))
    v1c, v2c = np.sqrt(MU / 7000), np.sqrt(MU / 14000)
    naive = abs(vp - v1c) + np.sqrt(va**2 + v2c**2 - 2 * va * v2c * np.cos(di))
    print(f"   7000->14000 + 20 deg: search {r['dv']:.5f}  all-plane-change-at-apogee {naive:.5f}  saving {100*(naive-r['dv'])/naive:.1f} %")

    print("-- (e) dt monotonic in p (40 values; sign of successive differences)")
    mu = 3.986e5; r1 = np.array([9000.0, 0, 0])
    for deg, rng_ in ((120, (5800, 40000)), (240, (500, 16000))):
        t = np.radians(deg); r2 = 15000 * np.array([np.cos(t), np.sin(t), 0]); ts = []
        for p in np.linspace(*rng_, 40):
            v1, _, d = velocities_from_p(r1, r2, p, mu); ts.append(time_of_flight(9000, 15000, d, p, v1, r1, mu))
        sg = np.unique(np.sign(np.diff(ts)))
        print(f"   dtheta = {deg} deg: unique sign(diff dt) = {sg}  (-1 = dt falls as p rises)")


# ---------------------------------------------------------------- 3. unimodality
def unimodal(n_geoms=300, seed=1):
    rng = np.random.default_rng(seed); bad = n = 0
    for _ in range(n_geoms):
        o1 = Orbit(rng.uniform(7000, 9000), rng.uniform(6800, 7000), rng.uniform(0, 1), rng.uniform(0, 6), rng.uniform(0, 6))
        o2 = Orbit(rng.uniform(10000, 30000), rng.uniform(9000, 10000), rng.uniform(0, 1), rng.uniform(0, 6), rng.uniform(0, 6))
        n1, n2 = rng.uniform(0, 6.28, 2)
        a, _ = state_at(o1, n1, MU); b, _ = state_at(o2, n2, MU)
        d = transfer_angle(a, b, o1.normal())
        if abs(np.sin(d)) < 1e-3:
            continue
        lo, hi, _ = p_bounds(np.linalg.norm(a), np.linalg.norm(b), d)
        ss = np.linspace(-10, 10, 400)
        v = np.array([delta_v(o1, o2, n1, n2, _p_from_s(s, lo, hi), MU) for s in ss])
        n += 1; bad += int(np.sum((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:])) > 1)
    print(f"geometries tested: {n}; with >1 interior local minimum of dv(p) on the 400-point sweep: {bad}")


# ---------------------------------------------------------------- 4. nodal parametrisation
def nodal_dv(r1v, r2v, psi, q, mu):
    """Two-impulse dv for ANTIPARALLEL r1, r2 (dtheta = pi). Free parameters: psi, q = e sin(nu1)."""
    r1n, r2n = np.linalg.norm(r1v), np.linalg.norm(r2v)
    return r1n, r2n


def nodal_velocities(r1v, r2v, psi, q, mu):
    r1n, r2n = np.linalg.norm(r1v), np.linalg.norm(r2v)
    rh = r1v / r1n
    assert np.allclose(r2v / r2n, -rh, atol=1e-9), "points are not antiparallel"
    p = 2 * r1n * r2n / (r1n + r2n)
    k = (r2n - r1n) / (r1n + r2n)               # e cos(nu1)
    u = np.cross(rh, [0, 0, 1.0]); u = u / np.linalg.norm(u) if np.linalg.norm(u) > 1e-9 else np.array([1.0, 0, 0])
    n = np.cos(psi) * u + np.sin(psi) * np.cross(rh, u)   # plane normal, rotated about the nodal line
    t1 = np.cross(n, rh)
    s = np.sqrt(mu / p)
    return s * (q * rh + (1 + k) * t1), s * (q * rh - (1 - k) * t1)


def nodal_min(o1, o2, mu):
    d = np.cross(o1.normal(), o2.normal()); d /= np.linalg.norm(d)
    best = (1e9, None)
    for sgn in (+1, -1):                                  # burn 1 at +d / -d, burn 2 at the opposite node
        nus = []
        for o, direction in ((o1, sgn * d), (o2, -sgn * d)):
            loc = o.dcm.T @ direction
            nus.append(np.arctan2(loc[1], loc[0]))
        ra, vc1 = state_at(o1, nus[0], mu); rb, vc2 = state_at(o2, nus[1], mu)
        def f(x):
            v1, v2 = nodal_velocities(ra, rb, x[0], x[1], mu)
            return np.linalg.norm(v1 - vc1) + np.linalg.norm(vc2 - v2)
        for psi0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            x, fx, _ = nelder_mead(f, [[psi0, 0.0], [psi0 + 0.3, 0.0], [psi0, 0.3]], xtol=1e-11, ftol=1e-14)
            if fx < best[0]:
                best = (fx, (sgn, np.degrees(nus[0]) % 360, np.degrees(nus[1]) % 360, x[0] % (2 * np.pi), x[1]))
    return best


def nodal():
    cases = {
        "circ 7000 -> circ 14000, di = 20 deg": (Orbit(7000, 7000), Orbit(14000, 14000, np.radians(20))),
        "circ 7000 -> circ 14000, di = 20 deg, RAAN2 = 30 deg": (Orbit(7000, 7000), Orbit(14000, 14000, np.radians(20), np.radians(30))),
        "elliptic 7000x9000 -> elliptic 12000x20000, tilted": (Orbit(9000, 7000, 0.2, 0.3, 0.4), Orbit(20000, 12000, 0.7, 1.1, 2.0)),
    }
    for name, (o1, o2) in cases.items():
        fx, info = nodal_min(o1, o2, MU)
        g = min_delta_v(o1, o2, MU, n_grid=20)["dv"]
        print(f"{name}\n   nodal 2-parameter search: {fx:.6f} km/s   general 3-variable search: {g:.6f} km/s   diff {fx-g:+.2e}")
        print(f"   (burn1 nu1 = {info[1]:.2f} deg, burn2 nu2 = {info[2]:.2f} deg, psi = {np.degrees(info[3]):.2f} deg, q = e sin(nu1) = {info[4]:+.4f})")


if __name__ == "__main__":
    sections = {"bench": bench, "validate": validate, "unimodal": unimodal, "nodal": nodal}
    which = sys.argv[1:] or list(sections)
    for w in which:
        print(f"\n===== {w} =====")
        sections[w]()
