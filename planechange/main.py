"""Entry point used by the web app (via Pyodide) and by the tests.

`solve(params)` takes plain numbers (km, deg, km^3/s^2) so it can be called from
JavaScript with a single JSON object.
"""
import numpy as np
from .frames import Orbit, state_at
from .kepler import propagate
from .min_dv import min_delta_v
from .nodal import min_delta_v_nodal

def conic_for_plot(r1, v1, r2, mu):
    h = np.cross(r1, v1); hn = np.linalg.norm(h); hh = h / hn
    evec = np.cross(v1, h) / mu - r1 / np.linalg.norm(r1)
    e = np.linalg.norm(evec)
    P = evec / e if e > 1e-9 else r1 / np.linalg.norm(r1)   # circular: periapsis direction arbitrary
    Q = np.cross(hh, P)
    nu1 = np.arctan2(r1 @ Q, r1 @ P)
    dth = np.arctan2(np.cross(r1, r2) @ hh, r1 @ r2) % (2 * np.pi)
    return dict(p=float(hn**2 / mu), e=float(e), P=P.tolist(), Q=Q.tolist(),
                nuStart=float(nu1), nuEnd=float(nu1 + dth))

def min_radius_on_arc(c):
    two_pi = 2 * np.pi
    k = np.ceil(c["nuStart"] / two_pi)                 # first periapsis passage at/after the start
    if k * two_pi <= c["nuEnd"]:
        return c["p"] / (1 + c["e"])
    r = lambda nu: c["p"] / (1 + c["e"] * np.cos(nu))
    return min(r(c["nuStart"]), r(c["nuEnd"]))

def solve(params):
    d = np.radians
    o1 = Orbit(params["r1a"], params["r1p"], d(params["i1"]), d(params["RAAN1"]), d(params["w1"]))
    o2 = Orbit(params["r2a"], params["r2p"], d(params["i2"]), d(params["RAAN2"]), d(params["w2"]))
    mu = params["mu"]
    candidates = []
    if np.linalg.norm(np.cross(o1.normal(), o2.normal())) > 1e-6:   # planes differ -> line of nodes exists
        nod = min_delta_v_nodal(o1, o2, mu); nod["method"] = "nodal"; candidates.append(nod)
    gen = min_delta_v(o1, o2, mu, n_grid=int(params.get("n_grid", 36))); gen["method"] = "general"
    candidates.append(gen)
    res = min(candidates, key=lambda c: c["dv"])
    a = res.get("a")
    if a is None:
        a = res["p"] / (1 - res["e"] ** 2) if abs(res["e"] - 1) > 1e-9 else None
    r1, _ = state_at(o1, res["nu1"], mu)
    r2, _ = state_at(o2, res["nu2"], mu)
    conic = conic_for_plot(r1, res["v1"], r2, mu)
    return {
        "totalDeltaV": float(res["dv"]), "transferTime": float(res["tof"]),
        "nu1Deg": float(np.degrees(res["nu1"])), "nu2Deg": float(np.degrees(res["nu2"])),
        "p": float(res["p"]), "a": None if a is None else float(a), "e": float(res["e"]),
        "conic": conic, "minRadius": min_radius_on_arc(conic), "method": res["method"],
    }
