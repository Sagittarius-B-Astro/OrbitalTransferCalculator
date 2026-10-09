"""Entry point used by the web app (via Pyodide) and by the tests.

`solve(params)` takes plain numbers (km, deg, km^3/s^2) so it can be called from
JavaScript with a single JSON object.
"""
import numpy as np
from .frames import Orbit, state_at
from .kepler import propagate
from .min_dv import min_delta_v
from .nodal import min_delta_v_nodal


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
    r1, _ = state_at(o1, res["nu1"], params["mu"])
    n_arc = 80
    arc = []
    for k in range(n_arc + 1):
        r, _ = propagate(r1, res["v1"], res["tof"] * k / n_arc, params["mu"])
        arc.append([float(x) for x in r])
    return {
        "totalDeltaV": float(res["dv"]), "transferTime": float(res["tof"]),
        "nu1Deg": float(np.degrees(res["nu1"])), "nu2Deg": float(np.degrees(res["nu2"])),
        "p": float(res["p"]), "a": None if a is None else float(a), "e": float(res["e"]),
        "arc": arc, "transferPeriapsis": float(res["p"] / (1 + res["e"])), "method": res["method"],
    }
