#!/usr/bin/env python3
"""Generate tests/js/fixtures/planechange_reference.json from the PYTHON planechange package.
Independent reference for the JS port: (a) F(nu1, nu2) = min over p and both ways, sampled at random points;
(b) best two-impulse delta-v from a high-effort Python search (n_grid=48, 12 polish starts) and from the Python nodal solver.
Run from the repo root:  python3 scripts/gen_planechange_fixtures.py        (takes ~1 min)"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from planechange.frames import Orbit, state_at
from planechange.pparam import p_bounds, transfer_angle
from planechange.optimizers import golden_min, nelder_mead
from planechange.nodal import min_delta_v_nodal

MU = 3.986004418e5
BIG = 1e12
R = np.radians

def p_of_s(s, lo, hi, par):                       # scale-aware mapping (same idea as the JS)
    s = np.clip(s, -40, 40)
    return hi / (1 + np.exp(-s)) if np.isfinite(hi) else lo + (par - lo) * np.exp(s)

def F_both(o1, o2, nu1, nu2):
    """min over p and over both ways (retrograde flag False/True == short/long way)."""
    r1, vc1 = state_at(o1, nu1, MU); r2, vc2 = state_at(o2, nu2, MU)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    best = BIG
    for retro in (False, True):
        th = transfer_angle(r1, r2, o1.normal(), retro); lo, hi, par = p_bounds(r1n, r2n, th); cd, sd = np.cos(th), np.sin(th)
        def dv(p):
            g = r1n * r2n * sd / np.sqrt(MU * p); fc = 1 - r2n / p * (1 - cd); gd = 1 - r1n / p * (1 - cd)
            return float(np.linalg.norm((r2 - fc * r1) / g - vc1) + np.linalg.norm(vc2 - (gd * r2 - r1) / g))
        f = lambda s: dv(p_of_s(s, lo, hi, par)); ss = np.linspace(-10, 10, 41); v = [f(s) for s in ss]; i = int(np.argmin(v))
        best = min(best, golden_min(f, ss[max(i - 1, 0)], ss[min(i + 1, 40)], tol=1e-11)[1])
    return best

def search(o1, o2, n_grid=48, n_refine=12):
    grid = np.linspace(0, 2 * np.pi, n_grid, endpoint=False) + 0.0137
    cells = sorted(((F_both(o1, o2, a, b), a, b) for a in grid for b in grid), key=lambda c: c[0])
    best = BIG
    for _, a, b in cells[:n_refine]:
        x, fx, _ = nelder_mead(lambda x: F_both(o1, o2, x[0], x[1]), [[a, b], [a + .05, b], [a, b + .05]], xtol=1e-9, ftol=1e-12, max_iter=500)
        best = min(best, fx)
    return best

def orbit(p, k):   # k = '1' or '2'
    return Orbit(p[f"r{k}a"], p[f"r{k}p"], R(p[f"i{k}"]), R(p[f"RAAN{k}"]), R(p[f"w{k}"]))

rng = np.random.default_rng(2024)
samples = []
for _ in range(14):
    p = dict(r1a=rng.uniform(7500, 12000), r1p=0, i1=rng.uniform(0, 40), RAAN1=rng.uniform(0, 360), w1=rng.uniform(0, 360),
             r2a=rng.uniform(13000, 40000), r2p=0, i2=rng.uniform(0, 70), RAAN2=rng.uniform(0, 360), w2=rng.uniform(0, 360))
    p["r1p"] = rng.uniform(6900, p["r1a"]); p["r2p"] = rng.uniform(8000, p["r2a"])
    n1, n2 = rng.uniform(0, 360, 2)
    samples.append(dict(params=p, nu1Deg=n1, nu2Deg=n2, dv=F_both(orbit(p, "1"), orbit(p, "2"), R(n1), R(n2))))

cases = {
 "circular 7000->14000, di=20": dict(r1a=7000, r1p=7000, i1=0, RAAN1=0, w1=0, r2a=14000, r2p=14000, i2=20, RAAN2=0, w2=0),
 "LEO->GEO 28.5": dict(r1a=7000, r1p=7000, i1=0, RAAN1=0, w1=0, r2a=42164, r2p=42164, i2=28.5, RAAN2=0, w2=0),
 "E1 elliptic->elliptic": dict(r1a=9000, r1p=7000, i1=10, RAAN1=20, w1=30, r2a=15000, r2p=10000, i2=25, RAAN2=60, w2=100),
 "E2 Molniya->GEO": dict(r1a=46000, r1p=7000, i1=63.4, RAAN1=0, w1=270, r2a=42164, r2p=42164, i2=0, RAAN2=0, w2=0),
 "E3 elliptic->circular": dict(r1a=12000, r1p=7000, i1=28.5, RAAN1=0, w1=0, r2a=7000, r2p=7000, i2=51.6, RAAN2=40, w2=0),
 "E4 random pair": dict(r1a=10220, r1p=9020, i1=15, RAAN1=86, w1=145, r2a=24500, r2p=10490, i2=7, RAAN2=348, w2=77),
}
optima = []
for name, p in cases.items():
    o1, o2 = orbit(p, "1"), orbit(p, "2")
    gen = search(o1, o2); nod = min_delta_v_nodal(o1, o2, MU)["dv"]
    optima.append(dict(name=name, params=p, pythonGeneral=gen, pythonNodal=nod, pythonBest=min(gen, nod))); print(f"{name:32s} general {gen:.6f} nodal {nod:.6f}")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "js", "fixtures", "planechange_reference.json")
os.makedirs(os.path.dirname(out), exist_ok=True)
json.dump(dict(mu_km3=MU, fSamples=samples, optima=optima), open(out, "w"), indent=1); print("wrote", os.path.normpath(out))
