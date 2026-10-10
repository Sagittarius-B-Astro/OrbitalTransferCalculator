#!/usr/bin/env python3
"""Reproduces every numerical claim in the README's research section.  Usage: python3 scripts/validate_claims.py
Each block prints the measured numbers and PASS/FAIL against the stated tolerance. Takes ~1 min."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from planechange.frames import Orbit, state_at
from planechange.pparam import velocities_from_p, time_of_flight, p_bounds, transfer_angle
from planechange.lambert_izzo import lambert
from planechange.kepler import propagate
from planechange.min_dv import min_delta_v, delta_v, _p_from_s
from planechange.nodal import min_delta_v_nodal

fails = []
def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok: fails.append(name)

# ---- A. Blanco (2025) worked example: r1=9000, r2=15000, dtheta=120 deg, GM=3.986e5
MUB = 3.986e5; R1, R2 = 9000.0, 15000.0
r1 = np.array([R1, 0, 0]); r2 = R2 * np.array([np.cos(np.radians(120)), np.sin(np.radians(120)), 0])
lo, hi, par = p_bounds(R1, R2, np.radians(120))
vc = np.sqrt(MUB / R1); p = 1.3128 * R1
v1, v2, th = velocities_from_p(r1, r2, p, MUB)
dv = np.linalg.norm(v1 - [0, vc, 0]) / vc
tof = time_of_flight(R1, R2, th, p, v1, r1, MUB) / (2 * np.pi * R1 / vc)
check("A1 p_par,1 = 0.6317 r1 and p_par,2 = 1.8173 r1", abs(lo/R1-0.6317) < 5e-4 and abs(par/R1-1.8173) < 5e-4, f"{lo/R1:.4f} {par/R1:.4f}")
check("A2 dv at p=1.3128 r1 is 0.1563 vc", abs(dv-0.1563) < 1e-4, f"{dv:.5f}")
check("A3 transfer time is 0.4607 Tc", abs(tof-0.4607) < 1e-4, f"{tof:.5f}")

# ---- B. p-family vs Izzo Lambert (max velocity difference over many geometries)
worst = 0.0; count = 0
for deg in (30, 60, 90, 120, 150, 170, 190, 210, 240, 270, 300, 330):
    t = np.radians(deg); rr2 = R2 * np.array([np.cos(t), np.sin(t), 0])
    lo_, hi_, par_ = p_bounds(R1, R2, t)
    top = par_ * 3 if np.isinf(hi_) else hi_
    for p in np.linspace(lo_ + 0.02 * (top - lo_), top * 0.98, 15):
        va, vb, th = velocities_from_p(r1, rr2, p, MUB)
        tf = time_of_flight(R1, R2, th, p, va, r1, MUB)
        m = min(lambert(r1, rr2, tf, MUB), key=lambda s: np.linalg.norm(s["v1"] - va))
        worst = max(worst, np.linalg.norm(m["v1"]-va), np.linalg.norm(m["v2"]-vb)); count += 1
check("B p-family == Izzo Lambert", worst < 1e-6, f"max |dv| = {worst:.2e} km/s over {count} cases")

# ---- C. Izzo solutions vs independent propagator (incl. multi-rev, long way)
worst = 0.0; count = 0
for deg in (60, 120, 175, 240, 300):
    t = np.radians(deg); rr2 = R2 * np.array([np.cos(t), np.sin(t), 0])
    for tf in (1000., 3915., 8000., 30000., 60000.):
        for s in lambert(r1, rr2, tf, MUB):
            r, v = propagate(r1, s["v1"], tf, MUB); worst = max(worst, np.linalg.norm(r - rr2)); count += 1
check("C Izzo endpoint error (propagated)", worst < 1e-6, f"max = {worst:.2e} km over {count} solutions")

# ---- D. Delta-t monotonic in p, M = 0 (dense scan, both ways, several geometries)
ok = True
for ratio in (0.6, 1.0, 1.67, 3.0):
    for deg in (40, 90, 120, 170, 190, 240, 300, 340):
        t = np.radians(deg); a = np.array([R1, 0, 0]); b = ratio*R1*np.array([np.cos(t), np.sin(t), 0])
        lo_, hi_, par_ = p_bounds(R1, ratio*R1, t)
        ps = np.linspace(lo_ + 1e-3*(par_), par_ * (3 if np.isinf(hi_) else 0.999), 400) if t < np.pi else np.linspace(1e-3*par_, hi_*0.999, 400)
        ts = []
        for p in ps:
            va, _, th = velocities_from_p(a, b, p, MUB); ts.append(time_of_flight(R1, ratio*R1, th, p, va, a, MUB))
        sgn = -1 if t < np.pi else 1
        ok &= bool(np.all(sgn*np.diff(ts) > 0))
check("D dt(p) strictly decreasing (short way) / increasing (long way)", ok, "4 radius ratios x 8 angles x 400 p-values")

# ---- E. dt upper limit -> one-sided bound on p:  p >= p*  <=>  dt(p) <= Tmax   (short way)
ok = True
for Tmax in (1500., 3000., 6000.):
    star = lambert(r1, r2, Tmax, MUB, max_revs=0)[0]; pstar = np.linalg.norm(np.cross(r1, star["v1"]))**2 / MUB
    for p in np.linspace(lo*1.05, 5*par, 300):
        va, _, th = velocities_from_p(r1, r2, p, MUB)
        ok &= (p >= pstar) == (time_of_flight(R1, R2, th, p, va, r1, MUB) <= Tmax*(1+1e-9))
    print(f"      Tmax={Tmax:.0f}s -> p* = {pstar:.2f} km (from one Lambert solve: p* = |r1 x v1|^2/mu)")
check("E dt <= Tmax  <=>  p >= p*(Tmax)", ok)

# ---- F. coplanar / pure-plane-change / combined burn
MU = 3.986004418e5
H = lambda a, b: (lambda at, v: abs(v(a, at)-v(a, a)) + abs(v(b, b)-v(b, at)))((a+b)/2, lambda r, aa: np.sqrt(MU*(2/r-1/aa)))
res = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000), MU, n_grid=16)
check("F1 coplanar circular == Hohmann", abs(res["dv"]-H(7000, 14000)) < 1e-4, f"{res['dv']:.6f} vs {H(7000,14000):.6f}")
di = np.radians(20)
res = min_delta_v(Orbit(7000, 7000), Orbit(7000, 7000, di), MU, n_grid=16); ex = 2*np.sqrt(MU/7000)*np.sin(di/2)
check("F2 pure 20 deg plane change == 2 v sin(di/2)", abs(res["dv"]-ex) < 1e-6, f"{res['dv']:.6f} vs {ex:.6f}")
res = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000, di), MU, n_grid=16)
vp, va = (np.sqrt(MU*(2/r-1/10500)) for r in (7000, 14000)); v1c, v2c = np.sqrt(MU/7000), np.sqrt(MU/14000)
naive = abs(vp-v1c) + np.sqrt(va**2+v2c**2-2*va*v2c*np.cos(di))
check("F3 combined optimum < Hohmann + whole plane change at apogee", res["dv"] < naive,
      f"{res['dv']:.4f} vs {naive:.4f} km/s ({100*(1-res['dv']/naive):.1f}% less)")
nod = min_delta_v_nodal(Orbit(7000, 7000), Orbit(14000, 14000, di), MU)
check("F4 nodal-slice solver == general search", abs(nod["dv"]-res["dv"]) < 1e-5, f"{nod['dv']:.6f} vs {res['dv']:.6f}")

# ---- G. Is dv(p) unimodal for fixed (nu1, nu2)?  (empirical only)
rng = np.random.default_rng(1); bad = n = 0
for _ in range(300):
    o1 = Orbit(rng.uniform(7000, 9000), rng.uniform(6800, 7000), rng.uniform(0, 1), rng.uniform(0, 6), rng.uniform(0, 6))
    o2 = Orbit(rng.uniform(10000, 30000), rng.uniform(9000, 10000), rng.uniform(0, 1), rng.uniform(0, 6), rng.uniform(0, 6))
    n1, n2 = rng.uniform(0, 6.28, 2)
    ra, _ = state_at(o1, n1, MU); rb, _ = state_at(o2, n2, MU)
    d = transfer_angle(ra, rb, o1.normal())
    if abs(np.sin(d)) < 1e-3: continue
    l, h, _ = p_bounds(np.linalg.norm(ra), np.linalg.norm(rb), d)
    v = np.array([delta_v(o1, o2, n1, n2, _p_from_s(s, l, h), MU) for s in np.linspace(-10, 10, 400)])
    n += 1; bad += int(np.sum((v[1:-1] < v[:-2]) & (v[1:-1] < v[2:])) > 1)
check("G dv(p) has a single interior minimum (empirical)", bad == 0, f"{bad} of {n} random geometries had more than one")

print("\nALL PASS" if not fails else f"\nFAILED: {fails}"); sys.exit(1 if fails else 0)
