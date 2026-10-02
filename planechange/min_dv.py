"""Minimum delta-v two-impulse transfer between two orbits (free time of flight).

Decision variables: true anomaly nu1 on orbit 1, true anomaly nu2 on orbit 2 and
the semi-latus rectum p of the connecting conic (Acta Astronautica 11, 1984 uses
the same three variables). No Lambert iteration and no time-of-flight sampling is
required: extra revolutions only add whole periods and never change delta-v.
"""
import numpy as np
from .frames import state_at
from .pparam import velocities_from_p, p_bounds, time_of_flight
from .optimizers import golden_min, nelder_mead

BIG = 1e12  # penalty for infeasible / degenerate geometry


def _p_from_s(s, p_lo, p_hi):
    """Unconstrained s -> open interval (p_lo, p_hi) so Nelder-Mead never leaves it."""
    if np.isinf(p_hi):
        return p_lo + np.exp(np.clip(s, -40, 40))
    return p_lo + (p_hi - p_lo) / (1 + np.exp(-np.clip(s, -40, 40)))


def delta_v(orbit1, orbit2, nu1, nu2, p, mu, retrograde=False):
    """Total delta-v for the two-impulse transfer defined by (nu1, nu2, p)."""
    r1, vc1 = state_at(orbit1, nu1, mu)
    r2, vc2 = state_at(orbit2, nu2, mu)
    if np.linalg.norm(np.cross(r1, r2)) < 1e-9 * np.linalg.norm(r1) * np.linalg.norm(r2):
        return BIG  # collinear: transfer plane undefined (see README, "nodal case")
    v1, v2, _ = velocities_from_p(r1, r2, p, mu, plane_normal=orbit1.normal(), retrograde=retrograde)
    return float(np.linalg.norm(v1 - vc1) + np.linalg.norm(vc2 - v2))


def _best_p(orbit1, orbit2, nu1, nu2, mu, retrograde, n_p=24):
    """For fixed (nu1, nu2) minimise over p: coarse sweep in s, then golden-section polish."""
    r1, _ = state_at(orbit1, nu1, mu)
    r2, _ = state_at(orbit2, nu2, mu)
    from .pparam import transfer_angle
    th = transfer_angle(r1, r2, orbit1.normal(), retrograde)
    p_lo, p_hi, _ = p_bounds(np.linalg.norm(r1), np.linalg.norm(r2), th)
    f = lambda s: delta_v(orbit1, orbit2, nu1, nu2, _p_from_s(s, p_lo, p_hi), mu, retrograde)
    ss = np.linspace(-8, 8, n_p)
    vals = [f(s) for s in ss]
    i = int(np.argmin(vals))
    a, b = ss[max(i - 1, 0)], ss[min(i + 1, n_p - 1)]
    s_best, dv = golden_min(f, a, b, tol=1e-10)
    return _p_from_s(s_best, p_lo, p_hi), dv


def min_delta_v(orbit1, orbit2, mu, n_grid=36, n_refine=3, retrograde=False):
    """Coarse (nu1, nu2) grid with inner p-minimisation, then 3-D Nelder-Mead polish.

    Returns dict(dv, nu1, nu2, p, tof, a, e). Angles in radians.
    """
    grid = np.linspace(0, 2 * np.pi, n_grid, endpoint=False) + 1e-3  # avoid exact collinear points
    cells = []
    for i, n1 in enumerate(grid):
        for j, n2 in enumerate(grid):
            p, dv = _best_p(orbit1, orbit2, n1, n2, mu, retrograde)
            cells.append((dv, n1, n2, p))
    cells.sort(key=lambda c: c[0])

    best = None
    for dv0, n1, n2, p0 in cells[:n_refine]:
        r1, _ = state_at(orbit1, n1, mu)
        r2, _ = state_at(orbit2, n2, mu)
        from .pparam import transfer_angle
        th = transfer_angle(r1, r2, orbit1.normal(), retrograde)
        p_lo, p_hi, _ = p_bounds(np.linalg.norm(r1), np.linalg.norm(r2), th)
        # recover s0 from p0 (inverse of _p_from_s)
        if np.isinf(p_hi):
            s0 = np.log(max(p0 - p_lo, 1e-300))
        else:
            q = (p0 - p_lo) / (p_hi - p_lo)
            s0 = np.log(q / (1 - q))

        def obj(x):
            a1, a2, s = x
            r1_, _ = state_at(orbit1, a1, mu)
            r2_, _ = state_at(orbit2, a2, mu)
            th_ = transfer_angle(r1_, r2_, orbit1.normal(), retrograde)
            lo, hi, _ = p_bounds(np.linalg.norm(r1_), np.linalg.norm(r2_), th_)
            return delta_v(orbit1, orbit2, a1, a2, _p_from_s(s, lo, hi), mu, retrograde)

        h = 0.05
        simplex = [[n1, n2, s0], [n1 + h, n2, s0], [n1, n2 + h, s0], [n1, n2, s0 + 0.2]]
        x, fx, _ = nelder_mead(obj, simplex, xtol=1e-10, ftol=1e-13, max_iter=3000)
        if best is None or fx < best[0]:
            best = (fx, x)

    dv, (a1, a2, s) = best
    r1, vc1 = state_at(orbit1, a1, mu)
    r2, vc2 = state_at(orbit2, a2, mu)
    th = transfer_angle(r1, r2, orbit1.normal(), retrograde)
    lo, hi, _ = p_bounds(np.linalg.norm(r1), np.linalg.norm(r2), th)
    p = _p_from_s(s, lo, hi)
    v1, v2, _ = velocities_from_p(r1, r2, p, mu, orbit1.normal(), retrograde)
    tof = time_of_flight(np.linalg.norm(r1), np.linalg.norm(r2), th, p, v1, r1, mu)
    eps = np.linalg.norm(v1) ** 2 / 2 - mu / np.linalg.norm(r1)
    a = -mu / (2 * eps)
    e = np.sqrt(max(0.0, 1 - p / a)) if a > 0 else np.sqrt(1 - p / a)
    return {"dv": dv, "nu1": a1 % (2 * np.pi), "nu2": a2 % (2 * np.pi), "p": p, "tof": tof, "a": a, "e": e,
            "v1": v1, "v2": v2}
