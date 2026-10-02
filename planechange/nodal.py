"""Nodal (transfer angle = pi) two-impulse transfers.

When burn 1 and burn 2 are at OPPOSITE ends of the line of nodes, r1 and r2 are
anti-parallel, sin(dtheta) = 0, and the (nu1, nu2, p) parametrisation of pparam.py is
singular.  The family of conics through the two points is then

    p = 2 r1 r2 / (r1 + r2)                      (a single value)
    e cos(nu1) = (r2 - r1)/(r1 + r2) =: c0       (fixed)
    e sin(nu1) = q                               (FREE: dimensionless radial velocity at r1)
    phi                                          (FREE: orientation of the transfer plane
                                                  about the line of nodes)

with   v1 = sqrt(mu/p) q rhat1 + (sqrt(mu p)/r1) that(phi)
       v2 = sqrt(mu/p) q rhat1 - (sqrt(mu p)/r2) that(phi)
where that(phi) is a unit vector perpendicular to rhat1 and phi in [0, 2pi) sweeps every
transfer plane containing the line of nodes (both directions of motion).

For fixed (s, phi) the total delta-v is a sum of two norms of functions affine in q, hence
CONVEX in q: the inner minimisation is guaranteed unimodal.  Feasible q: e < 1
(|q| < 2 sqrt(r1 r2)/(r1+r2)), or e >= 1 with q < 0 (the arc from nu1 to nu1 + pi must not
cross the hyperbola's asymptote).
"""
import numpy as np
from .frames import state_at
from .pparam import time_of_flight
from .optimizers import golden_min

BIG = 1e12


def nodal_direction(o1, o2):
    """Unit vector along the line of nodes, or None if the planes coincide."""
    c = np.cross(o1.normal(), o2.normal())
    n = np.linalg.norm(c)
    return None if n < 1e-12 else c / n


def anomaly_of_direction(orbit, u):
    """True anomaly at which `orbit` crosses the (in-plane) direction u."""
    Q = orbit.dcm
    return float(np.arctan2(u @ Q[:, 1], u @ Q[:, 0]))


def q_max(r1n, r2n):
    """q at the parabola (e = 1): ellipses have |q| < q_max."""
    return 2 * np.sqrt(r1n * r2n) / (r1n + r2n)


def nodal_velocities(r1vec, r2n, mu, phi, q):
    """Departure/arrival velocity of the nodal transfer (see module docstring)."""
    r1n = np.linalg.norm(r1vec)
    rh = r1vec / r1n
    helper = np.array([1.0, 0, 0]) if abs(rh[0]) < 0.9 else np.array([0, 1.0, 0])
    ea = np.cross(rh, helper)
    ea /= np.linalg.norm(ea)
    eb = np.cross(rh, ea)
    that = np.cos(phi) * ea + np.sin(phi) * eb
    p = 2 * r1n * r2n / (r1n + r2n)
    vr = np.sqrt(mu / p) * q
    v1 = vr * rh + np.sqrt(mu * p) / r1n * that
    v2 = vr * rh - np.sqrt(mu * p) / r2n * that
    return v1, v2, p, np.cross(rh, that)  # last = unit angular-momentum direction


def feasible(q, r1n, r2n):
    return q < q_max(r1n, r2n) if q >= 0 else True


def nodal_delta_v(o1, o2, s, phi, q, mu):
    """Total delta-v for burn 1 at s*dhat on orbit 1 and burn 2 at -s*dhat on orbit 2."""
    d = nodal_direction(o1, o2)
    nu1, nu2 = anomaly_of_direction(o1, s * d), anomaly_of_direction(o2, -s * d)
    r1, vc1 = state_at(o1, nu1, mu)
    r2, vc2 = state_at(o2, nu2, mu)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    if not feasible(q, r1n, r2n):
        return BIG
    v1, v2, _, _ = nodal_velocities(r1, r2n, mu, phi, q)
    return float(np.linalg.norm(v1 - vc1) + np.linalg.norm(v2 - vc2))


def _best_q(o1, o2, s, phi, mu, Q=10.0):
    d = nodal_direction(o1, o2)
    r1n = np.linalg.norm(state_at(o1, anomaly_of_direction(o1, s * d), mu)[0])
    r2n = np.linalg.norm(state_at(o2, anomaly_of_direction(o2, -s * d), mu)[0])
    hi = q_max(r1n, r2n) * (1 - 1e-9)
    f = lambda q: nodal_delta_v(o1, o2, s, phi, q, mu)
    q, dv = golden_min(f, -Q, hi, tol=1e-12)
    return q, dv


def min_delta_v_nodal(o1, o2, mu, n_phi=72):
    """Best two-impulse transfer with burns at opposite ends of the line of nodes."""
    d = nodal_direction(o1, o2)
    if d is None:
        raise ValueError("orbital planes coincide: no unique line of nodes (use min_dv.min_delta_v)")
    best = None
    for s in (+1, -1):
        phis = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
        vals = [_best_q(o1, o2, s, ph, mu)[1] for ph in phis]
        for i in np.argsort(vals)[:2]:  # F(phi) can be multimodal: polish the two best cells
            h = 2 * np.pi / n_phi
            phi, _ = golden_min(lambda ph: _best_q(o1, o2, s, ph, mu)[1], phis[i] - h, phis[i] + h, tol=1e-11)
            q, dv = _best_q(o1, o2, s, phi, mu)
            if best is None or dv < best["dv"]:
                best = {"dv": dv, "s": s, "phi": phi % (2 * np.pi), "q": q}
    s, phi, q = best["s"], best["phi"], best["q"]
    nu1, nu2 = anomaly_of_direction(o1, s * d), anomaly_of_direction(o2, -s * d)
    r1, _ = state_at(o1, nu1, mu)
    r2, _ = state_at(o2, nu2, mu)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    v1, v2, p, nhat = nodal_velocities(r1, r2n, mu, phi, q)
    c0 = (r2n - r1n) / (r1n + r2n)
    best.update(nu1=nu1 % (2 * np.pi), nu2=nu2 % (2 * np.pi), p=p, e=float(np.hypot(c0, q)),
                tof=time_of_flight(r1n, r2n, np.pi, p, v1, r1, mu), v1=v1, v2=v2,
                plane_angle=float(np.arccos(np.clip(nhat @ o1.normal(), -1, 1))))
    return best
