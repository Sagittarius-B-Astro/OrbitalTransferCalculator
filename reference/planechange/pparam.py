"""Transfer conics parametrised by the semi-latus rectum p (no Kepler/Lambert solve).

For fixed end positions r1, r2 the conics connecting them form a ONE-parameter
family. Using p as that parameter gives terminal velocities in closed form via the
Lagrange coefficients,

    f = 1 - (r2/p)(1 - cos dtheta)      g = r1 r2 sin(dtheta) / sqrt(mu p)
    gdot = 1 - (r1/p)(1 - cos dtheta)
    v1 = (r2_vec - f r1_vec) / g        v2 = (gdot r2_vec - r1_vec) / g

so delta-v(nu1, nu2, p) is an explicit function. Time of flight is only needed
as an OUTPUT (or if a time constraint is imposed) and is obtained from Kepler's
equation. See README "Constraining the Lambert problem" for the derivation and
literature (Blanco 2025; Bate-Mueller-White p-iteration; Acta Astronautica 1984).
"""
import numpy as np


def transfer_angle(r1, r2, plane_normal=(0.0, 0.0, 1.0), retrograde=False):
    """Transfer angle in (0, 2pi) measured about plane_normal (flipped if retrograde)."""
    r1, r2 = np.asarray(r1, float), np.asarray(r2, float)
    n = np.asarray(plane_normal, float)
    c = np.clip(r1 @ r2 / (np.linalg.norm(r1) * np.linalg.norm(r2)), -1, 1)
    th = np.arccos(c)
    sign = np.sign(np.cross(r1, r2) @ n)
    if retrograde:
        sign = -sign
    return th if sign >= 0 else 2 * np.pi - th


def p_bounds(r1n, r2n, dtheta):
    """Open interval of p for which a connecting conic of the requested sense exists.

    Short way (dtheta < pi): p in (p_a, inf). The orbit is a (very large) ellipse as
    p -> p_a+, the connecting parabola sits at p_b, and hyperbolas lie beyond p_b.
    Long way (dtheta > pi): p in (0, p_a) and every connecting conic is an ellipse; the
    connecting parabola is the upper limit p_a.  Returns (p_lo, p_hi, p_parabola).
    """
    k = 2 * r1n * r2n * np.sin(dtheta / 2) ** 2
    l = r1n + r2n
    w = 2 * np.sqrt(r1n * r2n) * np.cos(dtheta / 2)
    p_a, p_b = k / (l + w), k / (l - w)
    if dtheta < np.pi:
        return p_a, np.inf, p_b
    return 0.0, p_a, p_a


def velocities_from_p(r1, r2, p, mu, plane_normal=(0.0, 0.0, 1.0), retrograde=False):
    """Departure/arrival velocity on the connecting conic with semi-latus rectum p."""
    r1, r2 = np.asarray(r1, float), np.asarray(r2, float)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    dth = transfer_angle(r1, r2, plane_normal, retrograde)
    cd, sd = np.cos(dth), np.sin(dth)
    f = 1 - r2n / p * (1 - cd)
    g = r1n * r2n * sd / np.sqrt(mu * p)
    gdot = 1 - r1n / p * (1 - cd)
    return (r2 - f * r1) / g, (gdot * r2 - r1) / g, dth


def time_of_flight(r1n, r2n, dtheta, p, v1, r1, mu, revolutions=0):
    """Time of flight on the conic fixed by p (Kepler's equation; +revolutions periods if bound)."""
    h = np.sqrt(mu * p)
    vr1 = (np.asarray(v1) @ np.asarray(r1)) / r1n
    ecos, esin = p / r1n - 1.0, vr1 * h / mu
    e = np.hypot(ecos, esin)
    nu1 = np.arctan2(esin, ecos)
    nu2 = nu1 + dtheta
    if abs(e - 1) < 1e-9:  # parabola (Barker)
        D = lambda nu: np.tan(nu / 2)
        return 0.5 * np.sqrt(p ** 3 / mu) * ((D(nu2) + D(nu2) ** 3 / 3) - (D(nu1) + D(nu1) ** 3 / 3))
    if e < 1:
        a = p / (1 - e * e)
        n = np.sqrt(mu / a ** 3)
        Ef = lambda nu: np.arctan2(np.sqrt(1 - e * e) * np.sin(nu), e + np.cos(nu))
        Mf = lambda nu: Ef(nu) - e * np.sin(Ef(nu))
        dM = (Mf(nu2) - Mf(nu1)) % (2 * np.pi)
        return (dM + 2 * np.pi * revolutions) / n
    a = p / (1 - e * e)  # negative
    n = np.sqrt(mu / (-a) ** 3)
    Hf = lambda nu: 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(nu / 2))
    Nf = lambda nu: e * np.sinh(Hf(nu)) - Hf(nu)
    return (Nf(nu2) - Nf(nu1)) / n
