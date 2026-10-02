"""Izzo's Lambert solver (Izzo 2015, "Revisiting Lambert's problem").

Iterates on the Lancaster-Blanchard variable x with Householder iterations and
reconstructs terminal velocities Gooding-style (Gooding 1990). Handles single-
and multi-revolution solutions. Time of flight is a REQUIRED input here; if you
want to treat time as free, use `pparam.py` instead (see README).
"""
import numpy as np
from .optimizers import halley, householder3


def _hyp2f1_3_1_52(z, tol=1e-11):
    """Gauss 2F1(3, 1; 5/2; z) by direct series (used by Battin's TOF form near x = 1)."""
    s, c, j = 1.0, 1.0, 0
    while True:
        c *= (3 + j) * (1 + j) / (2.5 + j) * z / (j + 1)
        s += c
        j += 1
        if abs(c) < tol or j > 500:
            return s


def x2tof(x, lam, M):
    """Dimensionless time of flight T(x) for revolutions M (Eq. 18; Battin/Lagrange near x = 1)."""
    dist = abs(x - 1.0)
    if 0.01 < dist < 0.2:  # Lagrange form: better conditioned in this window
        a = 1.0 / (1.0 - x * x)
        if a > 0:
            alfa = 2 * np.arccos(x)
            beta = 2 * np.arcsin(np.sqrt(lam * lam / a))
            beta = -beta if lam < 0 else beta
            return a * np.sqrt(a) * ((alfa - np.sin(alfa)) - (beta - np.sin(beta)) + 2 * np.pi * M) / 2
        alfa = 2 * np.arccosh(x)
        beta = 2 * np.arcsinh(np.sqrt(-lam * lam / a))
        beta = -beta if lam < 0 else beta
        return -a * np.sqrt(-a) * ((beta - np.sinh(beta)) - (alfa - np.sinh(alfa))) / 2
    K = lam * lam
    E = x * x - 1.0
    rho = abs(E)
    z = np.sqrt(1 + K * E)
    if dist <= 0.01:  # Battin series, Eq. (20)
        eta = z - lam * x
        S1 = 0.5 * (1 - lam - x * eta)
        Q = 4.0 / 3.0 * _hyp2f1_3_1_52(S1)
        return (eta ** 3 * Q + 4 * lam * eta) / 2 + (M * np.pi / rho ** 1.5 if M else 0.0)  # M=0 at x=1 would be 0/0
    y = np.sqrt(rho)
    g = x * z - lam * E
    if E < 0:  # elliptic
        d = M * np.pi + np.arccos(np.clip(g, -1, 1))
    else:      # hyperbolic
        d = np.log(y * (z - lam * x) + g)
    return (x - lam * z - d / y) / E


def _derivs(x, T, lam):
    """First three derivatives of T(x), Eq. (22)."""
    umx2 = 1.0 - x * x
    y = np.sqrt(1 - lam * lam * umx2)
    d1 = (3 * T * x - 2 + 2 * lam ** 3 * x / y) / umx2
    d2 = (3 * T + 5 * x * d1 + 2 * (1 - lam ** 2) * lam ** 3 / y ** 3) / umx2
    d3 = (7 * x * d2 + 8 * d1 - 6 * (1 - lam ** 2) * lam ** 5 * x / y ** 5) / umx2
    return d1, d2, d3


def _nudge(x):
    """The derivative formulas (22) are singular exactly at x = 1 (parabola); step off it."""
    return x if abs(x - 1.0) > 1e-9 else 1.0 + 1e-9


def _solve_x(T_target, lam, M, x0, tol=1e-11, max_iter=30):
    f = lambda x: x2tof(x, lam, M) - T_target
    df = lambda x: _derivs(_nudge(x), x2tof(_nudge(x), lam, M), lam)[0]
    d2f = lambda x: _derivs(_nudge(x), x2tof(_nudge(x), lam, M), lam)[1]
    d3f = lambda x: _derivs(_nudge(x), x2tof(_nudge(x), lam, M), lam)[2]
    return householder3(f, df, d2f, d3f, x0, tol=tol, max_iter=max_iter)


def find_xy(lam, T, max_revs=None):
    """Return list of (x, M) solutions: M = 0 first, then (left, right) pairs for M >= 1."""
    assert abs(lam) < 1, "|lambda| must be < 1 (zero-length chord / degenerate geometry)"
    assert T > 0, "dimensionless time of flight must be positive"
    M_max = int(np.floor(T / np.pi))
    T00 = np.arccos(lam) + lam * np.sqrt(1 - lam * lam)
    if T < T00 + M_max * np.pi and M_max > 0:
        # Locate Tmin(M_max) via Halley on dT/dx = 0 starting at x = 0
        d1 = lambda x: _derivs(x, x2tof(x, lam, M_max), lam)[0]
        d2 = lambda x: _derivs(x, x2tof(x, lam, M_max), lam)[1]
        xmin, _ = halley(d1, d2, lambda x: _derivs(x, x2tof(x, lam, M_max), lam)[2], 0.0)
        if x2tof(xmin, lam, M_max) > T:
            M_max -= 1
    if max_revs is not None:
        M_max = min(M_max, max_revs)

    T0 = T00
    T1 = 2.0 / 3.0 * (1 - lam ** 3)
    if T >= T0:
        x0 = (T0 / T) ** (2.0 / 3.0) - 1
    elif T < T1:
        x0 = 5.0 / 2.0 * T1 * (T1 - T) / (T * (1 - lam ** 5)) + 1
    else:
        x0 = (T0 / T) ** np.log2(T1 / T0) - 1
    sols = [(_solve_x(T, lam, 0, x0)[0], 0)]

    for M in range(1, M_max + 1):
        tl = ((M * np.pi + np.pi) / (8 * T)) ** (2.0 / 3.0)
        tr = (8 * T / (M * np.pi)) ** (2.0 / 3.0)
        sols.append((_solve_x(T, lam, M, (tl - 1) / (tl + 1))[0], M))  # x < 0 branch
        sols.append((_solve_x(T, lam, M, (tr - 1) / (tr + 1))[0], M))  # x > 0 branch
    return sols


def lambert(r1, r2, tof, mu, max_revs=None, plane_normal=(0.0, 0.0, 1.0), retrograde=False):
    """Solve Lambert's problem. Returns list of dicts {v1, v2, M, x}.

    `plane_normal` + `retrograde` decide the direction of motion: with the default
    the transfer is prograde about +z (short way if (r1 x r2).z > 0).
    """
    r1, r2 = np.asarray(r1, float), np.asarray(r2, float)
    c_vec = r2 - r1
    c, r1n, r2n = np.linalg.norm(c_vec), np.linalg.norm(r1), np.linalg.norm(r2)
    s = (c + r1n + r2n) / 2
    ir1, ir2 = r1 / r1n, r2 / r2n
    ih = np.cross(ir1, ir2)
    ih /= np.linalg.norm(ih)
    lam = np.sqrt(1 - c / s)
    if ih @ np.asarray(plane_normal, float) < 0:  # transfer angle > 180 deg about plane_normal
        lam = -lam
        it1, it2 = np.cross(ir1, ih), np.cross(ir2, ih)
    else:
        it1, it2 = np.cross(ih, ir1), np.cross(ih, ir2)
    if retrograde:
        lam, it1, it2 = -lam, -it1, -it2

    T = np.sqrt(2 * mu / s ** 3) * tof
    gamma = np.sqrt(mu * s / 2)
    rho = (r1n - r2n) / c
    sigma = np.sqrt(1 - rho ** 2)

    out = []
    for x, M in find_xy(lam, T, max_revs):
        y = np.sqrt(1 - lam ** 2 * (1 - x ** 2))
        vr1 = gamma * ((lam * y - x) - rho * (lam * y + x)) / r1n
        vr2 = -gamma * ((lam * y - x) + rho * (lam * y + x)) / r2n
        vt1 = gamma * sigma * (y + lam * x) / r1n
        vt2 = gamma * sigma * (y + lam * x) / r2n
        out.append({"v1": vr1 * ir1 + vt1 * it1, "v2": vr2 * ir2 + vt2 * it2, "M": M, "x": x})
    return out
