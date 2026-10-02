"""Two-body propagation with universal variables (Curtis, Alg. 3.4).

Used to *verify* Lambert solutions independently of the Lambert code, and later
to sample transfer arcs for plotting.
"""
import numpy as np


def stumpff(z):
    if z > 1e-8:
        s = np.sqrt(z)
        return (1 - np.cos(s)) / z, (s - np.sin(s)) / s ** 3
    if z < -1e-8:
        s = np.sqrt(-z)
        return (np.cosh(s) - 1) / (-z), (np.sinh(s) - s) / s ** 3
    return 0.5 - z / 24, 1 / 6 - z / 120


def propagate(r0, v0, dt, mu, tol=1e-11, max_iter=200):
    """Return (r, v) after dt (any sign handled by caller; dt may be large/multi-rev)."""
    r0 = np.asarray(r0, float)
    v0 = np.asarray(v0, float)
    r0n, v0n = np.linalg.norm(r0), np.linalg.norm(v0)
    vr0 = r0 @ v0 / r0n
    alpha = 2 / r0n - v0n ** 2 / mu
    sm = np.sqrt(mu)
    chi = sm * abs(alpha) * dt if alpha != 0 else sm * dt / r0n
    for _ in range(max_iter):
        z = alpha * chi ** 2
        C, S = stumpff(z)
        F = r0n * vr0 / sm * chi ** 2 * C + (1 - alpha * r0n) * chi ** 3 * S + r0n * chi - sm * dt
        dF = r0n * vr0 / sm * chi * (1 - z * S) + (1 - alpha * r0n) * chi ** 2 * C + r0n
        step = F / dF
        chi -= step
        if abs(step) < tol:
            break
    z = alpha * chi ** 2
    C, S = stumpff(z)
    f = 1 - chi ** 2 / r0n * C
    g = dt - chi ** 3 / sm * S
    r = f * r0 + g * v0
    rn = np.linalg.norm(r)
    fdot = sm / (rn * r0n) * (z * S - 1) * chi
    gdot = 1 - chi ** 2 / rn * C
    return r, fdot * r0 + gdot * v0
