"""Small, dependency-free (numpy only) root finders and minimizers.

Everything here is Pyodide friendly (no scipy). Each routine returns enough
information (iterations / value) for the tests to check convergence.
"""
import numpy as np


def halley(f, df, d2f, x0, tol=1e-12, max_iter=20):
    """Halley's method: x <- x - 2 f f' / (2 f'^2 - f f'')."""
    x = x0
    for it in range(max_iter):
        fx, d1, d2 = f(x), df(x), d2f(x)
        dx = 2 * fx * d1 / (2 * d1 ** 2 - fx * d2)
        x -= dx
        if abs(dx) < tol:
            return x, it + 1
    return x, max_iter


def householder3(f, df, d2f, d3f, x0, tol=1e-12, max_iter=20):
    """Third-order Householder iteration in the form used by Izzo (2015, Sec. 4.1)."""
    x = x0
    for it in range(max_iter):
        fx, d1, d2, d3 = f(x), df(x), d2f(x), d3f(x)
        num = d1 ** 2 - fx * d2 / 2
        den = d1 * (d1 ** 2 - fx * d2) + d3 * fx ** 2 / 6
        dx = fx * num / den
        x -= dx
        if abs(dx) < tol:
            return x, it + 1
    return x, max_iter


def brent_root(f, a, b, tol=1e-12, max_iter=200):
    """Brent's bracketing ROOT finder (not a minimizer). Requires f(a)*f(b) < 0."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("root not bracketed: f(a) and f(b) have the same sign")
    if abs(fa) < abs(fb):
        a, b, fa, fb = b, a, fb, fa
    c, fc = a, fa
    d = c
    mflag = True
    for it in range(max_iter):
        if fb == 0 or abs(b - a) < tol:
            return b, it
        if fa != fc and fb != fc:  # inverse quadratic interpolation
            s = (a * fb * fc / ((fa - fb) * (fa - fc))
                 + b * fa * fc / ((fb - fa) * (fb - fc))
                 + c * fa * fb / ((fc - fa) * (fc - fb)))
        else:  # secant
            s = b - fb * (b - a) / (fb - fa)
        lo, hi = sorted(((3 * a + b) / 4, b))
        if (not (lo < s < hi)
                or (mflag and abs(s - b) >= abs(b - c) / 2)
                or (not mflag and abs(s - b) >= abs(c - d) / 2)
                or (mflag and abs(b - c) < tol)
                or (not mflag and abs(c - d) < tol)):
            s = (a + b) / 2
            mflag = True
        else:
            mflag = False
        fs = f(s)
        d, c, fc = c, b, fb
        if fa * fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs
        if abs(fa) < abs(fb):
            a, b, fa, fb = b, a, fb, fa
    return b, max_iter


def golden_min(f, a, b, tol=1e-9, max_iter=200):
    """Golden-section MINIMIZER on [a, b]. Robust for unimodal f; returns (x, f(x))."""
    invphi = (np.sqrt(5) - 1) / 2
    c, d = b - invphi * (b - a), a + invphi * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(max_iter):
        if abs(b - a) < tol * (abs(c) + abs(d) + 1e-300):
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
    x = (a + b) / 2
    return x, f(x)


def nelder_mead(func, simplex, xtol=1e-9, ftol=1e-12, max_iter=2000,
                alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5):
    """n-dimensional Nelder-Mead. `simplex` is a list of n+1 points. Returns (x, f, iters)."""
    pts = [np.asarray(p, float) for p in simplex]
    fs = [func(p) for p in pts]
    for it in range(max_iter):
        order = np.argsort(fs)
        pts = [pts[i] for i in order]
        fs = [fs[i] for i in order]
        spread = max(np.linalg.norm(p - pts[0]) for p in pts[1:])
        if abs(fs[-1] - fs[0]) <= ftol and spread <= xtol:
            return pts[0], fs[0], it
        centroid = np.mean(pts[:-1], axis=0)
        xr = centroid + alpha * (centroid - pts[-1])
        fr = func(xr)
        if fs[0] <= fr < fs[-2]:
            pts[-1], fs[-1] = xr, fr
        elif fr < fs[0]:
            xe = centroid + gamma * (xr - centroid)
            fe = func(xe)
            pts[-1], fs[-1] = (xe, fe) if fe < fr else (xr, fr)
        else:
            if fr < fs[-1]:  # outside contraction
                xc = centroid + rho * (xr - centroid)
                fc = func(xc)
                ok = fc <= fr
            else:            # inside contraction
                xc = centroid + rho * (pts[-1] - centroid)
                fc = func(xc)
                ok = fc < fs[-1]
            if ok:
                pts[-1], fs[-1] = xc, fc
            else:            # shrink toward best
                for i in range(1, len(pts)):
                    pts[i] = pts[0] + sigma * (pts[i] - pts[0])
                    fs[i] = func(pts[i])
    best = int(np.argmin(fs))
    return pts[best], fs[best], max_iter
