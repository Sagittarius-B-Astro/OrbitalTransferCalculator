/** Dependency-free 1-D and n-D minimisers. */

/** Brent's minimiser (parabolic interpolation + golden section) on [ax, cx]. Returns { x, fx, evals }. */
export function brentMin(f, ax, cx, tol = 1e-10, maxIter = 100) {
  const CGOLD = 0.3819660112501051, ZEPS = 1e-12;
  let a = Math.min(ax, cx), b = Math.max(ax, cx);
  let x = a + CGOLD * (b - a), w = x, v = x;
  let fx = f(x), fw = fx, fv = fx, d = 0, e = 0, evals = 1;
  for (let it = 0; it < maxIter; it++) {
    const xm = 0.5 * (a + b), tol1 = tol * Math.abs(x) + ZEPS, tol2 = 2 * tol1;
    if (Math.abs(x - xm) <= tol2 - 0.5 * (b - a)) break;
    let golden = true;
    if (Math.abs(e) > tol1) {
      let r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) p = -p;
      q = Math.abs(q);
      const etemp = e;
      e = d;
      if (!(Math.abs(p) >= Math.abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))) {
        d = p / q;
        const u = x + d;
        if (u - a < tol2 || b - u < tol2) d = xm >= x ? tol1 : -tol1;
        golden = false;
      }
    }
    if (golden) {
      e = x >= xm ? a - x : b - x;
      d = CGOLD * e;
    }
    const u = Math.abs(d) >= tol1 ? x + d : x + (d >= 0 ? tol1 : -tol1);
    const fu = f(u);
    evals++;
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; w = x; x = u; fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w === x) { v = w; w = u; fv = fw; fw = fu; }
      else if (fu <= fv || v === x || v === w) { v = u; fv = fu; }
    }
  }
  return { x, fx, evals };
}

/** n-D Nelder-Mead from a simplex of n+1 points. Returns { x, fx, iters }. */
export function nelderMead(f, simplex, { xtol = 1e-9, ftol = 1e-12, maxIter = 600 } = {}) {
  const n = simplex[0].length;
  let pts = simplex.map((p) => p.slice());
  let fs = pts.map(f);
  const comb = (a, b, t) => a.map((ai, k) => ai + t * (b[k] - ai));   // a + t (b - a)
  for (let it = 0; it < maxIter; it++) {
    const idx = fs.map((_, i) => i).sort((i, j) => fs[i] - fs[j]);
    pts = idx.map((i) => pts[i]);
    fs = idx.map((i) => fs[i]);
    let spread = 0;
    for (let i = 1; i <= n; i++) for (let k = 0; k < n; k++) spread = Math.max(spread, Math.abs(pts[i][k] - pts[0][k]));
    if (Math.abs(fs[n] - fs[0]) <= ftol && spread <= xtol) return { x: pts[0], fx: fs[0], iters: it };
    const c = new Array(n).fill(0);
    for (let i = 0; i < n; i++) for (let k = 0; k < n; k++) c[k] += pts[i][k] / n;
    const xr = comb(c, pts[n], -1), fr = f(xr);          // reflect
    if (fr >= fs[0] && fr < fs[n - 1]) { pts[n] = xr; fs[n] = fr; continue; }
    if (fr < fs[0]) {                                     // expand
      const xe = comb(c, pts[n], -2), fe = f(xe);
      if (fe < fr) { pts[n] = xe; fs[n] = fe; } else { pts[n] = xr; fs[n] = fr; }
      continue;
    }
    const outside = fr < fs[n];                           // contract
    const xc = outside ? comb(c, pts[n], -0.5) : comb(c, pts[n], 0.5), fc = f(xc);
    if (outside ? fc <= fr : fc < fs[n]) { pts[n] = xc; fs[n] = fc; continue; }
    for (let i = 1; i <= n; i++) { pts[i] = comb(pts[0], pts[i], 0.5); fs[i] = f(pts[i]); }   // shrink
  }
  return { x: pts[0], fx: fs[0], iters: maxIter };
}
