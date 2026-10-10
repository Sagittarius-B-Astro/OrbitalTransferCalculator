/**
 * Exact solver for the one place the (nu1, nu2, p) search is singular: burns at OPPOSITE ends of the line of
 * nodes (transfer angle = pi). There p = 2 r1 r2/(r1+r2) is fixed and the free parameters are
 *   q   = e sin(nu1)  (dimensionless radial velocity at r1)      -> delta-v is CONVEX in q
 *   phi = orientation of the transfer plane about the line of nodes (covers both senses of motion).
 *   v1 = sqrt(mu/p) q rhat1 + (sqrt(mu p)/r1) that(phi),   v2 = sqrt(mu/p) q rhat1 - (sqrt(mu p)/r2) that(phi)
 */
import { TWO_PI, cross, dot, norm, scale, stateAt, anomalyOf } from './geometry.js';
import { brentMin } from './optimize.js';

export function nodalDirection(o1, o2) {
  const c = cross(o1.n, o2.n), n = norm(c);
  return n < 1e-12 ? null : scale(c, 1 / n);
}

const qMax = (r1n, r2n) => (2 * Math.sqrt(r1n * r2n)) / (r1n + r2n);

function velocities(r1, r2n, mu, phi, q) {
  const r1n = norm(r1), rh = scale(r1, 1 / r1n);
  const helper = Math.abs(rh[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
  let ea = cross(rh, helper); ea = scale(ea, 1 / norm(ea));
  const eb = cross(rh, ea);
  const that = [0, 1, 2].map((k) => Math.cos(phi) * ea[k] + Math.sin(phi) * eb[k]);
  const p = (2 * r1n * r2n) / (r1n + r2n), vr = Math.sqrt(mu / p) * q;
  const t1 = Math.sqrt(mu * p) / r1n, t2 = Math.sqrt(mu * p) / r2n;
  return { p, v1: [0, 1, 2].map((k) => vr * rh[k] + t1 * that[k]), v2: [0, 1, 2].map((k) => vr * rh[k] - t2 * that[k]) };
}

export function minDeltaVNodal(o1, o2, mu, { nPhi = 72 } = {}) {
  const d = nodalDirection(o1, o2);
  if (!d) return null;
  let best = null;
  for (const s of [1, -1]) {
    const nu1 = anomalyOf(o1, scale(d, s)), nu2 = anomalyOf(o2, scale(d, -s));
    const [r1, vc1] = stateAt(o1, nu1, mu), [r2, vc2] = stateAt(o2, nu2, mu);
    const r1n = norm(r1), r2n = norm(r2), qHi = qMax(r1n, r2n) * (1 - 1e-9);
    const dv = (phi, q) => { const { v1, v2 } = velocities(r1, r2n, mu, phi, q);
      return Math.hypot(v1[0] - vc1[0], v1[1] - vc1[1], v1[2] - vc1[2]) + Math.hypot(v2[0] - vc2[0], v2[1] - vc2[1], v2[2] - vc2[2]); };
    const bestQ = (phi) => brentMin((q) => dv(phi, q), -10, qHi, 1e-12);        // convex in q
    const h = TWO_PI / nPhi, vals = [];
    for (let i = 0; i < nPhi; i++) vals.push(bestQ(i * h).fx);
    const order = vals.map((v, i) => i).sort((i, j) => vals[i] - vals[j]).slice(0, 2);   // F(phi) can be multimodal
    for (const i of order) {
      const rp = brentMin((phi) => bestQ(phi).fx, i * h - h, i * h + h, 1e-11), q = bestQ(rp.x).x;
      if (!best || rp.fx < best.dv) best = { dv: rp.fx, s, phi: rp.x, q, nu1, nu2, r1, r2, vc1, vc2 };
    }
  }
  const { v1, v2, p } = velocities(best.r1, norm(best.r2), mu, best.phi, best.q);
  return { ...best, v1, v2, p, dtheta: Math.PI };
}
