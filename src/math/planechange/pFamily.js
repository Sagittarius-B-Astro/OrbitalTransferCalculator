/**
 * Transfers through two fixed points parametrised by the semi-latus rectum p (no Lambert iteration).
 * way = +1: "short way" (transfer angle theta in (0, pi));  way = -1: "long way" (2*pi - theta).
 * Using BOTH ways at every (nu1, nu2) removes the prograde/retrograde branch jump that a one-sense search has.
 */
import { TWO_PI, dot, cross, norm, stateAt } from './geometry.js';
import { brentMin } from './optimize.js';

export const BIG = 1e12;
const N_SWEEP = 14;
const S_RANGE = { short: [-10, 8], long: [-10, 10] };

/** Valid p interval for a transfer angle: { lo, hi, par } with par = the connecting parabola. */
export function pBounds(r1n, r2n, dtheta) {
  const k = 2 * r1n * r2n * Math.sin(dtheta / 2) ** 2, l = r1n + r2n, w = 2 * Math.sqrt(r1n * r2n) * Math.cos(dtheta / 2);
  const pa = k / (l + w), pb = k / (l - w);
  return dtheta < Math.PI ? { lo: pa, hi: Infinity, par: pb } : { lo: 0, hi: pa, par: pa };
}

/** Unconstrained s -> valid p. Short way: scaled by the gap to the parabola (s = 0 is the parabola), so it is unit-free. */
export function pFromS(s, b) {
  const t = Math.max(-40, Math.min(40, s));
  return Number.isFinite(b.hi) ? b.hi / (1 + Math.exp(-t)) : b.lo + (b.par - b.lo) * Math.exp(t);
}

/** Everything that depends only on (nu1, nu2): the two burn points and the orbit velocities there. */
export function makeCell(o1, o2, nu1, nu2, mu) {
  const [r1, vc1] = stateAt(o1, nu1, mu), [r2, vc2] = stateAt(o2, nu2, mu);
  const r1n = norm(r1), r2n = norm(r2), c = norm(cross(r1, r2)), d = dot(r1, r2);
  return { r1, r2, vc1, vc2, r1n, r2n, mu, theta: Math.atan2(c, d), cosT: d / (r1n * r2n), sinT: c / (r1n * r2n),
           degenerate: c < 1e-4 * r1n * r2n };   // within ~1e-4 rad of r1 || r2: the p-formula divides by sin(theta) and amplifies rounding
                                                   // error until the optimiser starts exploiting noise. nodal.js owns that region exactly.
}

export const wayAngle = (cell, way) => (way > 0 ? cell.theta : TWO_PI - cell.theta);

/** Departure/arrival velocities on the conic with semi-latus rectum p (Lagrange coefficients). */
export function velocitiesOfP(cell, way, p) {
  const { r1, r2, r1n, r2n, mu } = cell, cd = cell.cosT, sd = way * cell.sinT;
  const inv = Math.sqrt(mu * p) / (r1n * r2n * sd);
  const fc = 1 - (r2n / p) * (1 - cd), gd = 1 - (r1n / p) * (1 - cd);
  return [r2.map((x, k) => (x - fc * r1[k]) * inv), r2.map((x, k) => (gd * x - r1[k]) * inv)];
}

/** Total delta-v for the transfer with semi-latus rectum p (hot loop: scalar arithmetic, no allocation). */
export function dvOfP(cell, way, p) {
  const { r1, r2, vc1, vc2, r1n, r2n, mu } = cell, cd = cell.cosT, sd = way * cell.sinT;
  const inv = Math.sqrt(mu * p) / (r1n * r2n * sd);
  const fc = 1 - (r2n / p) * (1 - cd), gd = 1 - (r1n / p) * (1 - cd);
  const ax = (r2[0] - fc * r1[0]) * inv - vc1[0], ay = (r2[1] - fc * r1[1]) * inv - vc1[1], az = (r2[2] - fc * r1[2]) * inv - vc1[2];
  const bx = vc2[0] - (gd * r2[0] - r1[0]) * inv, by = vc2[1] - (gd * r2[1] - r1[1]) * inv, bz = vc2[2] - (gd * r2[2] - r1[2]) * inv;
  const v = Math.hypot(ax, ay, az) + Math.hypot(bx, by, bz);
  return Number.isFinite(v) ? v : BIG;
}

/** min over p for one way: coarse sweep in s to bracket, then Brent. */
export function bestPForWay(cell, way) {
  const b = pBounds(cell.r1n, cell.r2n, wayAngle(cell, way));
  if (!(b.par - b.lo > 0) && Number.isFinite(b.par) && way > 0) return { dv: BIG };
  const [s0, s1] = way > 0 ? S_RANGE.short : S_RANGE.long;
  const f = (s) => dvOfP(cell, way, pFromS(s, b));
  const ds = (s1 - s0) / (N_SWEEP - 1);
  let bi = 0, bv = Infinity;
  for (let i = 0; i < N_SWEEP; i++) { const v = f(s0 + i * ds); if (v < bv) { bv = v; bi = i; } }
  const r = brentMin(f, s0 + Math.max(bi - 1, 0) * ds, s0 + Math.min(bi + 1, N_SWEEP - 1) * ds, 1e-10);
  return { dv: r.fx, p: pFromS(r.x, b), way };
}

/** F(nu1, nu2) = min over p and over both ways. */
export function cellMin(cell) {
  if (cell.degenerate) return { dv: BIG };
  const a = bestPForWay(cell, +1), b = bestPForWay(cell, -1);
  return a.dv <= b.dv ? a : b;
}
