import { M3_TO_KM3, deg2rad } from './common.js';

/**
 * Single-impulse transfer between coplanar ellipses whose apse lines are rotated by A degrees.
 * Intersection of the two conics: p1/(1+e1 cos(th)) = p2/(1+e2 cos(th - eta)).
 *
 * Fixes vs. the original: (1) cos(A) used degrees instead of eta, (2) the radius used eta
 * where e1 belongs, (3) acos(...) % 2pi had the wrong operator precedence.
 */
export function computeApseRotate({ r1a, r1p, r2a, r2p, A }, muM3) {
  const mu = muM3 * M3_TO_KM3;
  const eta = deg2rad(A);
  const a1 = (r1a + r1p) / 2;
  const a2 = (r2a + r2p) / 2;
  const e1 = (r1a - r1p) / (r1a + r1p);
  const e2 = (r2a - r2p) / (r2a + r2p);
  const p1 = a1 * (1 - e1 ** 2);
  const p2 = a2 * (1 - e2 ** 2);

  const a = e1 * p2 - e2 * p1 * Math.cos(eta);
  const b = -e2 * p1 * Math.sin(eta);
  const c = p1 - p2;

  const alpha = Math.atan2(b, a);
  const arg = (c / a) * Math.cos(alpha);
  if (Math.abs(arg) > 1) return { error: 'Orbits do not intersect for this apse rotation.' };

  const TA1 = (alpha - Math.acos(arg)) % (2 * Math.PI); // true anomaly on orbit 1 (second root: alpha + acos)
  const TA2 = (TA1 - eta) % (2 * Math.PI);              // true anomaly on orbit 2 at the same point

  const h1 = Math.sqrt(mu * p1);
  const h2 = Math.sqrt(mu * p2);
  const r = p1 / (1 + e1 * Math.cos(TA1));
  const vp1 = h1 / r;
  const vp2 = h2 / r;
  const vr1 = (mu / h1) * e1 * Math.sin(TA1);
  const vr2 = (mu / h2) * e2 * Math.sin(TA2);
  const v1 = Math.hypot(vp1, vr1);
  const v2 = Math.hypot(vp2, vr2);
  const phi1 = Math.atan2(vr1, vp1);
  const phi2 = Math.atan2(vr2, vp2);

  const totalDeltaV = Math.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * Math.cos(phi2 - phi1));
  const gamma = Math.atan2(vr2 - vr1, vp2 - vp1);
  return { totalDeltaV, gamma, TA1, TA2, r };
}
