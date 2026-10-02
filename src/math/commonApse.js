import { M3_TO_KM3, deg2rad } from './common.js';

/** Elliptic transfer between coplanar orbits sharing an apse line (true anomalies A1, A2 in degrees). */
export function computeCommonApse({ r1a, r1p, r2a, r2p, A1, A2 }, muM3) {
  const mu = muM3 * M3_TO_KM3;
  const TA1 = deg2rad(A1);
  const TA2 = deg2rad(A2);
  const a1 = (r1a + r1p) / 2;
  const a2 = (r2a + r2p) / 2;
  const e1 = (r1a - r1p) / (r1a + r1p);
  const e2 = (r2a - r2p) / (r2a + r2p);
  const p1 = a1 * (1 - e1 ** 2);
  const p2 = a2 * (1 - e2 ** 2);
  const rA = p1 / (1 + e1 * Math.cos(TA1));
  const rB = p2 / (1 + e2 * Math.cos(TA2));

  const et = (rB - rA) / (rA * Math.cos(TA1) - rB * Math.cos(TA2));
  const pt = (rA * rB * (Math.cos(TA1) - Math.cos(TA2))) / (rA * Math.cos(TA1) - rB * Math.cos(TA2));
  const at = pt / (1 - et ** 2);

  const h1 = Math.sqrt(mu * p1);
  const ht = Math.sqrt(mu * pt);
  const vp1 = h1 / rA;
  const vpt = ht / rA;
  const vr1 = (mu / h1) * e1 * Math.sin(TA1);
  const vrt = (mu / ht) * et * Math.sin(TA1);
  const v1 = Math.hypot(vp1, vr1);
  const vt = Math.hypot(vpt, vrt);
  const phi1 = Math.atan2(vr1, vp1);
  const phit = Math.atan2(vrt, vpt);

  const totalDeltaV = Math.sqrt(v1 ** 2 + vt ** 2 - 2 * v1 * vt * Math.cos(phit - phi1));
  const gamma = Math.atan2(vrt - vr1, vpt - vp1);
  const transferTime = Math.PI * Math.sqrt(at ** 3 / mu); // NOTE: half a period; exact only for apse-to-apse arcs (see README)
  return { totalDeltaV, gamma, transferTime };
}
