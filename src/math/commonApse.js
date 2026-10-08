import { M3_TO_KM3, deg2rad } from './common.js';

const TWO_PI = 2 * Math.PI;
const wrap = (x) => ((x % TWO_PI) + TWO_PI) % TWO_PI;
function meanAnomaly(th, e) {   // signed e is fine, see above
  const E = Math.atan2(Math.sqrt(1 - e * e) * Math.sin(th), e + Math.cos(th));
  return E - e * Math.sin(E);
}

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

  if (!Number.isFinite(et) || Math.abs(et) >= 1)
    return { error: 'The common-apse conic through A and B is not an ellipse (|e_t| >= 1). Choose different true anomalies.' };
  const at = pt / (1 - et ** 2); 

  const h1 = Math.sqrt(mu * p1);
  const ht = Math.sqrt(mu * pt);
  const h2 = Math.sqrt(mu * p2);
  const vp1 = h1 / rA;
  const vptA = ht / rA;
  const vptB = ht / rB;
  const vp2 = h2 / rB;
  const vr1 = (mu / h1) * e1 * Math.sin(TA1);
  const vrtA = (mu / ht) * et * Math.sin(TA1);
  const vrtB = (mu / ht) * et * Math.sin(TA2);
  const vr2 = (mu / h2) * e2 * Math.sin(TA2);

  const deltaV1 = Math.sqrt((vp1 - vptA) ** 2 + (vr1 - vrtA) ** 2);
  const deltaV2 = Math.sqrt((vp2 - vptB) ** 2 + (vr2 - vrtB) ** 2);
  const totalDeltaV = deltaV1 + deltaV2;
  const gamma1 = Math.atan2(vrtA - vr1, vptA - vp1);
  const gamma2 = Math.atan2(vr2 - vrtB, vp2 - vptB);
  const transferTime = Math.sqrt(at ** 3 / mu) * wrap(meanAnomaly(TA2, et) - meanAnomaly(TA1, et));
  return { deltaV1, deltaV2, totalDeltaV, gamma1, gamma2, transferTime };
}