import { M3_TO_KM3 } from './common.js';

/** Three-impulse bi-elliptic transfer via intermediate apoapsis radius `ri` (km). */
export function computeBielliptic({ r1, r2, ri }, muM3) {
  const mu = muM3 * M3_TO_KM3;
  const a1i = (r1 + ri) / 2;
  const ai2 = (ri + r2) / 2;
  const ht1i = Math.sqrt((2 * mu * (r1 * ri)) / (r1 + ri));
  const hti2 = Math.sqrt((2 * mu * (ri * r2)) / (ri + r2));

  const v1 = Math.sqrt(mu / r1);
  const va1 = ht1i / r1;
  const vb1 = ht1i / ri;
  const vb2 = hti2 / ri;
  const vc2 = hti2 / r2;
  const v2 = Math.sqrt(mu / r2);

  const deltaV1 = Math.abs(va1 - v1);
  const deltaV2 = Math.abs(vb2 - vb1);
  const deltaV3 = Math.abs(v2 - vc2);
  return {
    deltaV1,
    deltaV2,
    deltaV3,
    totalDeltaV: deltaV1 + deltaV2 + deltaV3,
    transferTime: Math.PI * (Math.sqrt(a1i ** 3 / mu) + Math.sqrt(ai2 ** 3 / mu)),
  };
}
