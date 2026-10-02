import { M3_TO_KM3 } from './common.js';

/** Two-impulse Hohmann transfer between coplanar circular orbits. `muM3` in m^3/s^2, radii in km. */
export function computeHohmann({ r1, r2 }, muM3) {
  const mu = muM3 * M3_TO_KM3;
  const a = (r1 + r2) / 2;
  const ht = Math.sqrt((2 * mu * (r1 * r2)) / (r1 + r2));

  const v1 = Math.sqrt(mu / r1);
  const va = ht / r1;
  const vb = ht / r2;
  const v2 = Math.sqrt(mu / r2);

  const deltaV1 = Math.abs(va - v1);
  const deltaV2 = Math.abs(v2 - vb);
  return {
    deltaV1,
    deltaV2,
    totalDeltaV: deltaV1 + deltaV2,
    transferTime: Math.PI * Math.sqrt(a ** 3 / mu),
  };
}
