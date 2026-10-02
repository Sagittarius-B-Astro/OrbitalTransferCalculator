/** Gravitational constant, m^3 kg^-1 s^-2. */
export const G = 6.674e-11;

/** Convert mu from m^3/s^2 to km^3/s^2 (all transfer math is done in km, s). */
export const M3_TO_KM3 = 1e-9;

export function computeCommonOrbitProperties(M, radius) {
  const mu = G * M; // m^3/s^2
  return { mu, radius };
}

export const deg2rad = (d) => (Math.PI * d) / 180;
