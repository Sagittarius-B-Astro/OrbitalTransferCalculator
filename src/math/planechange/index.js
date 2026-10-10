/** Small 3-vector helpers and Keplerian orbit geometry (km, s, rad). */
export const TWO_PI = 2 * Math.PI;
const DEG = Math.PI / 180;

export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
export const mod2pi = (x) => ((x % TWO_PI) + TWO_PI) % TWO_PI;

/** Orbit from apoapsis/periapsis radii and orientation angles in RADIANS. P, Q = perifocal axes in the inertial frame, n = normal. */
export function makeOrbit({ ra, rp, inc = 0, raan = 0, argp = 0 }) {
  const a = (ra + rp) / 2;
  const e = (ra - rp) / (ra + rp);
  const p = a * (1 - e * e);
  const cO = Math.cos(raan), sO = Math.sin(raan), ci = Math.cos(inc), si = Math.sin(inc), cw = Math.cos(argp), sw = Math.sin(argp);
  // first two columns of Rz(RAAN) Rx(i) Rz(w)
  const P = [cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si];
  const Q = [-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si];
  return { ra, rp, a, e, p, P, Q, n: cross(P, Q) };
}

export const orbitFromDegrees = (ra, rp, iDeg, raanDeg, wDeg) =>
  makeOrbit({ ra, rp, inc: iDeg * DEG, raan: raanDeg * DEG, argp: wDeg * DEG });

/** Inertial [position, velocity] on `o` at true anomaly `nu`. */
export function stateAt(o, nu, mu) {
  const c = Math.cos(nu), s = Math.sin(nu);
  const r = o.p / (1 + o.e * c);
  const k = Math.sqrt(mu / o.p);
  const rx = r * c, ry = r * s, vx = -k * s, vy = k * (o.e + c);
  const { P, Q } = o;
  return [
    [rx * P[0] + ry * Q[0], rx * P[1] + ry * Q[1], rx * P[2] + ry * Q[2]],
    [vx * P[0] + vy * Q[0], vx * P[1] + vy * Q[1], vx * P[2] + vy * Q[2]],
  ];
}

/** True anomaly at which `o` crosses the in-plane direction u. */
export const anomalyOf = (o, u) => Math.atan2(dot(u, o.Q), dot(u, o.P));
