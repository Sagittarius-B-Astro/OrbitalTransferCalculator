import { TWO_PI, dot, cross, norm } from './geometry.js';

/** Time of flight on the conic of semi-latus rectum p from r1 (velocity v1) through transfer angle dtheta. */
export function timeOfFlight(p, dtheta, r1, v1, mu) {
  const r1n = norm(r1), h = Math.sqrt(mu * p), vr1 = dot(r1, v1) / r1n;
  const ecos = p / r1n - 1, esin = (vr1 * h) / mu, e = Math.hypot(ecos, esin);
  const nu1 = Math.atan2(esin, ecos), nu2 = nu1 + dtheta;
  if (Math.abs(e - 1) < 1e-9) {                                   // parabola (Barker)
    const D = (nu) => Math.tan(nu / 2);
    const B = (nu) => D(nu) + D(nu) ** 3 / 3;
    return 0.5 * Math.sqrt(p ** 3 / mu) * (B(nu2) - B(nu1));
  }
  if (e < 1) {
    const a = p / (1 - e * e), n = Math.sqrt(mu / a ** 3);
    const M = (nu) => { const E = Math.atan2(Math.sqrt(1 - e * e) * Math.sin(nu), e + Math.cos(nu)); return E - e * Math.sin(E); };
    return (((M(nu2) - M(nu1)) % TWO_PI) + TWO_PI) % TWO_PI / n;
  }
  const a = p / (1 - e * e), n = Math.sqrt(mu / (-a) ** 3);
  const H = (nu) => 2 * Math.atanh(Math.sqrt((e - 1) / (e + 1)) * Math.tan(nu / 2));
  const N = (nu) => e * Math.sinh(H(nu)) - H(nu);
  return (N(nu2) - N(nu1)) / n;
}

/** Closed-form description of the arc: r(nu) = p/(1+e cos nu) (cos nu P + sin nu Q), nu in [nuStart, nuEnd]. */
export function conicForPlot(r1, v1, r2, mu) {
  const h = cross(r1, v1), hn = norm(h), hh = h.map((x) => x / hn);
  const vxh = cross(v1, h), r1n = norm(r1);
  const evec = [vxh[0] / mu - r1[0] / r1n, vxh[1] / mu - r1[1] / r1n, vxh[2] / mu - r1[2] / r1n];
  const e = norm(evec);
  const P = e > 1e-9 ? evec.map((x) => x / e) : r1.map((x) => x / r1n);   // circular: periapsis direction arbitrary
  const Q = cross(hh, P);
  const nuStart = Math.atan2(dot(r1, Q), dot(r1, P));
  const dth = ((Math.atan2(dot(cross(r1, r2), hh), dot(r1, r2)) % TWO_PI) + TWO_PI) % TWO_PI;
  return { p: (hn * hn) / mu, e, P, Q, nuStart, nuEnd: nuStart + dth };
}

/** Smallest distance from the planet's centre ON THE TRAVELLED ARC (not the whole conic). */
export function minRadiusOnArc(c) {
  const r = (nu) => c.p / (1 + c.e * Math.cos(nu));
  const k = Math.ceil(c.nuStart / TWO_PI);
  if (k * TWO_PI <= c.nuEnd) return c.p / (1 + c.e);
  return Math.min(r(c.nuStart), r(c.nuEnd));
}
