import { M3_TO_KM3 } from '../common.js';
import { TWO_PI, norm, cross, sub, orbitFromDegrees } from './geometry.js';
import { makeCell, cellMin, velocitiesOfP, wayAngle } from './pFamily.js';
import { globalSearch } from './search.js';
import { minDeltaVNodal } from './nodal.js';
import { timeOfFlight, conicForPlot, minRadiusOnArc } from './kepler.js';

const deg = (x) => (x * 180) / Math.PI;

function describe(method, c, mu) {
  const dv1 = norm(sub(c.v1, c.vc1)), dv2 = norm(sub(c.vc2, c.v2));
  const conic = conicForPlot(c.r1, c.v1, c.r2, mu);
  const e = conic.e;
  return {
    method, totalDeltaV: dv1 + dv2, deltaV1: dv1, deltaV2: dv2,
    transferTime: timeOfFlight(c.p, c.dtheta, c.r1, c.v1, mu),
    nu1Deg: deg(c.nu1) % 360, nu2Deg: deg(c.nu2) % 360,
    p: c.p, e, a: Math.abs(e - 1) > 1e-9 ? c.p / (1 - e * e) : null,
    transferPeriapsis: c.p / (1 + e), minRadius: minRadiusOnArc(conic), conic,
  };
}

/**
 * Best two-impulse transfer between two orbits (free time of flight).
 * params: r1a, r1p, i1, RAAN1, w1, r2a, r2p, i2, RAAN2, w2  (km, degrees); muM3 in m^3/s^2 like the other solvers.
 */
export function computePlaneChange(params, muM3, { nGrid = 72 } = {}) {
  for (const k of ['1', '2']) {
    const ra = params[`r${k}a`], rp = params[`r${k}p`];
    if (!(rp > 0) || !(ra >= rp)) return { error: `Orbit ${k}: apoapsis radius must be >= periapsis radius, and both positive.` };
  }
  const mu = muM3 * M3_TO_KM3;
  const o1 = orbitFromDegrees(params.r1a, params.r1p, params.i1, params.RAAN1, params.w1);
  const o2 = orbitFromDegrees(params.r2a, params.r2p, params.i2, params.RAAN2, params.w2);

  const g = globalSearch(o1, o2, mu, { nGrid });
  const cell = makeCell(o1, o2, g.best.nu1, g.best.nu2, mu), pick = cellMin(cell);
  const [v1, v2] = velocitiesOfP(cell, pick.way, pick.p);
  let result = describe('general', { v1, v2, vc1: cell.vc1, vc2: cell.vc2, r1: cell.r1, r2: cell.r2, p: pick.p,
                                     dtheta: wayAngle(cell, pick.way), nu1: g.best.nu1, nu2: g.best.nu2 }, mu);

  if (norm(cross(o1.n, o2.n)) > 1e-6) {                       // planes differ -> exact nodal-slice candidate
    const n = minDeltaVNodal(o1, o2, mu);
    if (n && n.dv < result.totalDeltaV) result = describe('nodal', n, mu);
  }
  result.basins = g.basins.map((b) => ({ nu1Deg: deg(b.nu1), nu2Deg: deg(b.nu2), deltaV: b.dv }));
  return result;
}

export { globalSearch, minDeltaVNodal };
