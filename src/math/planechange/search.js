/**
 * REFERENCE globalSearch (stand-in written to match the call sites in index.js and the tests, because the
 * original search.js content was not in the files provided for review).
 * 1) nGrid x nGrid scan of F(nu1, nu2) on the torus;  2) every grid local minimum (8-neighbourhood) is a basin seed;
 * 3) the best few seeds are polished with 2-D Nelder-Mead;  4) polished points are de-duplicated and sorted.
 * Returns { best:{nu1,nu2,dv}, basins:[{nu1,nu2,dv}] sorted by dv, nGrid, gridMin }.
 */
import { TWO_PI } from './geometry.js';
import { makeCell, cellMin, BIG } from './pFamily.js';
import { nelderMead } from './optimize.js';

const wrap = (x) => ((x % TWO_PI) + TWO_PI) % TWO_PI;
const angDist = (a, b) => { const d = Math.abs(wrap(a) - wrap(b)); return Math.min(d, TWO_PI - d); };

export function globalSearch(o1, o2, mu, { nGrid = 72, nPolish = 8, dedupe = 0.05 } = {}) {
  const h = TWO_PI / nGrid, off = 0.0137;                       // offset keeps the grid off exactly collinear points
  const F = (a, b) => cellMin(makeCell(o1, o2, a, b, mu)).dv;
  const G = Array.from({ length: nGrid }, (_, i) => Array.from({ length: nGrid }, (_, j) => F(off + i * h, off + j * h)));
  let gridMin = Infinity;
  const seeds = [];
  for (let i = 0; i < nGrid; i++) for (let j = 0; j < nGrid; j++) {
    const v = G[i][j];
    if (v < gridMin) gridMin = v;
    let isMin = v < BIG;
    for (let di = -1; di <= 1 && isMin; di++) for (let dj = -1; dj <= 1; dj++) {
      if ((di || dj) && G[(i + di + nGrid) % nGrid][(j + dj + nGrid) % nGrid] < v) { isMin = false; break; }
    }
    if (isMin) seeds.push({ i, j, v });
  }
  seeds.sort((a, b) => a.v - b.v);
  const polished = [];
  for (const s of seeds.slice(0, nPolish)) {
    const a = off + s.i * h, b = off + s.j * h;
    const r = nelderMead((x) => F(x[0], x[1]), [[a, b], [a + h / 3, b], [a, b + h / 3]], { xtol: 1e-10, ftol: 1e-13, maxIter: 800 });
    polished.push({ nu1: wrap(r.x[0]), nu2: wrap(r.x[1]), dv: Math.min(r.fx, s.v) });
  }
  polished.sort((p, q) => p.dv - q.dv);
  const basins = [];
  for (const p of polished) if (!basins.some((b) => angDist(b.nu1, p.nu1) < dedupe && angDist(b.nu2, p.nu2) < dedupe)) basins.push(p);
  return { best: basins[0], basins, nGrid, gridMin };
}
