import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { computePlaneChange } from '../../src/math/planechange/index.js';
import { orbitFromDegrees, stateAt, dot, cross } from '../../src/math/planechange/geometry.js';
import { makeCell, cellMin, bestPForWay } from '../../src/math/planechange/pFamily.js';
import { brentMin, nelderMead } from '../../src/math/planechange/optimize.js';
import { timeOfFlight, minRadiusOnArc } from '../../src/math/planechange/kepler.js';
import { globalSearch } from '../../src/math/planechange/search.js';

const ref = JSON.parse(fs.readFileSync(new URL('./fixtures/planechange_reference.json', import.meta.url)));
const MU_KM = ref.mu_km3, MU_M3 = MU_KM * 1e9, RAD = Math.PI / 180;
const near = (a, b, tol, msg = '') => assert.ok(Math.abs(a - b) <= tol, `${msg} ${a} !~ ${b} (tol ${tol})`);
const orb = (p, k) => orbitFromDegrees(p[`r${k}a`], p[`r${k}p`], p[`i${k}`], p[`RAAN${k}`], p[`w${k}`]);

test('F(nu1, nu2) matches the Python implementation at 14 random points (min over p and both ways)', () => {
  for (const s of ref.fSamples) {
    const f = cellMin(makeCell(orb(s.params, 1), orb(s.params, 2), s.nu1Deg * RAD, s.nu2Deg * RAD, MU_KM)).dv;
    near(f, s.dv, 1e-7 * Math.max(1, s.dv), 'F');
  }
});

test('global optimum is never worse than the Python high-effort search (and prints any improvement)', () => {
  for (const c of ref.optima) {
    const r = computePlaneChange(c.params, MU_M3);
    const tol = c.nearNodal ? 2e-4 : 1e-6;   // near-nodal Python values are rounding-noise artefacts (see README), so allow slack
    assert.ok(r.totalDeltaV <= c.pythonBest + tol, `${c.name}: JS ${r.totalDeltaV} vs Python ${c.pythonBest}`);
    if (r.totalDeltaV < c.pythonBest - 1e-6) console.log(`  note: JS beat Python on "${c.name}" by ${(c.pythonBest - r.totalDeltaV).toExponential(2)} km/s`);
  }
});

test('pure plane change = 2 v sin(di/2), exact', () => {
  const di = 20 * RAD, v = Math.sqrt(MU_KM / 7000);
  const r = computePlaneChange({ r1a: 7000, r1p: 7000, i1: 0, RAAN1: 0, w1: 0, r2a: 7000, r2p: 7000, i2: 20, RAAN2: 0, w2: 0 }, MU_M3);
  near(r.totalDeltaV, 2 * v * Math.sin(di / 2), 1e-9);
});

test('circular -> circular equals the textbook Hohmann + optimally split plane change', () => {
  const [r1, r2, di] = [7000, 14000, 20 * RAD], at = (r1 + r2) / 2, vis = (r, a) => Math.sqrt(MU_KM * (2 / r - 1 / a));
  const [vp, va, v1, v2] = [vis(r1, at), vis(r2, at), vis(r1, r1), vis(r2, r2)];
  const f = (a1) => Math.sqrt(v1 ** 2 + vp ** 2 - 2 * v1 * vp * Math.cos(a1)) + Math.sqrt(va ** 2 + v2 ** 2 - 2 * va * v2 * Math.cos(di - a1));
  const best = brentMin(f, 0, di, 1e-13).fx;
  const r = computePlaneChange({ r1a: r1, r1p: r1, i1: 0, RAAN1: 0, w1: 0, r2a: r2, r2p: r2, i2: 20, RAAN2: 0, w2: 0 }, MU_M3);
  near(r.totalDeltaV, best, 1e-8);
  assert.equal(r.method, 'nodal');
  near(r.transferTime, Math.PI * Math.sqrt(at ** 3 / MU_KM), 1e-4);
});

test('coplanar circular orbits recover the Hohmann delta-v (no nodal line: general search only)', () => {
  const [r1, r2] = [7000, 14000], at = (r1 + r2) / 2, vis = (r, a) => Math.sqrt(MU_KM * (2 / r - 1 / a));
  const h = Math.abs(vis(r1, at) - vis(r1, r1)) + Math.abs(vis(r2, r2) - vis(r2, at));
  const r = computePlaneChange({ r1a: r1, r1p: r1, i1: 0, RAAN1: 0, w1: 0, r2a: r2, r2p: r2, i2: 0, RAAN2: 0, w2: 0 }, MU_M3);
  near(r.totalDeltaV, h, 1e-5);
  assert.equal(r.method, 'general');
});

test('result is self-consistent: burns sum to the total, conic passes through both burn points', () => {
  const p = ref.optima.find((o) => o.name.startsWith('E1')).params, r = computePlaneChange(p, MU_M3), c = r.conic;
  near(r.deltaV1 + r.deltaV2, r.totalDeltaV, 1e-12);
  const at = (nu) => { const k = c.p / (1 + c.e * Math.cos(nu)); return [0, 1, 2].map((i) => k * (Math.cos(nu) * c.P[i] + Math.sin(nu) * c.Q[i])); };
  const [r1] = stateAt(orb(p, 1), r.nu1Deg * RAD, MU_KM), [r2] = stateAt(orb(p, 2), r.nu2Deg * RAD, MU_KM);
  near(Math.hypot(...at(c.nuStart).map((x, i) => x - r1[i])), 0, 1e-6, 'start');
  near(Math.hypot(...at(c.nuEnd).map((x, i) => x - r2[i])), 0, 1e-6, 'end');
  assert.ok(c.nuEnd > c.nuStart && c.nuEnd - c.nuStart < 2 * Math.PI);
  assert.ok(r.minRadius > 0 && r.minRadius <= Math.min(Math.hypot(...r1), Math.hypot(...r2)) + 1e-6);
});

test('F is continuous across the old prograde/retrograde branch flip (the prograde-labelled value jumps by ~1.3 km/s there)', () => {
  const p = ref.optima.find((o) => o.name.startsWith('E1')).params, o1 = orb(p, 1), o2 = orb(p, 2);
  const cellAt = (n2) => makeCell(o1, o2, 100 * RAD, n2 * RAD, MU_KM);
  const flip = 352.327649, a = flip - 2e-4, b = flip + 2e-4;
  near(cellMin(cellAt(a)).dv, cellMin(cellAt(b)).dv, 1e-3, 'best of both ways');
  // the old code kept only the way that is prograde about orbit 1's normal: short if (r1 x r2).n1 >= 0, else long
  const prograde = (n2) => { const c = cellAt(n2), s = dot(cross(c.r1, c.r2), o1.n); return { s, dv: bestPForWay(c, s >= 0 ? +1 : -1).dv }; };
  const A = prograde(a), B = prograde(b);
  assert.ok(A.s * B.s < 0, 'sign of (r1 x r2).n1 must flip across this point');
  assert.ok(Math.abs(A.dv - B.dv) > 1.0, `prograde-only value should jump, got ${Math.abs(A.dv - B.dv)}`);
  // each geometric way on its own is continuous
  for (const w of [+1, -1]) near(bestPForWay(cellAt(a), w).dv, bestPForWay(cellAt(b), w).dv, 1e-3, `way ${w}`);
});

test('basins are enumerated, sorted, and the best one is returned', () => {
  const p = ref.optima.find((o) => o.name.startsWith('E1')).params, r = computePlaneChange(p, MU_M3);
  assert.ok(r.basins.length >= 2);
  for (let i = 1; i < r.basins.length; i++) assert.ok(r.basins[i].deltaV >= r.basins[i - 1].deltaV);
  near(r.basins[0].deltaV, r.totalDeltaV, 1e-9);
});

test('answer does not depend on grid resolution (36 / 72 / 144)', () => {
  for (const name of ['E1', 'E4']) {
    const p = ref.optima.find((o) => o.name.startsWith(name)).params, v = [36, 72, 144].map((n) => computePlaneChange(p, MU_M3, { nGrid: n }).totalDeltaV);
    near(v[0], v[1], 1e-6, name); near(v[1], v[2], 1e-6, name);
  }
});

test('runs in well under a few seconds (it needed ~7 s per solve under Pyodide)', () => {
  const p = ref.optima.find((o) => o.name.startsWith('E2')).params, t = performance.now();
  computePlaneChange(p, MU_M3);
  assert.ok(performance.now() - t < 3000, `took ${performance.now() - t} ms`);
});

test('optimiser and Kepler building blocks', () => {
  near(brentMin((x) => (x - 1.234) ** 2 + 5, -10, 10).x, 1.234, 1e-7);
  near(brentMin(Math.cos, 2, 4, 1e-12).x, Math.PI, 1e-6);
  const rosen = (v) => (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2;
  const r = nelderMead(rosen, [[-1.2, 1], [-1, 1], [-1.2, 1.2]], { maxIter: 2000 });
  near(r.x[0], 1, 1e-3); near(r.x[1], 1, 1e-3);
  const at = 10500, vp = Math.sqrt(MU_KM * (2 / 7000 - 1 / at));       // Hohmann ellipse, half period
  near(timeOfFlight(9333.333333333334, Math.PI, [7000, 0, 0], [0, vp, 0], MU_KM), Math.PI * Math.sqrt(at ** 3 / MU_KM), 1e-6);
  near(minRadiusOnArc({ p: 9333.33, e: 1 / 3, nuStart: 0.5, nuEnd: 2 }), 9333.33 / (1 + (1 / 3) * Math.cos(0.5)), 1e-9);   // arc excludes periapsis: min at an end
  near(minRadiusOnArc({ p: 9333.33, e: 1 / 3, nuStart: -0.5, nuEnd: 2 }), 9333.33 / (1 + 1 / 3), 1e-9);                   // arc contains periapsis
});

test('globalSearch exposes the landscape diagnostics', () => {
  const p = ref.optima.find((o) => o.name.startsWith('E4')).params, g = globalSearch(orb(p, 1), orb(p, 2), MU_KM);
  assert.ok(g.basins.length >= 2 && g.nGrid === 72 && g.gridMin >= g.best.dv - 1e-9);
});

test('invalid orbits return an error message instead of NaN', () => {
  const base = { r1a: 7000, r1p: 7000, i1: 0, RAAN1: 0, w1: 0, r2a: 14000, r2p: 14000, i2: 20, RAAN2: 0, w2: 0 };
  assert.match(computePlaneChange({ ...base, r1a: 6000 }, MU_M3).error, /Orbit 1/);
  assert.match(computePlaneChange({ ...base, r2p: -5 }, MU_M3).error, /Orbit 2/);
});
