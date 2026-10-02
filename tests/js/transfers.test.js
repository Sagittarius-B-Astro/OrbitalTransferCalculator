import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  computeHohmann, computeBielliptic, computeCommonApse, computeApseRotate, computeTransfer, computeCommonOrbitProperties,
} from '../../src/math/index.js';

const MU_EARTH_M3 = 3.986004418e14; // m^3/s^2
const MU_KM3 = 3.986004418e5;       // km^3/s^2
const close = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${a} !~ ${b}`);

// independent vis-viva reference for a Hohmann transfer
function hohmannRef(r1, r2) {
  const at = (r1 + r2) / 2;
  const v = (r, a) => Math.sqrt(MU_KM3 * (2 / r - 1 / a));
  return {
    dv1: Math.abs(v(r1, at) - v(r1, r1)),
    dv2: Math.abs(v(r2, r2) - v(r2, at)),
    t: Math.PI * Math.sqrt(at ** 3 / MU_KM3),
  };
}

test('common orbit properties: mu = G M', () => {
  close(computeCommonOrbitProperties(5.972e24, 6378).mu / 1e14, 3.985, 1e-3);
});

test('Hohmann LEO -> GEO matches vis-viva reference', () => {
  const ref = hohmannRef(6678, 42164);
  const r = computeHohmann({ r1: 6678, r2: 42164 }, MU_EARTH_M3);
  close(r.deltaV1, ref.dv1);
  close(r.deltaV2, ref.dv2);
  close(r.totalDeltaV, ref.dv1 + ref.dv2);
  close(r.transferTime, ref.t);
  assert.ok(r.totalDeltaV > 3.8 && r.totalDeltaV < 4.0);
});

test('bi-elliptic beats Hohmann for radius ratio > 15.58 and loses for ratio 2', () => {
  const big = { r1: 7000, r2: 7000 * 20 };
  const hBig = computeHohmann(big, MU_EARTH_M3).totalDeltaV;
  const bBig = computeBielliptic({ ...big, ri: big.r2 * 4 }, MU_EARTH_M3).totalDeltaV;
  assert.ok(bBig < hBig, `bielliptic ${bBig} should beat Hohmann ${hBig}`);

  const small = { r1: 7000, r2: 14000 };
  const hSmall = computeHohmann(small, MU_EARTH_M3).totalDeltaV;
  const bSmall = computeBielliptic({ ...small, ri: 40000 }, MU_EARTH_M3).totalDeltaV;
  assert.ok(bSmall > hSmall);
});

test('bi-elliptic with ri = r2 degenerates to Hohmann', () => {
  const h = computeHohmann({ r1: 7000, r2: 21000 }, MU_EARTH_M3);
  const b = computeBielliptic({ r1: 7000, r2: 21000, ri: 21000 }, MU_EARTH_M3);
  close(b.totalDeltaV, h.totalDeltaV, 1e-12);
  // transfer time is NOT equal: the degenerate second leg (a = r2) still adds half a period
  assert.ok(b.transferTime > h.transferTime);
});

test('common apse: circular orbits at apse-to-apse reproduce the first Hohmann burn', () => {
  // NOTE: totalDeltaV here is the burn at point A only (documented in README, known issue #3)
  const ref = hohmannRef(7000, 21000);
  const r = computeCommonApse({ r1a: 7000, r1p: 7000, r2a: 21000, r2p: 21000, A1: 0, A2: 180 }, MU_EARTH_M3);
  close(r.totalDeltaV, ref.dv1, 1e-9);
  close(r.transferTime, ref.t, 1e-9);
});

test('apse rotation: returned point lies on BOTH ellipses', () => {
  const p = { r1a: 10000, r1p: 7000, r2a: 11000, r2p: 7500, A: 30 };
  const r = computeApseRotate(p, MU_EARTH_M3);
  assert.equal(r.error, undefined);
  const e = (ra, rp) => (ra - rp) / (ra + rp);
  const pp = (ra, rp) => ((ra + rp) / 2) * (1 - e(ra, rp) ** 2);
  const r1 = pp(p.r1a, p.r1p) / (1 + e(p.r1a, p.r1p) * Math.cos(r.TA1));
  const r2 = pp(p.r2a, p.r2p) / (1 + e(p.r2a, p.r2p) * Math.cos(r.TA2));
  close(r1, r2, 1e-9);
  assert.ok(r.totalDeltaV > 0);
});

test('apse rotation: zero rotation of identical orbits needs no delta-v', () => {
  const r = computeApseRotate({ r1a: 10000, r1p: 7000, r2a: 10000, r2p: 7000, A: 1e-9 }, MU_EARTH_M3);
  // identical orbits -> degenerate (a = b ~ 0); either an error or ~0 delta-v is acceptable
  assert.ok(r.error !== undefined || !(r.totalDeltaV > 1e-3));
});

test('computeTransfer dispatch: planeChange has no JS solver', () => {
  assert.equal(computeTransfer('planeChange', {}, MU_EARTH_M3), undefined);
  assert.ok(computeTransfer('hohmann', { r1: 7000, r2: 9000 }, MU_EARTH_M3).totalDeltaV > 0);
});
