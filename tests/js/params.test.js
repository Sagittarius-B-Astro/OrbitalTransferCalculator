import { test } from 'node:test';
import assert from 'node:assert/strict';
import { TRANSFER_SCHEMAS, hasMissing } from '../../src/ui/params.js';
import { renderResults } from '../../src/ui/results.js';

test('every schema has unique DOM ids and keys', () => {
  for (const [name, s] of Object.entries(TRANSFER_SCHEMAS)) {
    const ids = s.fields.map((f) => f.id);
    const keys = s.fields.map((f) => f.key);
    assert.equal(new Set(ids).size, ids.length, `${name}: duplicate ids`);
    assert.equal(new Set(keys).size, keys.length, `${name}: duplicate keys`);
  }
});

test('plane-change schema carries the keys planechange.main.solve expects', () => {
  const keys = TRANSFER_SCHEMAS.planeChange.fields.map((f) => f.key).sort();
  assert.deepEqual(keys, ['RAAN1', 'RAAN2', 'i1', 'i2', 'r1a', 'r1p', 'r2a', 'r2p', 'w1', 'w2'].sort());
});

test('hasMissing detects NaN', () => {
  assert.equal(hasMissing({ a: 1, b: NaN }), true);
  assert.equal(hasMissing({ a: 1, b: 2 }), false);
});

test('results show a zero delta-v instead of hiding it', () => {
  assert.match(renderResults({ totalDeltaV: 0 }), /0\.000 km\/s/);
});
