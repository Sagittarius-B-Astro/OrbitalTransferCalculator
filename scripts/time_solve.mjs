import fs from 'node:fs';
import { computePlaneChange } from '../src/math/planechange/index.js';

const ref = JSON.parse(fs.readFileSync(new URL('../tests/js/fixtures/planechange_reference.json', import.meta.url)));
const MU_M3 = ref.mu_km3 * 1e9;

for (const c of ref.optima) {
  let t = performance.now();
  computePlaneChange(c.params, MU_M3);            // first call includes JIT warm-up
  const cold = performance.now() - t;
  const runs = [];
  for (let i = 0; i < 10; i++) {
    t = performance.now();
    computePlaneChange(c.params, MU_M3);
    runs.push(performance.now() - t);
  }
  runs.sort((a, b) => a - b);
  console.log(c.name.padEnd(32), `cold ${cold.toFixed(0)} ms | median ${runs[5].toFixed(0)} ms | min ${runs[0].toFixed(0)} | max ${runs[9].toFixed(0)}`);
}