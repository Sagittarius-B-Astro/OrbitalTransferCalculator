import { computeHohmann } from './hohmann.js';
import { computeBielliptic } from './bielliptic.js';
import { computeCommonApse } from './commonApse.js';
import { computeApseRotate } from './apseRotate.js';
import { computePlaneChange } from './planechange/index.js';

export * from './common.js';
export { computeHohmann, computeBielliptic, computeCommonApse, computeApseRotate };

const SOLVERS = {
  hohmann: computeHohmann,
  bielliptic: computeBielliptic,
  commonApse: computeCommonApse,
  apseRotate: computeApseRotate,
  planeChange: computePlaneChange,
};

/** Dispatch by transfer type. Returns undefined for types without a JS solver (e.g. planeChange). */
export function computeTransfer(type, params, muM3) {
  const solver = SOLVERS[type];
  return solver ? solver(params, muM3) : undefined;
}
