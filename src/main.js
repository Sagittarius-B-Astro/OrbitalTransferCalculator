import { computeCommonOrbitProperties, computeTransfer, M3_TO_KM3 } from './math/index.js';
import { renderTypeParams, readTypeParams, hasMissing } from './ui/params.js';
import { renderResults } from './ui/results.js';
import { transferExpressions } from './plot/desmos.js';

document.addEventListener('DOMContentLoaded', () => {
  const calculator = Desmos.Calculator3D(document.getElementById('calculator'), { expressionsCollapsed: true });

  const typeSelect = document.querySelector('#transferType');
  const paramsDiv = document.querySelector('.parametersDynamic');
  const infoBox = document.querySelector('#infobox');

  renderTypeParams(paramsDiv, typeSelect.value);
  typeSelect.addEventListener('change', (e) => renderTypeParams(paramsDiv, e.target.value));

  document.querySelector('#calculateBtn').addEventListener('click', async () => {
    const M = parseFloat(document.querySelector('#mass').value);
    const rad = parseFloat(document.querySelector('#object_radius').value);
    const type = typeSelect.value;
    const params = readTypeParams(type);

    if (Number.isNaN(M) || Number.isNaN(rad) || hasMissing(params)) {
      infoBox.innerHTML = 'Parameter(s) missing! Please fill in missing values';
      return;
    }

    const { mu } = computeCommonOrbitProperties(M, rad); // m^3/s^2
    let result;
    if (type === 'planeChange') {
      infoBox.innerHTML = 'Computing… (first run downloads Python, may take a while)';
      try {
        const { runPlaneChange } = await import('./pyodide/loader.js');
        result = await runPlaneChange({
          ...params,
          mu: mu * M3_TO_KM3,
          n_grid: 24,
        });
      } catch (err) {
        infoBox.innerHTML = `Plane change failed: ${err.message}`;
        return;
      }
    } else {
      result = computeTransfer(type, params, mu);
    }

    if (result) infoBox.innerHTML = renderResults(result);

    calculator.setExpression({ id: 'planet', latex: `x^2 + y^2 + z^2 = (${rad})^2`, color: '#88aaff' });
    transferExpressions(calculator, type, params, result);
  });
});
