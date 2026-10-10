import { computeCommonOrbitProperties, computeTransfer } from './math/index.js';
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
    const result = computeTransfer(type, params, mu);

    if (result) infoBox.innerHTML = renderResults(result);
    if (result?.minRadius < rad) infoBox.innerHTML += `<br>⚠ The transfer passes ${result.minRadius.toFixed(0)} km from the centre: inside the planet.`;

    calculator.setBlank(); 
    calculator.setExpression({ id: 'planet', latex: `x^2 + y^2 + z^2 = (${rad})^2`, color: '#88aaff' });
    transferExpressions(calculator, type, params, result);
  });
});
