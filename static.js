let typeParamsDiv;

document.addEventListener("DOMContentLoaded", () => {
  var elt = document.getElementById('calculator');
  var calculator = Desmos.Calculator3D(elt, { expressionsCollapsed: true });

  const transferTypeSelect = document.querySelector('#transferType');
  typeParamsDiv = document.querySelector('.parametersDynamic');
  const infoBox = document.querySelector('#infobox');

  // Initialize and listen for dropdown changes
  updateTypeParams(transferTypeSelect.value);
  transferTypeSelect.addEventListener('change', (e) => updateTypeParams(e.target.value));
  
  document.querySelector('#calculateBtn').addEventListener('click', () => {
    const M = parseFloat(document.querySelector('#mass').value);
    const rad = parseFloat(document.querySelector('#object_radius').value);
    const transferType = transferTypeSelect.value;

    const common = computeCommonOrbitProperties(M, rad);

    let invalid;
    let result;
    let params;
    if (transferType === "hohmann") {
      params = {
        r1: parseFloat(document.querySelector('#r1').value),
        r2: parseFloat(document.querySelector('#r2').value),
      };
      result = computeHohmann(params, common.mu);
    } else if (transferType === "bielliptic") {
      params = {
        r1: parseFloat(document.querySelector('#r1').value),
        r2: parseFloat(document.querySelector('#r2').value),
        ri: parseFloat(document.querySelector('#ri').value),
      };
      result = computeBielliptic(params, common.mu);
    } else if (transferType === "commonApse") {
      params = {
        r1a: parseFloat(document.querySelector('#r1a').value),
        r1p: parseFloat(document.querySelector('#r1p').value),
        r2a: parseFloat(document.querySelector('#r2a').value),
        r2p: parseFloat(document.querySelector('#r2p').value),
        A1: parseFloat(document.querySelector('#ta1').value),
        A2: parseFloat(document.querySelector('#ta2').value),
      };
      result = computeCommonApse(params, common.mu);
    } else if (transferType === "apseRotate") {
      params = {
        r1a: parseFloat(document.querySelector('#r1a').value),
        r1p: parseFloat(document.querySelector('#r1p').value),
        r2a: parseFloat(document.querySelector('#r2a').value),
        r2p: parseFloat(document.querySelector('#r2p').value),
        A: parseFloat(document.querySelector('#eta').value),
      };
      result = computeApseRotate(params, common.mu);
    }
    
    for (const key in params) {
      if (isNaN(params[key])) {
        console.log(key)
        invalid=1
      }
    }

    if (isNaN(M) || isNaN(rad) || invalid) {
      infoBox.innerHTML = "Parameter(s) missing! Please fill in missing values";
      return;
    }

    if (result) {
      let html = `<h3>Results:</h3>`;
      if (result.totalDeltaV) html += `Δv (total): ${result.totalDeltaV.toFixed(3)} km/s<br>`;
      if (result.deltaV1) html += `Δv₁: ${result.deltaV1.toFixed(3)} km/s<br>`;
      if (result.deltaV2) html += `Δv₂: ${result.deltaV2.toFixed(3)} km/s<br>`;
      if (result.deltaV3) html += `Δv₃: ${result.deltaV3.toFixed(3)} km/s<br>`;
      if (result.gamma) html += `γ: ${(result.gamma * 180 / Math.PI).toFixed(3)}°<br>`;
      if (result.transferTime) html += `Transfer Time: ${(result.transferTime / 3600).toFixed(3)} hours<br>`;

      infoBox.innerHTML = html;
    }
    
    calculator.setExpression({ id: 'planet', latex: `x^2 + y^2 + z^2 = (${rad})^2`, color: '#88aaff' });
    transferExpressions(transferType);

  });
});

function updateTypeParams(type) {
    if (type === "hohmann") {
        typeParamsDiv.innerHTML = `
        <h3>Hohmann Transfer Parameters</h3>
        <label for="r1">Initial Orbit Radius (km):</label>
        <input type="number" id="r1" name="r1"><br>
        <label for="r2">Target Orbit Radius (km):</label>
        <input type="number" id="r2" name="r2"><br>
        `;
    } else if (type === "bielliptic") {
        typeParamsDiv.innerHTML = `
        <h3>Bi-Elliptic Hohmann Transfer Parameters</h3>
        <label for="r1">Initial Orbit Radius (km):</label>
        <input type="number" id="r1" name="r1"><br>
        <label for="r2">Target Orbit Radius (km):</label>
        <input type="number" id="r2" name="r2"><br>
        <label for="ri">Intermediate Orbit Apogee Radius (km):</label>
        <input type="number" id="ri" name="ri"><br>
        `;
    } else if (type === "commonApse") {
        typeParamsDiv.innerHTML = `
        <h3>Elliptic Transfer on Common Apse Line Parameters</h3>
        <label for="r1a">Initial Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r1a" name="r1a"><br>
        <label for="r1p">Initial Orbit Periapsis Radius (km):</label>
        <input type="number" id="r1p" name="r1p"><br>
        <label for="r2a">Target Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r2a" name="r2a"><br>
        <label for="r2p">Target Orbit Periapsis Radius (km):</label>
        <input type="number" id="r2p" name="r2p"><br>
        <label for="ta1">Initial True Anomaly (degrees):</label>
        <input type="number" id="ta1" name="ta1"><br>
        <label for="ta2">Target True Anomaly (degrees):</label>
        <input type="number" id="ta2" name="ta2"><br>
        `;
    } else if (type === "apseRotate") {
        typeParamsDiv.innerHTML = `
        <h3>Elliptic Transfer with Rotated Apse Line Parameters</h3>
        <label for="r1a">Initial Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r1a" name="r1a"><br>
        <label for="r1p">Initial Orbit Periapsis Radius (km):</label>
        <input type="number" id="r1p" name="r1p"><br>
        <label for="r2a">Target Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r2a" name="r2a"><br>
        <label for="r2p">Target Orbit Periapsis Radius (km):</label>
        <input type="number" id="r2p" name="r2p"><br>
        <label for="eta">Apse Line Rotation (degrees):</label>
        <input type="number" id="eta" name="eta"><br>
        `;
    } else if (type === "planeChange") {
        typeParamsDiv.innerHTML = `
        <h3>Minimum Delta V Plane Change Transfer Parameters</h3>
        <label for="r1a">Initial Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r1a" name="r1a"><br>
        <label for="r1p">Initial Orbit Periapsis Radius (km):</label>
        <input type="number" id="r1p" name="r1p"><br>
        <label for="r2a">Target Orbit Apoapsis Radius (km):</label>
        <input type="number" id="r2a" name="r2a"><br>
        <label for="r2p">Target Orbit Periapsis Radius (km):</label>
        <input type="number" id="r2p" name="r2p"><br>
        <label for="eta">Apse Line Rotation (degrees):</label>
        <input type="number" id="eta" name="eta"><br>
        `;
    } else {
        typeParamsDiv.innerHTML = "";
    }
}
