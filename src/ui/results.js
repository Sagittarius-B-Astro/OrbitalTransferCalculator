const fmt = (x, d = 3) => x.toFixed(d);

/** Build the HTML for the result box. Zero-valued deltas are shown (the original hid them). */
export function renderResults(r) {
  if (r.error) return `<h3>Results:</h3>${r.error}`;
  const has = (k) => typeof r[k] === 'number' && !Number.isNaN(r[k]);
  let html = '<h3>Results:</h3>';
  if (has('totalDeltaV')) html += `Δv (total): ${fmt(r.totalDeltaV)} km/s<br>`;
  if (has('deltaV1')) html += `Δv₁: ${fmt(r.deltaV1)} km/s<br>`;
  if (has('deltaV2')) html += `Δv₂: ${fmt(r.deltaV2)} km/s<br>`;
  if (has('deltaV3')) html += `Δv₃: ${fmt(r.deltaV3)} km/s<br>`;
  if (has('gamma')) html += `γ: ${fmt((r.gamma * 180) / Math.PI)}°<br>`;
  if (has('gamma1')) html += `γ₁ (burn at A): ${fmt((r.gamma1 * 180) / Math.PI)}°<br>`;
  if (has('gamma2')) html += `γ₂ (burn at B): ${fmt((r.gamma2 * 180) / Math.PI)}°<br>`; 
  if (has('nu1Deg')) html += `Impulse 1 true anomaly: ${fmt(r.nu1Deg, 2)}°<br>`;
  if (has('nu2Deg')) html += `Impulse 2 true anomaly: ${fmt(r.nu2Deg, 2)}°<br>`;
  if (has('p')) html += `Transfer semi-latus rectum: ${fmt(r.p, 1)} km<br>`;
  if (has('transferTime')) html += `Transfer Time: ${fmt(r.transferTime / 3600)} hours<br>`;
  return html;
}
