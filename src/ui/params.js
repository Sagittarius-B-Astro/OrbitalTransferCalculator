/**
 * Single source of truth for each transfer type's input fields.
 * `key` is the name the math layer expects; `id` is the DOM id.
 */
const ORBIT_APSES = [
  { key: 'r1a', id: 'r1a', label: 'Initial Orbit Apoapsis Radius (km)', group: 'Initial orbit' },
  { key: 'r1p', id: 'r1p', label: 'Initial Orbit Periapsis Radius (km)', group: 'Initial orbit' },
  { key: 'r2a', id: 'r2a', label: 'Target Orbit Apoapsis Radius (km)', group: 'Target orbit' },
  { key: 'r2p', id: 'r2p', label: 'Target Orbit Periapsis Radius (km)', group: 'Target orbit' },
];

export const TRANSFER_SCHEMAS = {
  hohmann: {
    title: 'Hohmann Transfer Parameters',
    fields: [
      { key: 'r1', id: 'r1', label: 'Initial Orbit Radius (km)' },
      { key: 'r2', id: 'r2', label: 'Target Orbit Radius (km)' },
    ],
  },
  bielliptic: {
    title: 'Bi-Elliptic Hohmann Transfer Parameters',
    fields: [
      { key: 'r1', id: 'r1', label: 'Initial Orbit Radius (km)' },
      { key: 'r2', id: 'r2', label: 'Target Orbit Radius (km)' },
      { key: 'ri', id: 'ri', label: 'Intermediate Orbit Apogee Radius (km)' },
    ],
  },
  commonApse: {
    title: 'Elliptic Transfer on Common Apse Line Parameters',
    fields: [
      ...ORBIT_APSES,
      { key: 'A1', id: 'ta1', label: 'Initial True Anomaly (degrees)', group: 'Burn points' },
      { key: 'A2', id: 'ta2', label: 'Target True Anomaly (degrees)', group: 'Burn points' },
    ],
  },
  apseRotate: {
    title: 'Elliptic Transfer with Rotated Apse Line Parameters',
    fields: [...ORBIT_APSES, { key: 'A', id: 'eta', label: 'Apse Line Rotation (degrees)', group: 'Rotation' }],
  },
  planeChange: {
    title: 'Minimum Delta V Two-ImpulsePlane Change Transfer Parameters',
    fields: [
      ...ORBIT_APSES.slice(0, 2),
      { key: 'i1', id: 'i1', label: 'Initial Inclination (degrees)', group: 'Initial orbit' },
      { key: 'RAAN1', id: 'raan1', label: 'Initial RAAN (degrees)', group: 'Initial orbit' },
      { key: 'w1', id: 'w1', label: 'Initial Argument of Periapsis (degrees)', group: 'Initial orbit' },
      ...ORBIT_APSES.slice(2),
      { key: 'i2', id: 'i2', label: 'Target Inclination (degrees)', group: 'Target orbit' },
      { key: 'RAAN2', id: 'raan2', label: 'Target RAAN (degrees)', group: 'Target orbit' },
      { key: 'w2', id: 'w2', label: 'Target Argument of Periapsis (degrees)', group: 'Target orbit' },
    ],
  },
};

/** Render the dynamic parameter inputs for `type` into `container`. */
export function renderTypeParams(container, type) {
  const schema = TRANSFER_SCHEMAS[type];
  if (!schema) {
    container.innerHTML = '';
    return;
  }
  let lastGroup = null;
  const rows = schema.fields
    .map((f) => {
      const heading = f.group && f.group !== lastGroup ? `<h4 class="group-title">${f.group}</h4>\n` : '';
      lastGroup = f.group ?? lastGroup;
      return `${heading}<div class="field"><label for="${f.id}">${f.label}</label><input type="number" id="${f.id}" name="${f.id}" step="any"></div>`;
    })
    .join('\n');
  container.innerHTML = `<h3>${schema.title}</h3>\n${rows}`;
}

/** Read the inputs for `type` from the DOM. Missing/invalid values come back as NaN. */
export function readTypeParams(type, root = document) {
  const schema = TRANSFER_SCHEMAS[type];
  const params = {};
  for (const f of schema.fields) params[f.key] = parseFloat(root.querySelector(`#${f.id}`).value);
  return params;
}

export const hasMissing = (obj) => Object.values(obj).some((v) => Number.isNaN(v));
