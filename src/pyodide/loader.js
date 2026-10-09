/** Lazy Pyodide loader for the plane-change solver (planechange/*.py). Only fetched when needed. */
const PYODIDE_BASE = 'https://cdn.jsdelivr.net/pyodide/v0.26.0/full/';
const PY_MODULES = ['__init__', 'optimizers', 'kepler', 'frames', 'lambert_izzo', 'pparam', 'min_dv', 'nodal','main'];

let pyPromise = null;

function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement('script');
    s.src = src;
    s.onload = resolve;
    s.onerror = () => reject(new Error(`failed to load ${src}`));
    document.head.appendChild(s);
  });
}

export function loadPlaneChange() {
  if (!pyPromise) {
    pyPromise = (async () => {
      await loadScript(`${PYODIDE_BASE}pyodide.js`);
      const py = await globalThis.loadPyodide({ indexURL: PYODIDE_BASE });
      await py.loadPackage('numpy');
      try { py.FS.mkdir('/planechange'); } catch (_) { /* already exists */ }
      for (const m of PY_MODULES) {
        const res = await fetch(`planechange/${m}.py`);
        if (!res.ok) throw new Error(`could not fetch planechange/${m}.py (${res.status})`);
        py.FS.writeFile(`/planechange/${m}.py`, await res.text());
      }
      py.runPython("import sys\nif '/' not in sys.path: sys.path.insert(0, '/')");
      return py;
    })().catch((err) => { pyPromise = null; throw err; });
  }
  return pyPromise;
}

/** params: km, degrees, mu in km^3/s^2. Returns plain JSON from planechange.main.solve. */
export async function runPlaneChange(params) {
  const py = await loadPlaneChange();
  py.globals.set('params_json', JSON.stringify(params));
  const out = py.runPython(
    'import json\nfrom planechange.main import solve\njson.dumps(solve(json.loads(params_json)))'
  );
  return JSON.parse(out);
}
