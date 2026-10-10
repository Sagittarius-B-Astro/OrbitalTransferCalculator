# Orbital Transfer Calculator

Calculates and visualizes the orbital transfers from Astronautics: **Hohmann**, **bi-elliptic**, **common-apse-line**, **apse-line rotation**, and the minimum-Δv **two-impulse plane change** between two arbitrary orbits. Results are plotted in 3-D with Desmos.

Everything runs **client-side in plain JavaScript** (ES modules, no build step, no server). The plane-change solver was originally written in Python and run in the browser through Pyodide; it has been ported to native JS (see [Why Pyodide was removed](#why-pyodide-was-removed)). The Python version is kept under `reference/python/` as an independent implementation and as the source of the JS test fixtures.

## Repository layout

```
index.html                      page markup (Desmos script in <head>, one ES-module entry point)
OTCstyle.css                    'mission control' theme: CSS grid, glass panels, sticky 3-D plot, responsive
src/
  main.js                       entry point: DOM wiring only
  math/                         pure functions, no DOM (unit-tested in Node)
    common.js  hohmann.js  bielliptic.js  commonApse.js  apseRotate.js  index.js
    planechange/                two-impulse plane change (JS port + improvements)
      index.js                  computePlaneChange(): best of general and nodal searches -> result object
      geometry.js               3-vector helpers, Orbit construction, stateAt, anomalyOf
      pFamily.js                transfers parametrised by semi-latus rectum p; F(nu1,nu2) = min over p and both ways
      search.js                 globalSearch(): grid over (nu1,nu2), Nelder-Mead polish, basin list
      nodal.js                  exact solver for burns at opposite ends of the line of nodes (transfer angle = pi)
      kepler.js                 time of flight, conic description for plotting, min radius on the travelled arc
      optimize.js               Brent minimiser, n-D Nelder-Mead (dependency-free)
  ui/params.js                  ONE schema per transfer type -> renders inputs and reads them back
  ui/results.js                 result box HTML
  plot/desmos.js                one plotter per transfer type
tests/js/                       Node tests (node:test) + fixtures/planechange_reference.json (from the Python reference)
reference/python/               Python reference implementation (numpy only) - NOT used by the web app
  planechange/                  optimizers, kepler, frames, lambert_izzo, pparam, min_dv, nodal
  tests/                        pytest suite for the Python package
  scripts/                      gen_planechange_fixtures.py, validate_claims.py, bench_p_vs_lambert.py
```

## Running it

Open `index.html` through any static file server (ES modules need `http://`, not `file://`), e.g. `python3 -m http.server`. Deployment is any static host (GitHub Pages works as is).

```bash
npm test                        # JS tests (node --test tests/js/*.test.js)
node tests/js/verify_planechange.mjs   # end-to-end plane-change checks against closed forms and Python reference values
# optional, Python reference package
pip install numpy pytest && python3 -m pytest -q reference/python/tests
python3 reference/python/scripts/gen_planechange_fixtures.py   # regenerate tests/js/fixtures/planechange_reference.json (~1 min)
```

## Transfer types

| Type | Inputs | Output |
|---|---|---|
| Hohmann | r₁, r₂ | Δv₁, Δv₂, total, time |
| Bi-elliptic | r₁, r₂, intermediate apoapsis rᵢ | Δv₁–Δv₃, total, time |
| Common apse line | both orbits' apse radii, true anomalies A₁, A₂ | Δv₁, Δv₂, total, flight-path angles, time (Kepler's equation between A₁ and A₂) |
| Apse-line rotation | both orbits' apse radii, rotation η | single-impulse Δv at the orbit intersection, flight-path angle |
| Two-impulse plane change | apoapsis/periapsis radius, i, Ω, ω of both orbits | min total Δv, per-burn Δv, burn true anomalies, transfer p, e, time, minimum radius, search method, list of local-minimum basins |

The plane-change result also reports the minimum radius on the *travelled arc*; the UI warns if the transfer passes inside the central body.

## Two-impulse plane change

### Problem

Given two orbits, find the two-impulse transfer with minimum total Δv and **free** time of flight. A transfer is described by ν₁ (first burn on orbit 1), ν₂ (second burn on orbit 2) and the conic through **r₁(ν₁)** and **r₂(ν₂)**. For fixed endpoints the connecting conics form a one-parameter family (plus the short/long way and the number of revolutions).

### Parametrize by the semi-latus rectum p, not by time of flight

With transfer angle Δθ the terminal velocities follow from Lagrange's coefficients,

$$f = 1-\frac{r_2}{p}(1-\cos\Delta\theta),\quad g=\frac{r_1 r_2\sin\Delta\theta}{\sqrt{\mu p}},\quad \dot g = 1-\frac{r_1}{p}(1-\cos\Delta\theta)$$

$$\mathbf v_1=\frac{\mathbf r_2-f\,\mathbf r_1}{g},\qquad \mathbf v_2=\frac{\dot g\,\mathbf r_2-\mathbf r_1}{g}$$

so Δv(ν₁, ν₂, p) = |v₁ − v꜀₁| + |v꜀₂ − v₂| is an explicit function: no Kepler solve and no Lambert iteration. Properties:

* **Analytic bounds on p.** With k = 2r₁r₂sin²(Δθ/2), l = r₁+r₂, w = 2√(r₁r₂)cos(Δθ/2): pₐ = k/(l+w), p_b = k/(l−w). Short way (Δθ < π): p ∈ (pₐ, ∞), ellipses up to the parabola at p_b, hyperbolas beyond. Long way (Δθ > π): p ∈ (0, pₐ), all ellipses.
* **Revolutions do not change Δv.** An M-revolution solution is the same conic with M extra periods of coasting, so a free-time search never enumerates M.
* **Time of flight is an output** (Kepler / Barker / hyperbolic form).
* **Cost.** In the Python reference, closed-form velocities took ≈ 40 µs/call vs ≈ 187 µs for an Izzo Lambert solve (≈ 4.7×; per-call micro-benchmark only).

The reduced problem is a 2-D function on the torus,

$$F(\nu_1,\nu_2)=\min_{\text{way}\in\{\text{short},\text{long}\}}\ \min_{p}\ \Delta v(\nu_1,\nu_2,p),$$

where the inner minimum over p is a coarse sweep to bracket followed by Brent's method.

### The nodal slice (Δθ = π)

If burn 1 is at s·d̂ and burn 2 at −s·d̂ on the line of nodes (d̂ = n̂₁×n̂₂/|n̂₁×n̂₂|, s = ±1), then r₂ ∥ −r₁ and the (ν₁, ν₂, p) form is singular. Directly from the orbit equation, p = 2r₁r₂/(r₁+r₂) is fixed, e·cos ν₁ = (r₂−r₁)/(r₁+r₂), and two parameters are free: q = e·sin ν₁ (dimensionless radial velocity at r₁) and φ, the orientation of the transfer plane about the line of nodes:

$$\mathbf v_1=\sqrt{\tfrac{\mu}{p}}\,q\,\hat{\mathbf r}_1+\tfrac{\sqrt{\mu p}}{r_1}\hat{\mathbf t}(\varphi),\qquad \mathbf v_2=\sqrt{\tfrac{\mu}{p}}\,q\,\hat{\mathbf r}_1-\tfrac{\sqrt{\mu p}}{r_2}\hat{\mathbf t}(\varphi)$$

Δv is **convex in q** for fixed (s, φ), so the inner minimization is guaranteed unimodal. Ellipses need |q| < 2√(r₁r₂)/(r₁+r₂); hyperbolas only exist for q < 0. The solver scans 72 values of φ for each s and polishes the two best cells. This solver is exact where the general search is singular, and is the right answer for circular-to-circular problems.

### Search strategy

1. **General search** (`search.js`): `nGrid × nGrid` scan of F (default 72 × 72) → local minima of the grid are basin seeds → the best seeds are polished with 2-D Nelder–Mead → de-duplicated basins, sorted by Δv.
2. **Nodal search** (`nodal.js`), run whenever the orbital planes differ.
3. `computePlaneChange` keeps whichever is smaller and reports `method: 'nodal' | 'general'` and the basin list.

Safeguards: both ways are evaluated at every (ν₁, ν₂), so there is no prograde/retrograde discontinuity; cells within about 10⁻⁴ rad of r₁ ∥ r₂ are excluded from the general search (the 1/sin Δθ factor amplifies rounding error until the optimizer exploits noise; the nodal solver owns that region); p is mapped from an open interval to ℝ so the search cannot leave the valid range.

**Why Nelder–Mead** (vs CMA-ES or Lipschitz branch and bound): with p eliminated the problem is only 2-D and F is cheap, so the grid supplies global exploration and what remains is precise local convergence. Nelder–Mead needs no derivatives or smoothness (F is a minimum over branches, so it has kinks), is deterministic, and is a few dozen lines of dependency-free JS. CMA-ES is stochastic and aimed at harder, higher-dimensional landscapes. Lipschitz branch and bound could certify the optimum, but Δv diverges as sin Δθ → 0, so there is no finite global Lipschitz constant, and the excluded neighbourhood is exactly the nodal slice. Neither alternative was implemented or benchmarked here; the LaTeX project write-up has the full discussion. Global optimality is established empirically (grid-size independence, agreement with a high-effort independent Python search), not proved.

### Why Pyodide was removed

The first implementation ran the Python solver in the browser through Pyodide. It was removed because:

* **Too slow.** About 7 s per solve under Pyodide (and about 2–6 s in native CPython with a 12–16 point grid). The JS solver uses a 72 × 72 grid and finishes the six benchmark cases below in roughly 0.1 s each (Node 22, one sandbox CPU; re-measure on your machine).
* **Startup cost.** Pyodide + numpy is a multi-megabyte runtime downloaded and initialized before the first solve.
* **Complexity.** A lazy loader, a virtual file system into which modules were fetched one by one, a JSON bridge, and a path that was never verified in a browser.
* **The port is not a transliteration.** The JS search evaluates both ways, uses a larger grid, enumerates basins, and treats the near-collinear region through the nodal solver. The Python `min_dv.py` keeps only the way that is prograde about orbit 1's normal, which makes F jump by about 1.3 km/s across the branch flip (regression-tested in JS).

## Verification

`node tests/js/verify_planechange.mjs` checks, all passing at the time of writing:

| Check | Result |
|---|---|
| Pure 20° plane change at 7000 km | 2.6207168 km/s = 2v·sin(Δi/2), matches to 1e-9 |
| 7000 → 14000 km circular with 20° plane change | 2.957460 km/s = Hohmann with the plane change optimally split between burns (closed-form reference); whole change at apoapsis would be 3.1071 km/s (4.8 % more); method `nodal`; transfer time = Hohmann half period |
| Coplanar circular 7000 → 14000 km | 2.146528 km/s = Hohmann, method `general` |
| Six reference cases (circular, LEO→GEO 28.5°, E1–E4 elliptic/Molniya) vs. the high-effort Python search | JS ≤ Python in every case; E1–E4 agree to ≥ 6 digits |
| F(ν₁, ν₂) at three random points vs Python | agree to ≈ 1e-9 |
| Resulting conic passes through both burn points; Δv₁ + Δv₂ = total | end-point errors < 1e-10 km |
| Grid 36 / 72 / 144 on two elliptic cases | identical to 8 digits |

Note on the two near-nodal reference cases (circular → circular, LEO → GEO): the Python general search reports values about 7e-6 and 8e-5 km/s *below* the exact nodal optimum. Those are rounding-noise artefacts of the 1/sin Δθ singularity, not real improvements (the closed-form value for the first case is 2.957460), which is why the JS tests allow 2e-4 km/s slack on cases flagged `nearNodal` in the fixture.

The Python reference suite additionally checks the p-family against Izzo's Lambert solver (agreement ≈ 5e-12 km/s over 180 cases), Izzo solutions against an independent universal-variable propagator, and the worked example of Blanco (2025). Run `reference/python/scripts/validate_claims.py` to reproduce those numbers.

## Known limitations

1. **Two impulses are not always optimal.** For large plane changes or large radius ratios a three-impulse (bi-elliptic-style) solution can be cheaper. Lawden's primer-vector conditions would certify optimality or show when a mid-course impulse helps; not implemented.
2. **Nelder–Mead / grid is not a global guarantee.** A very narrow basin could be missed by the grid; the returned basin list shows when several local minima exist.
3. **Free time only.** There is no maximum or fixed time-of-flight option yet.
4. **Planet radius** is checked against the minimum radius on the transfer arc only for plane change; the other transfer types use the radius for drawing.
5. **Apse-line rotation** is a single-impulse intersection transfer; the plot shows the two ellipses but not the burn point.

## Roadmap

1. "Max time of flight" input → one-sided bound on p (p ≥ p\* for the short way, from one M = 0 Lambert solve per cell).
2. Fixed-time-of-flight mode (needs a JS port of the Izzo solver; the Python one is in `reference/python/planechange/lambert_izzo.py`).
3. Primer-vector check / third impulse.
4. Hybrid certificate: Lipschitz branch and bound away from the singular set to certify the grid result, Nelder–Mead to polish.
5. Show the basin list in the UI and let the user choose a basin to plot.

## References

* D. Izzo, *Revisiting Lambert's problem*, Celest. Mech. Dyn. Astron. 121, 1–15 (2015), doi:10.1007/s10569-014-9587-y.
* R. H. Gooding, *A procedure for the solution of Lambert's orbital boundary-value problem*, Celest. Mech. Dyn. Astron. 48, 145–165 (1990).
* P. R. Blanco, *Connect the dots… finding all possible orbits between two points*, Eur. J. Phys. 46, 045004 (2025), doi:10.1088/1361-6404/ade37d, [arXiv:2508.02695](https://arxiv.org/abs/2508.02695): p-parametrization, bounds on p, ΔT(p).
* *Optimal switching conditions for minimum fuel fixed time transfer between non coplanar elliptical orbits*, Acta Astronautica 11(10–11), 621–631 (1984) ([abstract](https://deepblue.lib.umich.edu/items/39c4ca47-898f-4f88-a655-4092df3dd9db)).
* *Optimal bi-impulse orbital transfer between coplanar orbits*, EUCASS 2017-151 ([pdf](https://www.eucass.eu/doi/EUCASS2017-151.pdf)); excerpts only.
* J. A. Nelder and R. Mead, *A simplex method for function minimization*, Comput. J. 7(4), 308–313 (1965); R. P. Brent, *Algorithms for Minimization without Derivatives* (1973).
* H. D. Curtis, *Orbital Mechanics for Engineering Students*; R. R. Bate, D. D. Mueller, J. E. White, *Fundamentals of Astrodynamics*; D. F. Lawden, *Optimal Trajectories for Space Navigation*.
