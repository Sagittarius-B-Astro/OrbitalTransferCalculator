# Orbital Transfer Calculator

Calculates and visualizes the orbital transfers from Astronautics: Hohmann, bi-elliptic, common-apse-line, apse-line rotation, and (in progress) the minimum-Δv **plane change**. Results are plotted in 3-D with Desmos.

* Hohmann, bi-elliptic, common apse and apse-rotation run entirely in JavaScript.
* Plane change runs a Python solver (`planechange/`) in the browser through [Pyodide](https://pyodide.org) (numpy only, loaded lazily the first time you pick that transfer).

> **Status of plane change:** the solver is implemented and verified natively (Python 3.12 + numpy, see [Verification](#verification)). The Pyodide wiring and the Desmos arc plot have **not yet been exercised in a browser** — treat that path as experimental.

## Repository layout

```
index.html                       page markup (Desmos script in <head>, one ES-module entry point)
OTCstyle.css                     'mission control' theme: CSS grid layout, glass panels sticky 3-D plot, responsive
src/
  main.js                        entry point: DOM wiring only
  math/                          pure functions, no DOM  (unit-tested in Node)
    common.js  hohmann.js  bielliptic.js  commonApse.js  apseRotate.js  index.js
  ui/params.js                   ONE schema per transfer type -> renders inputs and reads them back
  ui/results.js                  result box HTML
  plot/desmos.js                 one plotter function per transfer type
  pyodide/loader.js              lazy Pyodide loader + runPlaneChange()
planechange/                     plane-change solver (Python, numpy only)
  optimizers.py                  Halley, Householder-3, Brent root, golden-section, n-D Nelder-Mead
  kepler.py                      universal-variable propagator (independent check + arc sampling)
  frames.py                      Orbit dataclass, perifocal -> ECI
  lambert_izzo.py                Izzo (2015) Lambert solver, multi-revolution
  pparam.py                      transfers parametrised by semi-latus rectum p (closed form)
  min_dv.py                      free-time minimum-Δv search over (ν1, ν2, p)
  nodal.py                       Δθ = π slice (burns at opposite ends of the line of nodes): (φ, q) parametrization
  main.py                        solve(params) -> JSON for the web app
tests/js, tests/py               58 Python + 12 JS tests
scripts/                         bench_p_vs_lambert.py, validate_claims.py (reproduce every number below)
```

Everything that used to live in `static.js`, `orbitMath.js`, `desmosSetup.js`, `PlaneChange.js` and the single 400-line `PlaneChange.py` is now in the files above; those five old files can be deleted.

---

## Plane change: problem statement

Given two orbits (apoapsis/periapsis radii, inclination, RAAN, argument of periapsis), find the two-impulse transfer with minimum total Δv. A two-impulse transfer is fully described by

* ν₁ – true anomaly of the first burn on orbit 1,
* ν₂ – true anomaly of the second burn on orbit 2,
* the transfer conic through the two resulting points **r₁(ν₁), r₂(ν₂)**.

For fixed r₁, r₂ the connecting conics form a **one-parameter family** (plus discrete choices: short/long way, number of revolutions). The question is how to label that one parameter. Lambert's problem labels it with the time of flight; there is a second natural label, the semi-latus rectum.

---

## Research: constraining the Lambert problem for minimum Δv

### 1. Two ways to close the problem

| | Constrain **ΔT** (Lambert) | Constrain **p** (semi-latus rectum) |
|---|---|---|
| Unknowns searched | ν₁, ν₂ (+ ΔT) | ν₁, ν₂, p |
| Per-point work | Solve Kepler/Lambert iteration | Closed-form (Lagrange coefficients) |
| Natural for | Fixed-time / rendezvous transfers | Free-time orbit-to-orbit transfers |
| Search domain of the 3rd variable | ΔT ∈ (0, ∞), no natural scale | p ∈ (p_lo, p_hi), explicit analytic bounds |
| Revolutions M | Must be enumerated (Mmax) | Irrelevant to Δv (see §3) |

Both describe the same set of conics, so the minimum Δv is the same; they differ in cost and in how awkward the bounds are. (A third equivalent label is Izzo's Lancaster–Blanchard variable **x**: x = 0 is the minimum-energy ellipse, x = 1 the parabola, and Izzo's velocity reconstruction takes (x, y) directly — the root-find that Izzo/Gooding spend their effort on is only the map T → x. Using x as the search variable would skip that map in the same way p does.)

### 2. Option A — constrain ΔT (what the original code did)

The original design samples ΔT on a coarse grid, calls the Izzo solver, then refines. Things that fall out of the papers:

* **Izzo's variables.** Lambert problems with the same c/s are "L-similar" (Gooding 1990); the solution depends on the geometry only through λ² = 1 − c/s, and the dimensionless time is T = √(2μ/s³)·Δt. Two landmarks are known in closed form (Izzo Eqs. 19, 21): the minimum-energy time T₀₀ = arccos λ + λ√(1−λ²) (x = 0) and the parabolic time T₁ = ⅔(1−λ³) (x = 1). Those, not a heuristic "1e-5 … outer period" window, are the natural scale for a ΔT grid.
* **Mmax.** Izzo computes M_max from T (⌊T/π⌋, corrected with T_min(M) via Halley iteration). It is a property of the *requested* ΔT, so it cannot be used to *choose* the ΔT range without circularity (the concern in the old code comments). Separately, the user-facing "max revolutions" is just a cap passed to the solver (`lambert(..., max_revs=)`).
* **Cost.** Every grid cell needs ≥ 1 Lambert solve per sampled ΔT per M, and the old plan also ran a 1-D refinement inside every cell.
* **Bracketing.** Δv(ΔT) can have several local minima (one family per M, two branches per M ≥ 1), so a coarse sample + local refine can land in the wrong basin. The refinement must be a *minimizer* (golden-section / Brent-minimize); Brent's method as implemented is a root finder.

Option A is still the right tool when ΔT is a real constraint (below).

### 3. Option B — constrain p

For fixed endpoints with transfer angle Δθ the velocities on the connecting conic follow from Lagrange's f and g coefficients:

$$f = 1-\frac{r_2}{p}(1-\cos\Delta\theta),\qquad g=\frac{r_1 r_2\sin\Delta\theta}{\sqrt{\mu p}},\qquad \dot g = 1-\frac{r_1}{p}(1-\cos\Delta\theta)$$

$$\mathbf v_1=\frac{\mathbf r_2-f\,\mathbf r_1}{g},\qquad \mathbf v_2=\frac{\dot g\,\mathbf r_2-\mathbf r_1}{g}$$

so **Δv(ν₁, ν₂, p) = |v₁ − v꜀₁| + |v꜀₂ − v₂| is an explicit function** — no Kepler solve, no Lambert iteration (`planechange/pparam.py`). This is the same three-variable formulation used for optimal-switching conditions in the literature (Acta Astronautica 11, 1984; the EUCASS 2017 coplanar bi-impulse paper; see References). Properties worth knowing:

* **Bounds on p are analytic.** With k = 2r₁r₂ sin²(Δθ/2), l = r₁+r₂, w = 2√(r₁r₂) cos(Δθ/2): p_a = k/(l+w), p_b = k/(l−w). For the short way (Δθ < π) valid p ∈ (p_a, ∞): a huge ellipse as p → p_a⁺ (energy → 0⁻, ΔT → ∞; checked numerically), the connecting parabola at p_b, hyperbolas beyond. For the long way (Δθ > π) valid p ∈ (0, p_a), all ellipses, with the parabola at the upper limit. (`p_bounds`)
* **ΔT is monotonic in p** for a given way: decreasing for Δθ < π, increasing for Δθ > π (Blanco 2025, Fig. 4; checked numerically in `test_tof_monotonic_in_p`). So p ↔ ΔT is one-to-one per branch.
* **Revolutions do not affect Δv.** An M-revolution solution is the *same conic* (same p, same velocities) with M extra periods of coasting. In a free-time search M is therefore irrelevant and `Mmax` never needs to be computed — one of the open questions in the old code comments disappears.
* **Time of flight becomes an output.** `time_of_flight` evaluates it from Kepler's equation for elliptic, parabolic (Barker) and hyperbolic cases, optionally adding whole periods.
* **Cost.** Per call, measured with `scripts/bench_p_vs_lambert.py` (best of 5 repeats × 2000 calls, CPython 3.12 / numpy 2.4, one sandbox CPU): Izzo `lambert()` ≈ 187 µs vs `velocities_from_p()` ≈ 40 µs (≈ 4.7×), or ≈ 54 µs (≈ 3.4×) if the time of flight is also evaluated. An earlier single-shot run gave ≈ 6×; that figure was noise-sensitive and should not be quoted. This is a per-call comparison on one geometry; **no end-to-end search comparison was measured** (the original ΔT-grid search never ran, so there is no baseline).
* **A convenient inner problem.** For fixed (ν₁, ν₂), Δv(p) is smooth and, empirically, unimodal (no second interior local minimum in 300 random non-coplanar geometries I tried; not proven). For a circular start orbit the minimum is a root of a quartic (Blanco 2025, §4.3); `min_dv._best_p` just does a coarse sweep + golden section, so the 3-D problem collapses to a 2-D grid over (ν₁, ν₂) with a cheap inner minimization, followed by a 3-D Nelder–Mead polish (p is mapped from an open interval to ℝ so the simplex can never leave the valid range).

### 4. When ΔT really is a constraint

* **Fixed ΔT (rendezvous / phasing):** p is no longer free: p = p(ν₁, ν₂, ΔT) from Lambert (Izzo), and the search is 2-D over (ν₁, ν₂) — but ν₂ is then also tied to the target's phase, so the dimension drops further. Use `lambert_izzo.lambert`.
* **Upper limit ΔT ≤ T_max:** for fixed (ν₁, ν₂) the M = 0 time of flight is monotonic in p, and M = 0 is the *shortest* time on any given conic (extra revolutions only add whole periods). So "ΔT ≤ T_max" is equivalent to p ≥ p*(T_max) for the short way and p ≤ p*(T_max) for the long way, where p* = |r₁ × v₁|²/μ with v₁ from one M = 0 Lambert solve at T_max. Caveat: p* depends on the geometry, so it is a **per-(ν₁, ν₂) bound, not a global one**; it costs one Lambert solve per grid cell and turns the time limit into a box constraint on the inner 1-D problem rather than a new dimension. Verified in `scripts/validate_claims.py` block E.

### 5. Things that bite near the optimum

* **Δθ = π (the nodal case).** See the next subsection; implemented in `planechange/nodal.py`.
* **Two impulses are not always optimal.** For large plane changes or large radius ratios, three-impulse (bi-elliptic-style) solutions can win. Lawden's primer-vector conditions (|primer| = 1 at an impulse, d|primer|/dt = 0 for free-time interior impulses, |primer| ≤ 1 on the arc) give a certificate of optimality and tell you when a mid-course impulse would help; not implemented. (Background knowledge, not from the supplied papers.)
* **Grid resolution.** The coarse grid only needs to land in the right basin; the Nelder–Mead polish does the rest. The three best grid cells are polished, since Δv(ν₁, ν₂) is multi-modal.


### 5a. The nodal case (Δθ = π) parametrization

If burn 1 is at **s·d̂** on orbit 1 and burn 2 at **−s·d̂** on orbit 2 (d̂ = n̂₁ × n̂₂ / |n̂₁ × n̂₂|, s = ±1), then r₂ ∥ −r₁, sin Δθ = 0 and g → 0, so the (ν₁, ν₂, p) form is singular. Work it out directly from the orbit equation with ν₂ = ν₁ + π:

r₁ = p/(1 + e cos ν₁), r₂ = p/(1 − e cos ν₁) ⇒ 1/r₁ + 1/r₂ = 2/p, so

* **p = 2 r₁ r₂ / (r₁ + r₂)** — a single value (no free p),
* **e cos ν₁ = c₀ = (r₂ − r₁)/(r₁ + r₂)** — fixed,
* **e sin ν₁ = q** — *free*: the dimensionless radial velocity at r₁ (v_r1 = √(μ/p)·q), giving e = √(c₀² + q²),
* **φ ∈ [0, 2π)** — *free*: orientation of the transfer plane about the line of nodes. With r̂₁ = r₁/r₁ and t̂(φ) a unit vector ⟂ r̂₁ (t̂ = cos φ ê_a + sin φ ê_b for any orthonormal pair ⟂ r̂₁), the transfer's angular momentum is √(μp)·(r̂₁ × t̂), and

  **v₁ = √(μ/p)·q·r̂₁ + (√(μp)/r₁)·t̂,  v₂ = √(μ/p)·q·r̂₁ − (√(μp)/r₂)·t̂**

  (the radial velocity *vector* is the same at both ends; the transverse components are opposite). φ and φ + π are the two senses of motion, so no separate "retrograde" flag is needed.

Cost: Δv(φ, q) = |v₁ − v_c1| + |v₂ − v_c2|, where v_c1, v_c2 are the circular/orbital velocities of orbits 1 and 2 at the two nodes. Properties:

* **Convex in q** for fixed (s, φ): each term is the norm of a function affine in q. The inner minimization is therefore a guaranteed-unimodal 1-D problem (stronger than the empirical unimodality in p for the general case). At an interior optimum the stationarity condition is (û₁·r̂₁) + (û₂·r̂₁) = 0, where û_i are the unit impulse directions (v_i − v_ci)/|v_i − v_ci| — tested.
* **Feasibility:** ellipse if |q| < q_max = 2√(r₁r₂)/(r₁ + r₂) (parabola at q_max); hyperbolic transfers exist only for **q < 0** (for q > 0 the forward arc from ν₁ to ν₁ + π would cross the hyperbola's asymptote).
* **Time of flight:** Kepler's equation from ν₁ = atan2(q, c₀) to ν₁ + π (`pparam.time_of_flight` with Δθ = π). For q = 0 this is half a period, i.e. the Hohmann-like apse-to-apse case.
* **Search:** 2 values of s × (72-point grid in φ, inner golden-section in q, polish the two best φ cells). `min_delta_v_nodal(o1, o2, mu)`; raises if the planes coincide (then no line of nodes exists and the general search applies).
* **Scope warning:** this is the *Δθ = π slice*, not proof that the global optimum lies on it. It does for circular orbits (verified against the general search, F4 below). For general ellipses with different apse lines the optimum can sit off the nodal line; compare `min_delta_v_nodal` against `min_delta_v` and keep the smaller. It is **not yet wired into `min_dv.py` or the UI**.

### 6. Recommendation

Use **p** as the third variable for the free-time plane-change calculator, **ΔT/Izzo** when the user supplies a time of flight or a time cap, and keep `lambert_izzo.py` for the latter and as an independent cross-check of the p-family (done in the tests). The UI can later expose an optional "max time of flight" box that simply tightens the p bound.

---

## Verification

All of the following run in `tests/py` and pass:

| Check | Result |
|---|---|
| Blanco 2025 worked example (r₁ = 9000 km, r₂ = 15000 km, Δθ = 120°) | p bounds 0.6317 r₁ and 1.8173 r₁; min-Δv orbit at p = 1.3128 r₁ gives Δv = 0.1563 v_c and ΔT = 0.4607 T_c — all reproduced |
| p-family vs Izzo Lambert, short & long way, elliptic & hyperbolic | velocities agree to ~1e-14 km/s in my runs (tests assert 1e-8) |
| Izzo solutions vs independent universal-variable propagator, M = 0, 1, 2, long & short way, 3-D geometry | final position error < 1e-6 km (observed ~1e-10) |
| Coplanar circular 7000 → 14000 km | search returns the Hohmann Δv (2.14653 km/s), a = 10500 km, e = 1/3 |
| Pure 20° plane change at 7000 km | Δv = 2v sin(Δi/2) = 2.62072 km/s, exact to 6 digits |
| 7000 → 14000 km **with** a 20° plane change | 2.9575 km/s, vs 3.1071 km/s for a Hohmann with the whole plane change at apogee (4.8 % less); optimum is the Hohmann ellipse with the plane change shared between burns |

Runtime (native CPython + numpy, laptop-class sandbox): a 12–16 point grid takes about 2–6 s. Pyodide will be slower; I have not measured it. Cheap next steps if it is too slow: vectorize the p sweep over numpy arrays, drop the grid to 8–12 points, or exploit that the first burn can often be fixed at a node.

---

## Reproducing the numbers

```bash
pip install numpy pytest
python3 scripts/bench_p_vs_lambert.py 2000 5   # speed claim (per-call micro-benchmark)
python3 scripts/validate_claims.py             # blocks A–G: every numerical claim, PASS/FAIL (~1 min)
python3 -m pytest -q tests/py                  # 58 unit tests;  node --test tests/js/*.test.js  # 12 JS tests
```

`validate_claims.py` blocks: **A** Blanco worked example · **B** p-family vs Izzo (max |Δv| ≈ 5e-12 km/s over 180 cases) · **C** Izzo vs independent propagator (≈ 3e-9 km over 91 solutions incl. multi-rev) · **D** ΔT(p) monotonic (4 radius ratios × 8 angles × 400 p) · **E** ΔT ≤ T_max ⇔ p ≥ p* · **F** Hohmann / pure plane change / combined burn / nodal-vs-general · **G** single interior minimum of Δv(p) in 300 random geometries (empirical; not a proof). Known wart: block B/`lambert()` can emit a `RuntimeWarning: invalid value in log` from the hyperbolic branch of `x2tof` during intermediate Householder steps; the converged results are unaffected, but the branch should be clamped.

## Changes beyond a pure refactor

The refactor moved code into modules, but I also **changed behaviour** in these places:

* `apseRotate`: three formula/precedence fixes + an "orbits do not intersect" error (see below). Results for this transfer type will differ from before.
* Result box: zero-valued deltas are now shown; the apse-rotation result also carries TA1/TA2/r.
* Inputs are validated before any math runs; unused mathjs script and the dead `computerPlaneChange` stub were removed; the plane-change dropdown label changed from "(coming soon)" to "(experimental)" and now calls the Python solver.
* Desmos: Hohmann / bi-elliptic transfer arcs now get an explicit color; the plane-change plotter is new.
* `PlaneChange.py` was **rewritten, not just split**: Izzo solver re-implemented (hyperbolic branch, Battin series near x = 1, corrected T < T₁ starter), Nelder–Mead generalized to n-D, golden-section added as the actual minimizer, `trajectoryCurve` replaced by sampling the arc with a propagator, and the old ΔT-grid search (`findTOFrange`, `loopOverOrbits`, the per-TOF `minDeltaV`) **replaced** by the (ν₁, ν₂, p) search. The Lambert/ΔT path survives only as `lambert_izzo.py`.
* Math for Hohmann, bi-elliptic and common-apse is numerically unchanged.

---

## Review of the original code

### Fixed in this refactor

JavaScript
* `calculator` was declared with `var` inside the `DOMContentLoaded` callback, but `desmosSetup.js` used it as a global → `ReferenceError` on every plot. It is now passed explicitly.
* `computeApseRotate`: used the angle in **degrees** inside `cos(A)`; used `eta` where `e1` belongs in the radius; `acos(...) % 2π` had the wrong precedence; no check for non-intersecting orbits.
* `if (result.totalDeltaV)` etc. hid legitimate zero values.
* Parameters were validated *after* the math ran (NaN in, NaN out); `planeChange` fell through to the plot with `params` undefined and `r` undefined in its bounds code.
* `PlaneChange.js` attached a handler to `#computeBtn`, which does not exist (the button is `#calculateBtn`) → crash at load. Unused mathjs `<script>` removed. Input HTML is no longer duplicated between `updateTypeParams` and the readers.

Python (`PlaneChange.py`) — the file could not run as written
* Degrees passed to `np.cos/np.sin`; grid indexed with float angles (`grid[TA1][TA2]`); `Rz`/`Rx` defined *after* use; `np.sqrt(mu/p) * [[...]]` (list × scalar).
* `hunit = ..., Lambda = ...` (comma instead of newline), `for x, y in xSols, ySols` (needs `zip`), `solutions[0] = ...` on an empty list, `y = y(x)` shadowing the function, `TOF(x, M=Mmax)` bound before `Mmax` exists, `x0l, x0r = a, \n b` creating a tuple.
* `Halley1d` called with the wrong arity; both Halley and Householder never updated their convergence variable.
* Initial guess for T < T₁ used `− 1`; Izzo's improved guess is `+ 1` (it enforces the parabola value and slope at x = 1).
* Time of flight used only `arccos`, i.e. elliptic; hyperbolic (x > 1) is needed for short flight times.
* `TOFs[TOFi - 1]` wraps around at index 0 and `TOFs[TOFi + 1]` overruns at the end.
* `minDeltaV` is both a function and a local variable inside `minCoords`.
* Brent's method was used as a minimizer; the Brent test bracket (0, 2) is not a bracket (both ends are negative); the Nelder–Mead tests called `f(x, y)` with two arguments while the optimizer passes one array.
* After the fix, `x = 1` exactly gave 0/0 in Battin's series (M = 0) and a singular derivative; both guarded and covered by a continuity test.

### Known issues still open (behaviour intentionally unchanged)

1. `computeCommonApse` returns the Δv of the burn at point A only; I could not find the burn at B in the code, so `totalDeltaV` looks like it is missing a term. The test documents the current behaviour (it equals the first Hohmann burn for circular orbits).
2. `computeCommonApse.transferTime` is `π√(aₜ³/μ)` (half a period), which is only right for an apse-to-apse arc. For general ν₁, ν₂ the time follows from Kepler's equation between the two anomalies (this repo already has that in `pparam.time_of_flight`).
3. The planet radius is only used for drawing; r₁, r₂ are not checked against it.
4. `OTCstyle.css` has `padding-top: -20px` (invalid, ignored by browsers).
5. Plane change: two-impulse only, nodal solver not yet wired in, no fixed-time option in the UI, Pyodide path unverified in a browser.

---

## Roadmap

1. Wire `nodal.py` into the plane-change entry point (take the better of nodal and general search) and add the UI hook.
2. "Max time of flight" input → bound on p; "fixed time" mode using `lambert_izzo`.
3. Primer-vector check / third impulse.
4. Speed: vectorized p sweep; measure under Pyodide.
5. Fix the two commonApse issues once confirmed against the textbook formulas.

---

## References

* D. Izzo, *Revisiting Lambert's problem*, Celestial Mechanics and Dynamical Astronomy 121, 1–15 (2015), doi:10.1007/s10569-014-9587-y — the solver in `lambert_izzo.py` (Eqs. 18–22, 30–31, Algorithms 1–2).
* R. H. Gooding, *A procedure for the solution of Lambert's orbital boundary-value problem*, Celestial Mechanics and Dynamical Astronomy 48, 145–165 (1990) — L-similarity, velocity reconstruction.
* P. R. Blanco, *Connect the dots… finding all possible orbits between two points*, Eur. J. Phys. 46, 045004 (2025), doi:10.1088/1361-6404/ade37d, [arXiv:2508.02695](https://arxiv.org/abs/2508.02695) — p-parametrization, bounds on p, minimum-Δv-from-circular, ΔT(p).
* *Optimal switching conditions for minimum fuel fixed time transfer between non coplanar elliptical orbits*, Acta Astronautica 11(10–11), 621–631 (1984) — transfer elements, Δv and ΔT as explicit functions of (p, ν₁, ν₂); free-time and fixed-time optimality conditions. ([abstract](https://deepblue.lib.umich.edu/items/39c4ca47-898f-4f88-a655-4092df3dd9db))
* *The minimum delta-V Lambert's problem*, Brazilian Society of Automation journal, v7 n2 — [pdf](https://www.sba.org.br/revista/volumes/v7n2/v7n2a04.pdf). Skimmed via search excerpts only, not read in full.
* *Optimal bi-impulse orbital transfer between coplanar orbits*, EUCASS 2017-151 — [pdf](https://www.eucass.eu/doi/EUCASS2017-151.pdf). Optimizes over the two true anomalies and the transfer semi-latus rectum. Excerpts only.
* NASA NTRS 19710029163, *minimal two-impulse orbital transfer* — [pdf](https://ntrs.nasa.gov/api/citations/19710029163/downloads/19710029163.pdf). Table of contents / excerpts only.
* H. D. Curtis, *Orbital Mechanics for Engineering Students* (Lagrange coefficients, universal-variable propagation, apse-line transfers) and D. F. Lawden, *Optimal Trajectories for Space Navigation* (primer vector) — cited from background knowledge.
