/** Desmos 3D plotting for each transfer type. `calc` is the Desmos Calculator3D instance. */
import { deg2rad } from '../math/common.js';

const BLUE = '#00ccff';
const ORANGE = '#ff6600';
const GREEN = '#00ff88';

function setBounds(calc, r) {
  const b = 1.1 * r;
  calc.setMathBounds({ xmin: -b, xmax: b, ymin: -b, ymax: b, zmin: -b, zmax: b });
}

const expr = (calc, id, latex, extra = {}) => calc.setExpression({ id, latex, ...extra });

const curve = (calc, id, latex, color, min, max) =>
  calc.setExpression({ id, expressionType: 'parametric3d', latex, color, parametricDomain: { min, max } });

const circles = (calc) => {
  curve(calc, 'initial_orbit', '(r_1 \\cos(t),r_1 \\sin(t),0)', BLUE, '0', '2 \\pi');
  curve(calc, 'target_orbit', '(r_2 \\cos(t),r_2 \\sin(t),0)', ORANGE, '0', '2 \\pi');
};

/** Set up the two coplanar ellipses (initial/target) with names a_i, b_i, e_i. */
function ellipseVars(calc) {
  expr(calc, 'a1', 'a_1=(r_{1p}+r_{1a})/2');
  expr(calc, 'b1', 'b_1=\\sqrt{r_{1p} r_{1a}}');
  expr(calc, 'e1', 'e_1=(r_{1a}-r_{1p})/(r_{1a}+r_{1p})');
  expr(calc, 'a2', 'a_2=(r_{2p}+r_{2a})/2');
  expr(calc, 'b2', 'b_2=\\sqrt{r_{2p} r_{2a}}');
  expr(calc, 'e2', 'e_2=(r_{2a}-r_{2p})/(r_{2a}+r_{2p})');
}

function apseVars(calc, p) {
  expr(calc, 'r1a', `r_{1a}=${p.r1a}`);
  expr(calc, 'r1p', `r_{1p}=${p.r1p}`);
  expr(calc, 'r2a', `r_{2a}=${p.r2a}`);
  expr(calc, 'r2p', `r_{2p}=${p.r2p}`);
}

const PLOTTERS = {
  hohmann(calc, p) {
    setBounds(calc, Math.max(p.r1, p.r2));
    expr(calc, 'r1', `r_1=${p.r1}`);
    expr(calc, 'r2', `r_2=${p.r2}`);
    expr(calc, 'a', 'a=(r_1+r_2)/2');
    expr(calc, 'b', 'b=\\sqrt{r_1 r_2}');
    expr(calc, 'e', 'e_c=(r_2-r_1)/(r_2+r_1)');
    circles(calc);
    curve(calc, 'trajectory', '(a(\\cos(t)-e_c),b\\sin(t),0)', GREEN, '0', '\\pi');
  },

  bielliptic(calc, p) {
    setBounds(calc, p.ri);
    expr(calc, 'r1', `r_1=${p.r1}`);
    expr(calc, 'r2', `r_2=${p.r2}`);
    expr(calc, 'ri', `r_i=${p.ri}`);
    expr(calc, 'a1', 'a_1=(r_1+r_i)/2');
    expr(calc, 'b1', 'b_1=\\sqrt{r_1 r_i}');
    expr(calc, 'e1', 'e_1=(r_i-r_1)/(r_i+r_1)');
    expr(calc, 'a2', 'a_2=(r_i+r_2)/2');
    expr(calc, 'b2', 'b_2=\\sqrt{r_i r_2}');
    expr(calc, 'e2', 'e_2=(r_i-r_2)/(r_2+r_i)');
    circles(calc);
    curve(calc, 'trajectory1', '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)', GREEN, '0', '\\pi');
    curve(calc, 'trajectory2', '(a_2(\\cos(t)-e_2),b_2\\sin(t),0)', GREEN, '\\pi', '2 \\pi');
  },

  commonApse(calc, p) {
    setBounds(calc, Math.max(p.r1a, p.r2a));
    apseVars(calc, p);
    expr(calc, 'TA1', `T_1=${p.A1}\\pi/180`);
    expr(calc, 'TA2', `T_2=${p.A2}\\pi/180`);
    ellipseVars(calc);
    expr(calc, 'p1', 'p_1=a_1(1-e_1^2)');
    expr(calc, 'p2', 'p_2=a_2(1-e_2^2)');
    expr(calc, 'rA', 'r_A=p_1/(1+e_1\\cos(T_1))');
    expr(calc, 'rB', 'r_B=p_2/(1+e_2\\cos(T_2))');
    expr(calc, 'et', 'e_t=(r_B-r_A)/(r_A\\cos(T_1)-r_B\\cos(T_2))');
    expr(calc, 'pt', 'p_t=r_A r_B(\\cos(T_1)-\\cos(T_2))/(r_A\\cos(T_1)-r_B\\cos(T_2))');
    expr(calc, 'at', 'a_t=p_t/(1-e_t^2)');
    expr(calc, 'bt', 'b_t=a_t\\sqrt{1-e_t^2}');
    const A2u = p.A2 > p.A1 ? p.A2 : p.A2 + 360;
    expr(calc, 'TA2', `T_2=${A2u}\\pi/180`);
    curve(calc, 'initial_orbit', '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)', BLUE, '0', '2 \\pi');
    curve(calc, 'target_orbit', '(a_2(\\cos(t)-e_2),b_2\\sin(t),0)', ORANGE, '0', '2 \\pi');
    curve(calc, 'trajectory', '((p_t\\cos(t))/(1+e_t\\cos(t)),(p_t\\sin(t))/(1+e_t\\cos(t)),0)', GREEN, 'T_1', 'T_2');
  },

  apseRotate(calc, p) {
    setBounds(calc, Math.max(p.r1a, p.r2a));
    apseVars(calc, p);
    expr(calc, 'eta', `\\eta=${p.A}\\pi/180`);
    ellipseVars(calc);
    curve(calc, 'initial_orbit', '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)', BLUE, '0', '2 \\pi');
    curve(
      calc,
      'target_orbit',
      '(((a_2(\\cos(t)-e_2))\\cos(\\eta)-(b_2\\sin(t))\\sin(\\eta)),((a_2(\\cos(t)-e_2))\\sin(\\eta)+(b_2\\sin(t))\\cos(\\eta)),0)',
      ORANGE,
      '0',
      '2 \\pi'
    );
  },

  /** Plane change: both orbits are drawn from numeric rotation matrices; the arc comes from Python. */
  planeChange(calc, p, result) {
    setBounds(calc, Math.max(p.r1a, p.r2a));
    const draw = (id, ra, rp, inc, raan, w, color) => {
      const a = (ra + rp) / 2;
      const e = (ra - rp) / (ra + rp);
      const b = a * Math.sqrt(1 - e * e);
      const [ci, si, cO, sO, cw, sw] = [inc, inc, raan, raan, w, w].map((v, k) => (k % 2 ? Math.sin(deg2rad(v)) : Math.cos(deg2rad(v))));
      // perifocal -> inertial: Rz(RAAN) Rx(i) Rz(w), first two columns
      const q11 = cO * cw - sO * sw * ci, q12 = -cO * sw - sO * cw * ci;
      const q21 = sO * cw + cO * sw * ci, q22 = -sO * sw + cO * cw * ci;
      const q31 = sw * si, q32 = cw * si;
      const X = `${a}(\\cos(t)-${e})`;
      const Y = `${b}\\sin(t)`;
      curve(calc, id, `(${q11}(${X})+${q12}(${Y}),${q21}(${X})+${q22}(${Y}),${q31}(${X})+${q32}(${Y}))`, color, '0', '2 \\pi');
    };
    draw('initial_orbit', p.r1a, p.r1p, p.i1, p.RAAN1, p.w1, BLUE);
    draw('target_orbit', p.r2a, p.r2p, p.i2, p.RAAN2, p.w2, ORANGE);
    if (result && result.arc) {
      const col = (k) => `[${result.arc.map((pt) => pt[k]).join(',')}]`;
      expr(calc, 'transfer_arc', `(${col(0)},${col(1)},${col(2)})`, { color: GREEN });
    }
  },
};

export function transferExpressions(calc, type, params, result) {
  const plot = PLOTTERS[type];
  if (plot) plot(calc, params, result);
}
