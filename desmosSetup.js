function transferExpressions(transferType) {
    if (transferType === "hohmann"){
      const r = Math.max(params['r1'], params['r2']);
      calculator.setMathBounds({
        xmin: -1.1*r,
        xmax: 1.1*r,
        ymin: -1.1*r,
        ymax: 1.1*r,
        zmin: -1.1*r,
        zmax: 1.1*r
      });
      calculator.setExpression({ id: 'r1', latex: `r_1=${params['r1']}` });
      calculator.setExpression({ id: 'r2', latex: `r_2=${params['r2']}` });
      calculator.setExpression({ id: 'a', latex: 'a=(r_1+r_2)/2' });
      calculator.setExpression({ id: 'b', latex: 'b=\\sqrt{r_1 r_2}' });
      calculator.setExpression({ id: 'e', latex: `e_c=(r_2-r_1)/(r_2+r_1)` });
      calculator.setExpression({ 
        id: 'initial_orbit', 
        expressionType: "parametric3d", 
        latex: '(r_1 \\cos(t),r_1 \\sin(t),0)', 
        color: '#00ccff', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'target_orbit', 
        expressionType: "parametric3d", 
        latex: '(r_2 \\cos(t),r_2 \\sin(t),0)', 
        color: '#ff6600', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'trajectory', 
        expressionType: "parametric3d", 
        latex: '(a(\\cos(t)-e_c),b\\sin(t),0)', 
        parametricDomain: { min: "0", max: "\\pi"} 
      });
    } else if (transferType === "bielliptic"){
      calculator.setMathBounds({
        xmin: -1.1*params['ri'],
        xmax: 1.1*params['ri'],
        ymin: -1.1*params['ri'],
        ymax: 1.1*params['ri'],
        zmin: -1.1*params['ri'],
        zmax: 1.1*params['ri']
      });
      calculator.setExpression({ id: 'r1', latex: `r_1=${params['r1']}` });
      calculator.setExpression({ id: 'r2', latex: `r_2=${params['r2']}` });
      calculator.setExpression({ id: 'ri', latex: `r_i=${params['ri']}` });
      calculator.setExpression({ id: 'a1', latex: 'a_1=(r_1+r_i)/2' });
      calculator.setExpression({ id: 'b1', latex: 'b_1=\\sqrt{r_1 r_i}' });
      calculator.setExpression({ id: 'e1', latex: `e_1=(r_i-r_1)/(r_i+r_1)` });
      calculator.setExpression({ id: 'a2', latex: 'a_2=(r_i+r_2)/2' });
      calculator.setExpression({ id: 'b2', latex: 'b_2=\\sqrt{r_i r_2}' });
      calculator.setExpression({ id: 'e2', latex: `e_2=(r_i-r_2)/(r_2+r_i)` });
      calculator.setExpression({ 
        id: 'initial_orbit', 
        expressionType: "parametric3d", 
        latex: '(r_1 \\cos(t),r_1 \\sin(t),0)',
        color: '#00ccff', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'target_orbit', 
        expressionType: "parametric3d", 
        latex: '(r_2 \\cos(t),r_2 \\sin(t),0)',
        color: '#ff6600', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'trajectory1', 
        expressionType: "parametric3d", 
        latex: '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)', 
        parametricDomain: { min: "0", max: "\\pi"} 
      });
      calculator.setExpression({ 
        id: 'trajectory2', 
        expressionType: "parametric3d", 
        latex: '(a_2(\\cos(t)-e_2),b_2\\sin(t),0)', 
        parametricDomain: { min: "\\pi", max: "2 \\pi"} 
      });
    } else if (transferType === "commonApse"){
      const r = Math.max(params['r1a'], params['r2a']);
      calculator.setMathBounds({
        xmin: -1.1*r,
        xmax: 1.1*r,
        ymin: -1.1*r,
        ymax: 1.1*r,
        zmin: -1.1*r,
        zmax: 1.1*r
      });
      calculator.setExpression({ id: 'r1a', latex: `r_{1a}=${params['r1a']}` });
      calculator.setExpression({ id: 'r1p', latex: `r_{1p}=${params['r1p']}` });
      calculator.setExpression({ id: 'r2a', latex: `r_{2a}=${params['r2a']}` });
      calculator.setExpression({ id: 'r2p', latex: `r_{2p}=${params['r2p']}` });
      calculator.setExpression({ id: 'TA1', latex: `T_1=${params['A1']}\\pi/180` });
      calculator.setExpression({ id: 'TA2', latex: `T_2=${params['A2']}\\pi/180` });
      calculator.setExpression({ id: 'a1', latex: 'a_1=(r_{1p}+r_{1a})/2' });
      calculator.setExpression({ id: 'b1', latex: 'b_1=\\sqrt{r_{1p} r_{1a}}' });
      calculator.setExpression({ id: 'e1', latex: `e_1=(r_{1a}-r_{1p})/(r_{1a}+r_{1p})` });
      calculator.setExpression({ id: 'a2', latex: 'a_2=(r_{2p}+r_{2a})/2' });
      calculator.setExpression({ id: 'b2', latex: 'b_2=\\sqrt{r_{2p} r_{2a}}' });
      calculator.setExpression({ id: 'e2', latex: `e_2=(r_{2a}-r_{2p})/(r_{2a}+r_{2p})` });
      calculator.setExpression({ id: 'p1', latex: `p_1=a_1(1-e_1^2)` });
      calculator.setExpression({ id: 'p2', latex: `p_2=a_2(1-e_2^2)` });
      calculator.setExpression({ id: 'rA', latex: `r_A=p_1/(1+e_1\\cos(T_1))` });
      calculator.setExpression({ id: 'rB', latex: `r_B=p_2/(1+e_2\\cos(T_2))` });
      calculator.setExpression({ id: 'et', latex: `e_t=(r_B-r_A)/(r_A\\cos(T_1)-r_B\\cos(T_2))` });
      calculator.setExpression({ id: 'pt', latex: `p_t=r_A r_B(\\cos(T_1)-\\cos(T_2))/(r_A\\cos(T_1)-r_B\\cos(T_2))` });
      calculator.setExpression({ id: 'at', latex: `a_t=p_t/(1-e_t^2)` });
      calculator.setExpression({ id: 'bt', latex: `b_t=a_t\\sqrt{1-e_t^2}` });
      calculator.setExpression({ 
        id: 'initial_orbit', 
        expressionType: "parametric3d", 
        latex: '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)',
        color: '#00ccff', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'target_orbit', 
        expressionType: "parametric3d", 
        latex: '(a_2(\\cos(t)-e_2),b_2\\sin(t),0)',
        color: '#ff6600', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({
        id: 'trajectory',
        expressionType: "parametric3d",
        latex: `((p_t*\\cos(t))/(1+e_t*\\cos(t)),(p_t*\\sin(t))/(1+e_t*\\cos(t)),0)`,
        color: '#00ff88',
        parametricDomain: { min: "T_2", max: "T_1" }
      });
    } else if (transferType === "apseRotate"){
      const r = Math.max(params['r1a'], params['r2a']);
      calculator.setMathBounds({
        xmin: -1.1*r,
        xmax: 1.1*r,
        ymin: -1.1*r,
        ymax: 1.1*r,
        zmin: -1.1*r,
        zmax: 1.1*r
      });
      calculator.setExpression({ id: 'r1a', latex: `r_{1a}=${params['r1a']}` });
      calculator.setExpression({ id: 'r1p', latex: `r_{1p}=${params['r1p']}` });
      calculator.setExpression({ id: 'r2a', latex: `r_{2a}=${params['r2a']}` });
      calculator.setExpression({ id: 'r2p', latex: `r_{2p}=${params['r2p']}` });
      calculator.setExpression({ id: 'eta', latex: `\\eta=${params['A']}\\pi/180` });
      calculator.setExpression({ id: 'a1', latex: 'a_1=(r_{1p}+r_{1a})/2' });
      calculator.setExpression({ id: 'b1', latex: 'b_1=\\sqrt{r_{1p} r_{1a}}' });
      calculator.setExpression({ id: 'e1', latex: `e_1=(r_{1a}-r_{1p})/(r_{1a}+r_{1p})` });
      calculator.setExpression({ id: 'a2', latex: 'a_2=(r_{2p}+r_{2a})/2' });
      calculator.setExpression({ id: 'b2', latex: 'b_2=\\sqrt{r_{2p} r_{2a}}' });
      calculator.setExpression({ id: 'e2', latex: `e_2=(r_{2a}-r_{2p})/(r_{2a}+r_{2p})` });
      calculator.setExpression({ 
        id: 'initial_orbit', 
        expressionType: "parametric3d", 
        latex: '(a_1(\\cos(t)-e_1),b_1\\sin(t),0)',
        color: '#00ccff', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
      calculator.setExpression({ 
        id: 'target_orbit', 
        expressionType: "parametric3d", 
        latex: '(((a_2(\\cos(t)-e_2))\\cos(\\eta)-(b_2\\sin(t))\\sin(\\eta)),((a_2(\\cos(t)-e_2))\\sin(\\eta)+(b_2\\sin(t))\\cos(\\eta)),0)',
        color: '#ff6600', 
        parametricDomain: { min: "0", max: "2 \\pi" }
      });
    } else if (transferType === "planeChange"){
      calculator.setMathBounds({
        xmin: -1.1*r,
        xmax: 1.1*r,
        ymin: -1.1*r,
        ymax: 1.1*r,
        zmin: -1.1*r,
        zmax: 1.1*r
      });
    }
}