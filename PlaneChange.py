import numpy as np

def PlaneChange():
    orbit1params = (r1a, r1p, i1, RAAN1, w1)
    orbit2params = (r2a, r2p, i2, RAAN2, w2)
    grid = loopOverOrbits(orbit1params, orbit2params)
    coarse1, coarse2 = minCoords(grid)
    point1, point2, deltaV = NelderMead2d(coarse1, coarse2)
    plotTrajectory(point1, point2)

    def loopOverOrbits(init, final):
        numOrbitSamples = 50
        TA_array = np.linspace(0, 360, numOrbitSamples)
        minDeltaVgrid = [[0] * n for _ in range(n)]

        for TA1 in TA_array:
            r1, v1 = PFtoECIframe(orbit1params, TA1)
            for TA2 in TA_array:
                r2, v2 = PFtoECIframe(orbit2params, TA2)
                minDeltaVgrid[TA1][TA2] = minDeltaVTrajectory(r1, r2)

        return minDeltaVgrid

    def PFtoECIframe(params, TA):
        {ra, rp, i, RAAN, w} = params
        a = (ra + rp) / 2
        e = (ra - rp) / (ra + rp)
        p = a * (1 - e ** 2)
        r = p / (1 + e * Math.cos(TA))

        rpf = [[r * Math.cos(TA)], [r * Math.sin(TA)], [0]];
        vpf = Math.sqrt(mu / p) * [[-Math.sin(TA)], [e + Math.cos(TA)], [0]];
        reci = Rz(RAAN) @ Rx(i) @ Rz(w) @ rpf
        veci = Rz(RAAN) @ Rx(i) @ Rz(w) @ vpf
        
        def Rz(A):
            return np.array([[Math.cos(A), -Math.sin(A),0], [Math.sin(A), Math.cos(A),0], [0, 0, 1]])

        def Rx(A):
            return np.array([[1,0,0], [0, Math.cos(A), -Math.sin(A)], [0, Math.sin(A), Math.cos(A)]])
        
        return [reci, veci] 

    def minDeltaVTrajectory(r1, r2):
        revolutions = 1
        TOF_range = findTOFrange(r1, r2)
        deltaV = float('inf')

        for TOF in TOF_range:
            v2t, vt1 = lambertIzzoMethod(r1, r2, TOF, revolutions)
            deltaV = min(deltaV, Math.abs(v2 - v2t) + Math.abs(vt1 - v1))

        return deltaV

    def plotTrajectory(TA1, TA2):
        # Placeholder, should be in desmosSetup.js instead

    def minCoords(grid):
        minDeltaV = float('inf')
        minDVCoords = (0, 0)
        numOrbitSamples = len(grid)

        for coord1 in range(numOrbitSamples):
            for coord2 in range(numOrbitSamples):
                currentDeltaV = grid[coord1][coord2]
                if currentDeltaV < minDeltaV: 
                    minDeltaV = currentDeltaV
                    minDVCoords = (coord1, coord2)

        return minDVCoords

    def findTOFrange()
    
    def lambertIzzoMethod(r1, r2):
        # The following are the steps for the Izzo method
        # 1: Normalize geometry so that mu is 1, r's lie in R2, c = abs(r1-r2), s = (r1+r2+c)/2 and find theta
        # 2: Define Izzo param y so that the Lambert geometry can be expressed smoothly where x = y^2 / (1+y^2)
        # 3: TOF(y) = function of y, theta, r1, r2, s, c so that it's monotonic on each branch
        # 4: Find a good guess that approximiates the inverse TOF mapping y0
        # 5: Use Householder iteration of order 3 (y_(n+1) = y_n - f/f' * (1-f*f''/2f'^2) / (1-f*f''/2f'^2+f^2f'''/6f'^3))
        # 6: Once y is found, compute orbit using semimajor axis, f and g, and vt1 and v2t

    
    def Brent1d(a, b, func, tol = 1e-5, max_steps = 1000):
        fa, fb = func(a), func(b)

        if fa * fb > 0: return "Interval not bracketed properly"

        if np.abs(fa) < np.abs(fb): 
            a, b = b, a
            fa, fb = fb, fa

        c, fc = a, fa 
        d = 0
        
        mflag = True
        steps = 0

        while steps < max_steps and np.abs(a-b) > tol:
            fd = func(d)

            if (fa != fc) and (fb != fc):
                s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (b * fa * fc) / ((fb - fa) * (fb - fc)) + (c * fa * fb) / ((fc - fa) * (fc - fb))
            else: s = b - fb * (b - a) / (fb - fa)

            if (((s - (3 * a + b) / 4) * (s - b)) > 0) or (mflag and np.abs(s - b) >= np.abs(b - c) / 2) or (not mflag and np.abs(s - b) >= np.abs(c - d) / 2) or 
                (mflag and np.abs(b - c) < np.abs(tol)) or (not mflag and np.abs(c - d) < np.abs(tol)): 
                s = (a + b) / 2
                mflag = True
            else: mflag = False

            fs = func(s)
            d, c, fc = c, b, fb

            if fa * fs < 0: b, fb = s, fs
            else: a, fa = s, fs 

            if np.abs(fa) < np.abs(fb): 
                a, b = b, a
                fa, fb = fb, fa

            steps += 1

        return b, steps  
 
    def NelderMead2d(simplex, func, xtol = 1e-4, ftol = 1e-4, max_steps = 1000):
        # The function should expect a numpy array representing the points 

        alpha = 1 # alpha > 1
        gamma = 2 # gamma > 1
        rho = 0.5 # rho in [0, 0.5]
        sigma = 0.5 # sigma
        
        func_array = [(x, func(x)) for x in simplex] # an element is this array is formatted as (nparray representing point, scalar representing minDeltaV)
        steps = 0
        
        while steps < max_steps and (np.std([func_array[i][1] for i in range(3)]) > ftol or np.std([func_array[i][0] for i in range(3)]) > xtol):
            # Order
            func_array = sorted(func_array, key=lambda item: item[1])
            (x1, fx1), (x2, fx2), (x3, fx3) = func_array
            x0 = np.mean(np.array([x1, x2]), axis = 0) # Simpler to do (x1 + x2) / 2 but wanted it do be generalizable. To be fully generalizable it would be from 1 to n-1 node

            # Reflect
            xr = x0 + alpha * (x0 - x3)
            fxr = func(xr)
            if fxr >= fx1 and fxr < fx2: 
                func_array[2] = (xr, fxr)
                steps += 1
                continue

            # Expansion
            if fxr < fx1:
                xe = x0 + gamma * (xr - x0)
                fxe = func(xe)
                if fxe < fxr: 
                    func_array[2] = (xe, fxe)
                    steps += 1
                    continue
                else: 
                    func_array[2] = (xr, fxr)
                    steps += 1
                    continue

            # Contraction
            if fxr < fx3:
                xc = x0 + rho * (xr - x0)
                fxc = func(xc)
                if fxc < fxr:
                    func_array[2] = (xc, fxc)
                    steps += 1
                    continue
                else: 
                    #Shrink
                    for i, (xi, fxi) in enumerate(func_array[1:], start = 1):
                        xi = x1 + sigma * (xi - x1)
                        fxi = func(xi)
                        func_array[i] = (xi, fxi)
            else:
                xc = x0 + rho * (x3 - x0)
                fxc = func(xc)
                if fxc < fx3:
                    func_array[2] = (xc, fxc)
                    steps += 1 
                    continue
                else: 
                    #Shrink
                    for i, (xi, fxi) in enumerate(func_array[1:], start = 1):
                        xi = x1 + sigma * (xi - x1)
                        fxi = func(xi)
                        func_array[i] = (xi, fxi)
            
            steps += 1
        
        func_array = sorted(func_array, key=lambda item: item[1])
        return func_array[0] # the minimum point and the corresponding value of the function
    
    return deltaV
