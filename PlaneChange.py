import numpy as np

def PlaneChange(r1a, r1p, i1, RAAN1, w1, r2a, r2p, i2, RAAN2, w2, mu):
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
            vt1, v2t = lambertIzzoMethod(r1, r2, TOF, revolutions, mu) # Returns two arrays for vt1 and v2t
            deltaV = min(deltaV, Math.abs(v2 - v2t) + Math.abs(vt1 - v1))

        return deltaV

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
    
    def lambertIzzoMethod(r1vec, r2vec, TOF, revolutions, mu):
        # The following are the steps for the Izzo method
        # Use Householder iteration of order 3 (y_(n+1) = y_n - f/f' * (1-f*f''/2f'^2) / (1-f*f''/2f'^2+f^2f'''/6f'^3))

        cvec = r2vec - r1vec
        c, r1, r2 = np.linalg.norm(cvec), np.linalg.norm(r1vec), np.linalg.norm(r2vec)
        s = (c + r1 + r2) / 2

        r1unit, r2unit = r1vec / r1, r2vec / r2
        hunit = np.cross(r1unit, r2unit), Lambda = np.sqrt(1 - c / s)

        if (r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0]) < 0:
            Lambda = -Lambda
            vt1unit, v2tunit = np.cross(r1unit, hunit), np.cross(r2unit, r2unit)
        else: vt1unit, v2tunit = np.cross(hunit, r1unit), np.cross(hunit, hunit)

        T = np.sqrt(2 * mu / s ** 3) * TOF
        xSols, ySols = findxy(Lambda, T)
        gamma, rho, sigma = np.sqrt(mu * s / 2), (r1 - r2) / c, np.sqrt(1 - rho ** 2) 

        vt1, v2t = [], []

        for x, y in xSols, ySols:
            Vr1, Vr2 = gamma * ((Lambda * y - x) - rho * (Lambda * y + x)) / r1, -gamma * ((Lambda * y - x) + rho * (Lambda * y + x)) / r2
            Vt1, Vt2 = gamma * sigma * (y + Lambda * x) / r1, gamma * sigma * (y + Lambda * x) / r2
            vt1vec, v2tvec = Vr1 * r1unit + Vt1 * vt1unit, Vr2 * r2unit + Vt2 * vt2unit
            vt1.append(vt1vec)
            vt2.append(v2tvec)

        return vt1, v2t

    def findxy(Lambda, T):
        assert np.abs(Lambda) < 1, "Magnitude of lambda must be less than 1"
        assert T < 0, "T must be less than 0"

        Mmax = np.floor(T / np.pi)
        T00 = np.arccos(Lambda) + Lambda * np.sqrt(1 - Lambda ** 2)

        if (T < T00 + Mmax * np.pi) and (Mmax > 0):
            Halley1d(0, ) # solve Halley iterations from x = 0, T = T0 and find Tmin(Mmax)
            if Tmin > T:
                Mmax -= 1

        T1 = 2 / 3 * (1 - Lambda ** 3)

        if T >= T0: x0 = (T0 / T) ** (2 / 3) - 1
        elif T < T1: x0 = 5 / 2 * T1 * (T1 - T) / (T * (1 - Lambda ** 5)) - 1
        else: x0 = (T0 / T) ** np.log2(T1 / T0) - 1

        x, y = Householder1d(x0, func)

        while Mmax > 0:
            x0l, x0r = (((Mmax + 1) * np.pi / (8 * T)) ** (2 / 3) - 1) / (((Mmax + 1) * np.pi / (8 * T)) ** (2 / 3) + 1), 
                ((8 * T) / (Mmax * np.pi) ** (2 / 3) - 1) / ((8 * T) / (Mmax * np.pi) ** (2 / 3) + 1)
            xr, yr = Householder1d(x0l, func)
            xl, yl = Householder1d(x0r, func)
            Mmax -= 1

        return [(x, y), (xr, yr), (xl, yl)]

    def Halley1d(x0, func, tol = 1e-5, max_steps = 10):

    def Householder1d(x0, func, tol = 1e-5, max_steps = 10):
        # Use Householder iteration of order 3 (y_(n+1) = y_n - f/f' * (1-f*f''/2f'^2) / (1-f*f''/2f'^2+f^2f'''/6f'^3))

    def GradientDescent2d():
    
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


