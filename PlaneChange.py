import numpy as np

def PlaneChange(r1a, r1p, i1, RAAN1, w1, r2a, r2p, i2, RAAN2, w2, Mmax, mu): # Requires radii in m, not km!
    TOF_range = findTOFrange(Mmax, (r1a + r1p)/2, (r2a + r2p)/2) # TOF range will be dependent on semi major axes
    # How do I choose an appropriate time range to sample and run Lambert problems on? Izzo's algorithm uses the TOF
    # to calculate the max revolutions (Mmax), and it feels questionable to use Mmax to get the range of TOFs. 

    # Upon reflection, thsi seems alright, but the following two distinctions should be made: 1) the Mmax given by the user is 
    # not the same as Mmax in findxy, and 2) I should verify that the Mmax calculated in findxy should never exceed the user input.

    orbit1params = (r1a, r1p, i1, RAAN1, w1)
    orbit2params = (r2a, r2p, i2, RAAN2, w2)

    grid = loopOverOrbits(orbit1params, orbit2params, TOF_range, Mmax, mu)
    minDVCoords = minCoords(grid, orbit1params, orbit2params, TOF_range, Mmax, mu) # Returns minimum delta V point
    bestParams = trajectoryParams(orbit1params, orbit2params, minDVCoords, mu) # Returns best possible trajectory params for plotting
    plotTrajectory(bestParams) # Placeholder

    return ans

def findTOFrange(Mmax, a1, a2): # Gets appropriate TOF range depending on Mmax input. The first element is for M = 0, second is for M > 0
    TOFs = [(0, 0)] * 2
    TOrbitSmall, TOrbitBig = 2 * np.pi * np.sqrt(a1 ** 3 / mu), 2 * np.pi * np.sqrt(a2 ** 3 / mu)
    TOFsamples = 9 # Used for coarse sampling, may give user the option later
    
    TOFs[0] = np.linspace(1e-5, TOrbitBig, TOFsamples)
    TOFs[1] = np.linspace(Mmax * TOrbitSmall, (Mmax + 1) * TOrbitBig, TOFsamples)

    return TOFs

def loopOverOrbits(init, final, TOF_range, Mmax, mu): # Creates a grid of n = numOrbitSamples uniformly sampled minDeltaV points for every (TA1, TA2) point
    numOrbitSamples = 50
    TA_array = np.linspace(0, 360, numOrbitSamples)
    minDeltaVgrid = [[0] * numOrbitSamples for _ in range(numOrbitSamples)]

    for TA1 in TA_array:
        r1, v1 = PFtoECIframe(init, TA1)
        for TA2 in TA_array:
            r2, v2 = PFtoECIframe(final, TA2)
            IzzoParams = (r1, v1, r2, v2, Mmax, mu) # Everything needed for Izzo except TOF, passed separately to distinguish and because it's an array
            minDeltaVgrid[TA1][TA2] = minDeltaV(TOF_range, IzzoParams)

    return minDeltaVgrid

def PFtoECIframe(params, TA): # Converts from PF frame to ECI frame
    ra, rp, i, RAAN, w = params
    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)
    p = a * (1 - e ** 2)
    r = p / (1 + e * np.cos(TA))

    rpf = np.array([[r * np.cos(TA)], [r * np.sin(TA)], [0]])
    vpf = np.array(np.sqrt(mu / p) * [[-np.sin(TA)], [e + np.cos(TA)], [0]])
    reci = Rz(RAAN) @ Rx(i) @ Rz(w) @ rpf
    veci = Rz(RAAN) @ Rx(i) @ Rz(w) @ vpf
    
    def Rz(A):
        return np.array([[np.cos(A), -np.sin(A),0], [np.sin(A), np.cos(A),0], [0, 0, 1]])

    def Rx(A):
        return np.array([[1,0,0], [0, np.cos(A), -np.sin(A)], [0, np.sin(A), np.cos(A)]])
    
    return [reci, veci] 

def minDeltaV(TOF_range, IzzoParams): # Determines the minimum delta V trajectory by looking at all the possible Ms and branches
    r1, v1, r2, v2, Mmax, mu = IzzoParams
    minDeltaV = float('inf')
    minDeltaVTOFleft, minDeltaVTOFright = 0, 0 # arbitrary localized range of TOF array containing minDeltaV TOF
    minDeltaVM = 0 # arbitrary M

    for M in range(0, Mmax + 1):
        # The following pair of if statements are because I created the TOF ranges to be dependent on current M; when m = 0
        # the TOF is essentially 0 to outer orbit period, and when m > 0 the TOF is from M * inner orbit period to 
        # (M + 1) * outer orbit period.  This loop currently finds the smallest delta V for a given r1, r2 and the goal is to
        # define the smallest delta V trajectory. 

        if M == 0: TOFs = TOF_range[0]
        else: TOFs = TOF_range[1]

        # To clarify, this is refining over the rough TOF range that I initialized. 
        for TOFi in range(len(TOFs)):
            # This looks basically the same as the lambertIzzoMinimizer, and that's because it is. The difference is that this one
            # takes note of which of the randomly sampled TOF points has the lowest minDeltaV. The Minimizer does not save this infomration
            # because its purpose is to return only the minimum delta V. 
            vt1, v2t = lambertIzzoMethod(r1, r2, TOFs[TOFi], M, mu) # Returns two arrays for vt1 and v2t
        
            for i in range(len(vt1)):
                currentDeltaV = np.abs(v2 - v2t[i]) + np.abs(vt1[i] - v1)

                if currentDeltaV < minDeltaV: 
                    minDeltaV = currentDeltaV
                    minDeltaVTOFleft, minDeltaVTOFright = TOFs[TOFi - 1], TOFs[TOFi + 1]

        # Figure out how to call lambertIzzo so that it's minimizable; should return minimum delta V, not two arrays
        minDeltaV = Brent1d(minDeltaVTOFleft, minDeltaVTOFright, Izzo = lambda TOF: lambertIzzoMinimizer(TOF, IzzoParams))

    return minDeltaV

def minCoords(grid, init, final, TOF_range, Mmax, mu): # Finds the location of the min delta V trajectory in the 3d plot created by loopOverOrbits
    minDeltaV = float('inf')
    minDVCoords = (0, 0)
    numOrbitSamples = len(grid)
    
    # This is simply finding the coarse minimum of the grid before refining with Nelder-Mead
    for coord1 in range(numOrbitSamples):
        for coord2 in range(numOrbitSamples):
            currentDeltaV = grid[coord1][coord2]
            if currentDeltaV < minDeltaV: 
                minDeltaV = currentDeltaV
                minDVCoords = (coord1, coord2)

    # Equilateral triangle guess of side length 2 (degrees) 
    coord1, coord2 = minDVCoords
    simplex0 = ((coord1, coord2 + 2), (coord1 + np.sqrt(3) / 2, coord2 - 1), (coord1 - np.sqrt(3) / 2, coord2 - 1)) 

    # IzzoParams must now convert from PF to ECI with TAs corresponding to coordinates of orbits instead.  
    minDVCoords, deltaV = NelderMead2d(simplex0, 
        minDeltaV = lambda TAs: minDeltaV(TOF_range, (PFtoECIframe(orbit1params, TAs[0]), PFtoECIframe(orbit1params, TAs[1], Mmax, mu)))) 

    return minDVCoords, deltaV

def trajectoryParams(orbit1params, orbit2params, TAs, mu):
    TA1, TA2 = TAs
    r1, v1 = PFtoECIframe(orbit1params, TA1)
    r2, v2 = PFtoECIframe(orbit2params, TA2)
    k = np.array([0, 0, 1])

    E1, E2 = eccentricAnomaly(r1, v1), eccentricAnomaly(r2, v2)

    rOrbitE = [x = lambda E: a * (np.cos(E) - e), y = lambda E: a * np.sqrt(1 - e ** 2) * np.sin(E), 0]

    return rOrbitE, [E1, E2]

    def eccentricAnomaly(r, v):
        hv, h = np.cross(r, v), np.linalg.norm(hv)
        ev, e = ((np.linalg.norm(v) ** 2 - mu / np.linalg.norm(r)) * r1 - np.dot(r, v) * v) / mu, np.linalg.norm(ev)
        nv, n = np.cross(k, hv), np.linalg.norm(nv)
        a = (2 / np.linalg.norm(r) - np.linalg.norm(v) ** 2 / mu) ** -1

        i = np.arccos(hv[2] / h)
        if nv[1] >= 0: RAAN = np.arccos(nv[0] / n) else RAAN = 2 * np.pi - np.arccos(nv[0]/n)
        if ev[2] >= 0: w = np.arccos(np.dot(nv, ev) / (n * e)) else w = 2 * np.pi - np.arccos(np.dot(nv, ev) / (n * e))
        if np.dot(r, v) >= 0: w = np.arccos(np.dot(ev, r) / (e * np.linalg.norm(r))) else w = 2 * np.pi - np.arccos(np.dot(ev, r) / (e * np.linalg.norm(r)))

        return 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * tan(TA / 2))

def lambertIzzoMinimizer(TOF, IzzoParams): # Finds the minimum solution to lambertIzzo method. Used for Brent and Nelder-Mead
    r1, v1, r2, v2, Mmax, mu = IzzoParams
    vt1, v2t = lambertIzzoMethod(r1, r2, TOF, M, mu)

    for i in range(len(vt1)):
        currentDeltaV = np.abs(v2 - v2t[i]) + np.abs(vt1[i] - v1)
        minDeltaV = min(currentDeltaV, minDeltaV)

    return minDeltaV

def lambertIzzoMethod(r1vec, r2vec, TOF, revolutions, mu): # Izzo's method for solving Lambert's Problem, returns velocities instead of x, y
    cvec = r2vec - r1vec
    c, r1, r2 = np.linalg.norm(cvec), np.linalg.norm(r1vec), np.linalg.norm(r2vec)
    s = (c + r1 + r2) / 2

    r1unit, r2unit = r1vec / r1, r2vec / r2
    hunit = np.cross(r1unit, r2unit), Lambda = np.sqrt(1 - c / s)

    if (r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0]) < 0:
        Lambda = -Lambda
        vt1unit, v2tunit = np.cross(r1unit, hunit), np.cross(r2unit, hunit)
    else: vt1unit, v2tunit = np.cross(hunit, r1unit), np.cross(hunit, r2unit)

    T = np.sqrt(2 * mu / s ** 3) * TOF # Seconds to dimensionless time
    xSols, ySols = findxy(Lambda, T)
    gamma, rho, sigma = np.sqrt(mu * s / 2), (r1 - r2) / c, np.sqrt(1 - rho ** 2) 

    vt1, v2t = [], []

    for x, y in xSols, ySols: # Gooding's method of converting x, y to velocities
        Vr1, Vr2 = gamma * ((Lambda * y - x) - rho * (Lambda * y + x)) / r1, -gamma * ((Lambda * y - x) + rho * (Lambda * y + x)) / r2
        Vt1, Vt2 = gamma * sigma * (y + Lambda * x) / r1, gamma * sigma * (y + Lambda * x) / r2
        vt1vec, v2tvec = Vr1 * r1unit + Vt1 * vt1unit, Vr2 * r2unit + Vt2 * vt2unit
        vt1.append(vt1vec)
        vt2.append(v2tvec)
    
    return vt1, v2t

    def findxy(Lambda, T): # Helper function for LambertIzzo solver, returns 2 * Mmax + 1 solutions
        assert np.abs(Lambda) < 1 # Magnitude of lambda must be less than 1
        assert T > 0 # T must be more than 0. There was a typo in the paper

        solutions = []

        Mmax = np.floor(T / np.pi)
        # I omitted T00 = TOF(0, 0) because I may as well calculate it in the if statement since I never use it again.

        if (T < TOF(0)) and (Mmax > 0):
            Tmin = TOF(Halley1d(0, dTOFdx, d2TOFdx2, d3TOFdx3)) # solve Halley iterations from x = 0, T = T0 and find Tmin(Mmax)
            # Note: Final steps are to calculate Tmin(Mmax) by taking the Halley's method for dT/dx = 0. This actaully makes sense since Householder
            # for dT/dx would require the fourth derivative, which is not worth calculating.  
            if Tmin > T:
                Mmax -= 1

        T1 = TOF(1)

        if T >= T0: x0 = (T0 / T) ** (2 / 3) - 1
        elif T < T1: x0 = 5 / 2 * T1 * (T1 - T) / (T * (1 - Lambda ** 5)) - 1
        else: x0 = (T0 / T) ** np.log2(T1 / T0) - 1

        x, y = Householder1d(x0, TOFroot = lambda x: TOF(x, 0) - T, dTOFdx, d2TOFdx2, d3TOFdx3), y(x)
        solutions[0] = (x, y)

        # Renamed variable because it's a little confusing on first glance and it also stores them in reverse with the exception of the M = 0 case.
        for M in range(1, Mmax + 1):
            x0l, x0r = (((M + 1) * np.pi / (8 * T)) ** (2 / 3) - 1) / (((M + 1) * np.pi / (8 * T)) ** (2 / 3) + 1), 
                ((8 * T) / (M * np.pi) ** (2 / 3) - 1) / ((8 * T) / (M * np.pi) ** (2 / 3) + 1)
            xr, yr = Householder1d(x0l, TOFroot = lambda x: TOF(x, M) - T, dTOFdx, d2TOFdx2, d3TOFdx3), y(xr)
            xl, yl = Householder1d(x0r, TOFroot = lambda x: TOF(x, M) - T, dTOFdx, d2TOFdx2, d3TOFdx3), y(xl)
            solutions[2 * M - 1] = (xr, yr)
            solutions[2 * M] = (xl, yl)

        return solutions
        
        # TOF, y equations from Izzo's paper for Lambert Method

        def y(x): return np.sqrt(1 - Lambda ** 2 * (1 - x ** 2))

        def TOF(x, M = Mmax):
            if x == 1: return 2 / 3 * (1 - Lambda ** 3)
            if x == 0: return np.arccos(Lambda) + Lambda * np.sqrt(1 - Lambda ** 2) + M * np.pi
            y = y(x)
            psi = np.arccos(x * y + Lambda * (1 - x ** 2))
            return ((psi + M * pi) / np.sqrt(abs(1 - x ** 2)) - x + Lambda * y) / (1 - x ** 2)

        def dTOFdx(x):
            if x == 1: return 2 / 5 * (Lambda ** 5 - 1)
            if x == 0: return -2
            return ((3 * TOF(x) * x) - 2 + 2 * Lambda ** 3 * x / y(x)) / (1 - x ** 2)

        def d2TOFdx2(x):
            return ((3 * TOF(x) + 5 * x * dTOFdx(x) + 2 * (1 - Lambda ** 2) * (Lambda / y(x)) ** 3)) / (1 - x ** 2)

        def d3TOFdx3(x):
            return ((7 * x * d2TOFdx2(x) + 8 * dTOFdx(x) - 6 * (1 - Lambda ** 2) * (Lambda / y(x)) ** 5 * x)) / (1 - x ** 2)

# 1d Rootfinding Algorithms, used to find minimums or roots of functions for minimum delta V trajectory call stack

# THESE BOTH CURRENTLY  ARE ONE ITERATION
def Halley1d(xn, f, dfdx, d2fdx2, tol = 1e-5, max_steps = 10):
    # Use Householder iteration of order 2 (y_(n+1)) = y_n - (2f*f') / (2f'^2 - f*f'')
    xp = xn + 1e-4
    steps = 0

    while steps < max_steps and np.abs(xn - xp) > tol:
        xn = xn - (2 * f(xn) * dfdx(xn)) / (2 * dfdx(xn) ** 2 - f(xn) * d2fdx2(xn))

    return xn

def Householder1d(xn, f, dfdx, d2fdx2, d3fdx3, tol = 1e-5, max_steps = 10):
    # Use Householder iteration of order 3 (y_(n+1) = y_n - f/f' * (1-f*f''/2f'^2) / (1-f*f''/2f'^2+f^2f'''/6f'^3))
    xp = xn + 1e-4
    steps = 0

    while steps < max_steps and np.abs(xn - xp) > tol:
        xn = xn - f(xn) * (dfdx(xn) ** 2 - f(xn) * d2fdx2(xn) / 2) / ((dfdx(xn) * (dfdx(xn) ** 2 - f(xn) * d2fdx2(xn))) + d3fdx3(xn) * f(xn) ** 2 / 6)

    return xn
    
def Brent1d(a, b, func, tol = 1e-5, max_steps = 1000): # Brent's bracketing method
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

# 2d Minimizing Algorithms, tested to find which is most efficient in finding smallest delta V

def NelderMead2d(simplex, func, xtol = 1e-4, ftol = 1e-4, max_steps = 1000): # Nelder-Mead local optimization
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