def PlaneChange():
    orbit1params = (r1a, r1p, i1, RAAN1, w1)
    orbit2params = (r2a, r2p, i2, RAAN2, w2)
    grid = loopOverOrbits(orbit1params, orbit2params)
    [coarse1, coarse2] = min(grid) //should find min coordinates, not min deltaV
    [point1, point2, deltaV] = corasetoFine(coarse1, coarse2)
    plotTrajectory(point1, point2);

    return deltaV;

    def loopOverOrbits(init, final):
        const numOrbitSamples = 50;
        const TA_array // should be an array of numOrbitSample elements equally spacing out 360
        const minDeltaVgrid = [];

        for (const TA1 of TA_array) {
            const { r1, v1 } = PFtoECIframe(orbit1params, TA1);
            for (const TA2 of TA_array) {
                const { r2, v2 } = PFtoECIframe(orbit2params, TA2);
                minDeltaVgrid[TA1][TA2] = minDeltaVTrajectory(r1, r2);
            }
        }

        return minDeltaVgrid;

    def PFtoECIframe(params, TA):
        const {ra, rp, i, RAAN, w} = params;
        const a = (ra + rp) / 2;
        const e = (ra - rp) / (ra + rp);
        const p = a * (1 - e ** 2);
        const r = p / (1 + e * Math.cos(TA));

        const rpf = [[r * Math.cos(TA)], [r * Math.sin(TA)], [0]];
        const vpf = Math.sqrt(mu / p) * [[-Math.sin(TA)], [e + Math.cos(TA)], [0]];
        const reci = Math.multiply(Rz(RAAN), Rx(i), Rz(w), rpf);
        const veci = Math.multiply(Rz(RAAN), Rx(i), Rz(w), vpf);

        function Rz(A) {
            return rot = [[Math.cos(A), -Math.sin(A),0], [Math.sin(A), Math.cos(A),0], [0, 0, 1]];
        }

        function Rx(A) {
            return rot = [[1,0,0], [0, Math.cos(A), -Math.sin(A)], [0, Math.sin(A), Math.cos(A)]];
        }
        
        return { reci, veci }

    def minDeltaVTrajectory(r1, r2) {
        var minP = float('inf');
        const revolutions = 1;
        const TOF_range = findTOFrange(r1, r2);
        const deltaV;
        for (const TOF of TOF_range) {
            const { v2t, vt1 } = lambertIzzoMethod(r1, r2, TOF, revolutions);
            deltaV = min(deltaV, Math.abs(v2 - v2t) + Math.abs(vt1 - v1));
        }

        return deltaV;
    }

    def plotTrajectory(TA1, TA2) {
        // Placeholder, should be in desmosSetup.js instead
    }

    def findTOFrange()
    
    def lambertIzzoMethod()
    
    def brent1d()
    
    def coarsetoFine()
