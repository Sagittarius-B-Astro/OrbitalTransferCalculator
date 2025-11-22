function computeCommonOrbitProperties(M, radius) {
    const G = 6.674e-11;
    const mu = G * M;

    return { mu, radius };
}

function computeHohmann(params, common) {
    const { r1, r2 } = params;
    const mu = common * 10 ** -9;
    const a = (r1 + r2) / 2;
    const ht = Math.sqrt(2 * mu * (r1 * r2) / (r1 + r2));

    const v1 = Math.sqrt(mu / r1);
    const va = ht / r1;
    const vb = ht / r2;
    const v2 = Math.sqrt(mu / r2);

    const deltaV1 = Math.abs(va - v1);
    const deltaV2 = Math.abs(v2 - vb);

    const totalDeltaV = deltaV1 + deltaV2;
    const transferTime = Math.PI * Math.sqrt(a ** 3 / mu)

    return { deltaV1, deltaV2, totalDeltaV, transferTime };
}

function computeBielliptic(params, common) {
    const { r1, r2, ri } = params;
    const mu = common * 10 ** -9;
    const a1i = (r1 + ri) / 2;
    const ai2 = (ri + r2) / 2;
    const ht1i = Math.sqrt(2 * mu * (r1 * ri) / (r1 + ri));
    const hti2 = Math.sqrt(2 * mu * (ri * r2) / (ri + r2));

    const v1 = Math.sqrt(mu / r1);
    const va1 = ht1i / r1;
    const vb1 = ht1i / ri;
    const vb2 = hti2 / ri;
    const vc2 = hti2 / r2;
    const v2 = Math.sqrt(mu / r2);

    const deltaV1 = Math.abs(va1 - v1);
    const deltaV2 = Math.abs(vb2 - vb1);
    const deltaV3 = Math.abs(v2 - vc2);

    const totalDeltaV = deltaV1 + deltaV2 + deltaV3;
    const transferTime = Math.PI * (Math.sqrt(a1i ** 3 / mu) + Math.sqrt(ai2 ** 3 / mu))

    return { deltaV1, deltaV2, deltaV3, totalDeltaV, transferTime};
}

function computeCommonApse(params, common) {
    const { r1a, r1p, r2a, r2p, A1, A2 } = params;
    const mu = common * 10 ** -9;
    const TA1 = Math.PI * A1 / 180;
    const TA2 = Math.PI * A2 / 180;
    const a1 = (r1a + r1p) / 2;
    const a2 = (r2a + r2p) / 2;
    const e1 = (r1a - r1p) / (r1a + r1p);
    const e2 = (r2a - r2p) / (r2a + r2p);
    const p1 = a1 * (1 - e1 ** 2);
    const p2 = a2 * (1 - e2 ** 2);
    const rA = p1 / (1 + e1 * Math.cos(TA1));
    const rB = p2 / (1 + e2 * Math.cos(TA2));

    const et = (rB - rA) / (rA * Math.cos(TA1) - rB * Math.cos(TA2));
    const pt = rA * rB * (Math.cos(TA1) - Math.cos(TA2)) / (rA * Math.cos(TA1) - rB * Math.cos(TA2));
    const at = pt / (1 - et ** 2);

    const h1 = Math.sqrt(mu * p1);
    const ht = Math.sqrt(mu * pt);
    const vp1 = h1 / rA;
    const vpt = ht / rA;
    const vr1 = mu / h1 * e1 * Math.sin(TA1);
    const vrt = mu / ht * et * Math.sin(TA1);
    const v1 = Math.sqrt(vp1 ** 2 + vr1 ** 2);
    const vt = Math.sqrt(vpt ** 2 + vrt ** 2);
    const phi1 = Math.atan2(vr1, vp1)
    const phit = Math.atan2(vrt, vpt)

    const totalDeltaV = Math.sqrt(v1 ** 2 + vt ** 2 - 2 * v1 * vt * Math.cos(phit - phi1));
    const gamma = Math.atan2((vrt - vr1), (vpt - vp1));
    const transferTime = Math.PI * Math.sqrt(at ** 3 / mu)

    return { totalDeltaV, gamma, transferTime};
}

function computeApseRotate(params, common) {
    const { r1a, r1p, r2a, r2p, A } = params;
    const mu = common * 10 ** -9;
    const eta = Math.PI * A / 180;
    const a1 = (r1a + r1p) / 2;
    const a2 = (r2a + r2p) / 2;
    const e1 = (r1a - r1p) / (r1a + r1p);
    const e2 = (r2a - r2p) / (r2a + r2p);
    const p1 = a1 * (1 - e1 ** 2);
    const p2 = a2 * (1 - e2 ** 2);

    const a = e1 * p2 - e2 * p1 * Math.cos(A);
    const b = -e2 * p1* Math.sin(eta);
    const c = p1 - p2;

    const alpha = Math.atan2(b, a);
    const RA1 = alpha - Math.acos(c / a * Math.cos(alpha)) % (2 * Math.PI)
    const RA2 = (RA1 - eta) % (2 * Math.PI)

    const h1 = Math.sqrt(mu * p1);
    const h2 = Math.sqrt(mu * p2);
    const r = p1 / (1 + eta * Math.cos(RA1))
    const vp1 = h1 / r;
    const vp2 = h2 / r;
    const vr1 = mu / h1 * e1 * Math.sin(RA1);
    const vr2 = mu / h2 * e2 * Math.sin(RA2);
    const v1 = Math.sqrt(vp1 ** 2 + vr1 ** 2);
    const v2 = Math.sqrt(vp2 ** 2 + vr2 ** 2);
    const phi1 = Math.atan2(vr1, vp1);
    const phi2 = Math.atan2(vr2, vp2);

    const totalDeltaV = Math.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * Math.cos(phi2 - phi1));
    const gamma = Math.atan2((vr2 - vr1), (vp2 - vp1));

    return { totalDeltaV, gamma };
}

function computerPlaneChange(params, common) {
    const { r1a, r1p, i1, RAAN1, w1, r2a, r2p, i2, RAAN2, w2 } = params;
    const mu = common * 10 ** -9;

    const orbit1params = {r1a, r1p, i1, RAAN1, w1};
    const orbit2params = {r2a, r2p, i2, RAAN2, w2};
    const grid = loopOverOrbits(orbit1params, orbit2params);
    const { coarse1, coarse2 } = min(grid); //should find min coordinates, not min deltaV
    const { point1, point2, deltaV } = corasetoFine(coarse1, coarse2);
    plotTrajectory(point1, point2);

    return deltaV;

    function loopOverOrbits(init, final) {
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
    }

    function PFtoECIframe(params, TA) {
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
    }

    function minDeltaVTrajectory(r1, r2) {
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

    function plotTrajectory(TA1, TA2) {
        // Placeholder, should be in desmosSetup.js instead
    }

    //findTOFrange()
    //lambertIzzoMethod()
    //brent1d()
    //coarsetoFine()

}