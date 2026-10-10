"""The p-parametrisation: closed-form velocities must match Lambert and the literature."""
import numpy as np
import pytest
from planechange.pparam import velocities_from_p, time_of_flight, p_bounds, transfer_angle
from planechange.lambert_izzo import lambert
from planechange.kepler import propagate

MU = 3.986e5  # value used by Blanco (2025)
R1, R2 = 9000.0, 15000.0


def setup(deg):
    t = np.radians(deg)
    return np.array([R1, 0, 0]), R2 * np.array([np.cos(t), np.sin(t), 0])


def test_blanco_2025_worked_example():
    r1, r2 = setup(120)
    lo, hi, par = p_bounds(R1, R2, np.radians(120))
    assert lo == pytest.approx(0.6317 * R1, rel=1e-3) and par == pytest.approx(1.8173 * R1, rel=1e-3)
    vc = np.sqrt(MU / R1)
    p = 1.3128 * R1
    v1, v2, th = velocities_from_p(r1, r2, p, MU)
    assert np.linalg.norm(v1 - [0, vc, 0]) / vc == pytest.approx(0.1563, abs=1e-4)
    tof = time_of_flight(R1, R2, th, p, v1, r1, MU)
    assert tof / (2 * np.pi * R1 / vc) == pytest.approx(0.4607, abs=1e-4)


@pytest.mark.parametrize("deg,ps", [(120, [7000, 9000, 11816, 16000, 25000]), (240, [2000, 5000, 9000, 16000])])
def test_p_family_agrees_with_izzo_lambert(deg, ps):
    r1, r2 = setup(deg)
    for p in ps:
        v1, v2, th = velocities_from_p(r1, r2, p, MU)
        tof = time_of_flight(R1, R2, th, p, v1, r1, MU)
        match = min(lambert(r1, r2, tof, MU), key=lambda s: np.linalg.norm(s["v1"] - v1))
        assert np.linalg.norm(match["v1"] - v1) < 1e-8
        assert np.linalg.norm(match["v2"] - v2) < 1e-8


def test_p_velocities_hit_the_target_by_propagation():
    r1, r2 = setup(75)
    for p in [6000, 12000, 40000]:  # elliptic, near-parabolic, hyperbolic
        v1, v2, th = velocities_from_p(r1, r2, p, MU)
        tof = time_of_flight(R1, R2, th, p, v1, r1, MU)
        r, v = propagate(r1, v1, tof, MU)
        assert np.linalg.norm(r - r2) < 1e-6


def test_extra_revolutions_change_time_but_not_velocity():
    r1, r2 = setup(120)
    v1, _, th = velocities_from_p(r1, r2, 9000.0, MU)
    t0 = time_of_flight(R1, R2, th, 9000.0, v1, r1, MU, revolutions=0)
    t2 = time_of_flight(R1, R2, th, 9000.0, v1, r1, MU, revolutions=2)
    eps = v1 @ v1 / 2 - MU / R1
    T = 2 * np.pi * np.sqrt((-MU / (2 * eps)) ** 3 / MU)
    assert t2 - t0 == pytest.approx(2 * T, rel=1e-10)


def test_tof_monotonic_in_p():
    short, long_ = setup(120), setup(240)
    for (r1, r2), sign in ((short, -1), (long_, +1)):
        ps = np.linspace(*( (5800, 40000) if sign < 0 else (500, 16000) ), 40)
        ts = []
        for p in ps:
            v1, _, th = velocities_from_p(r1, r2, p, MU)
            ts.append(time_of_flight(R1, R2, th, p, v1, r1, MU))
        assert np.all(np.sign(np.diff(ts)) == sign)  # short way decreasing, long way increasing


def test_transfer_angle_conventions():
    r1, r2 = setup(120)
    assert np.degrees(transfer_angle(r1, r2)) == pytest.approx(120)
    assert np.degrees(transfer_angle(r1, r2, retrograde=True)) == pytest.approx(240)
    assert np.degrees(transfer_angle(r1, r2, plane_normal=(0, 0, -1))) == pytest.approx(240)
