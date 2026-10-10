import numpy as np
import pytest
from planechange.kepler import propagate
from planechange.frames import Orbit, state_at

MU = 3.986004418e5


def test_circular_quarter_period():
    r0, v0 = np.array([7000.0, 0, 0]), np.array([0, np.sqrt(MU / 7000), 0])
    T = 2 * np.pi * np.sqrt(7000 ** 3 / MU)
    r, v = propagate(r0, v0, T / 4, MU)
    assert r == pytest.approx([0, 7000, 0], abs=1e-6)
    r, v = propagate(r0, v0, T, MU)
    assert r == pytest.approx(r0, abs=1e-6)


def test_energy_and_angular_momentum_conserved_elliptic_and_hyperbolic():
    for v_scale in (1.1, 1.5):  # elliptic, hyperbolic
        r0 = np.array([7000.0, 100, 50])
        v0 = v_scale * np.sqrt(MU / 7000) * np.array([0.0, 0.8, 0.6])
        r, v = propagate(r0, v0, 4000.0, MU)
        e = lambda r, v: v @ v / 2 - MU / np.linalg.norm(r)
        assert e(r, v) == pytest.approx(e(r0, v0), rel=1e-9)
        assert np.cross(r, v) == pytest.approx(np.cross(r0, v0), rel=1e-9)


def test_orbit_properties_and_state():
    o = Orbit(14000, 7000, np.radians(30), np.radians(40), np.radians(50))
    assert o.a == 10500 and o.e == pytest.approx(1 / 3) and o.p == pytest.approx(9333.333, rel=1e-6)
    r, v = state_at(o, 0.0, MU)
    assert np.linalg.norm(r) == pytest.approx(7000)               # periapsis radius
    assert abs(r @ v) < 1e-9                                         # velocity perpendicular at periapsis
    assert np.linalg.norm(np.cross(r, v)) == pytest.approx(np.sqrt(MU * o.p))
    assert o.normal()[2] == pytest.approx(np.cos(np.radians(30)))    # h_z = cos(i)
    r, v = state_at(o, np.pi, MU)
    assert np.linalg.norm(r) == pytest.approx(14000)
