import numpy as np
import pytest
from planechange.frames import Orbit, state_at
from planechange.nodal import (min_delta_v_nodal, nodal_delta_v, nodal_direction, nodal_velocities,
                               anomaly_of_direction, q_max)
from planechange.kepler import propagate
from planechange.min_dv import min_delta_v

MU = 3.986004418e5


def test_pure_plane_change_matches_closed_form():
    di = np.radians(20)
    r = min_delta_v_nodal(Orbit(7000, 7000), Orbit(7000, 7000, di), MU)
    assert r["dv"] == pytest.approx(2 * np.sqrt(MU / 7000) * np.sin(di / 2), rel=1e-9)


def test_agrees_with_general_search_on_hohmann_plus_plane_change():
    o1, o2 = Orbit(7000, 7000), Orbit(14000, 14000, np.radians(20))
    nod = min_delta_v_nodal(o1, o2, MU)
    gen = min_delta_v(o1, o2, MU, n_grid=12)
    assert nod["dv"] == pytest.approx(gen["dv"], abs=1e-5)
    assert nod["q"] == pytest.approx(0, abs=1e-6) and nod["e"] == pytest.approx(1 / 3, rel=1e-6)
    assert nod["tof"] == pytest.approx(np.pi * np.sqrt(10500 ** 3 / MU), rel=1e-6)  # apse to apse


def test_nodal_velocities_reach_target_and_match_p_and_e():
    o1 = Orbit(9000, 7000, np.radians(10), np.radians(30), np.radians(40))
    o2 = Orbit(16000, 11000, np.radians(35), np.radians(70), np.radians(10))
    res = min_delta_v_nodal(o1, o2, MU)
    r1, _ = state_at(o1, res["nu1"], MU)
    r2, _ = state_at(o2, res["nu2"], MU)
    assert r1 @ r2 / (np.linalg.norm(r1) * np.linalg.norm(r2)) == pytest.approx(-1)  # anti-parallel
    r, v = propagate(r1, res["v1"], res["tof"], MU)
    assert np.linalg.norm(r - r2) < 1e-6 and np.linalg.norm(v - res["v2"]) < 1e-9
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    assert res["p"] == pytest.approx(2 * r1n * r2n / (r1n + r2n))
    h = np.cross(r1, res["v1"])
    assert h @ h / MU == pytest.approx(res["p"])


def test_q_stationarity_condition_at_interior_optimum():
    o1 = Orbit(9000, 7000, np.radians(10), np.radians(30), np.radians(40))
    o2 = Orbit(16000, 11000, np.radians(35), np.radians(70), np.radians(10))
    res = min_delta_v_nodal(o1, o2, MU)
    r1, vc1 = state_at(o1, res["nu1"], MU)
    r2, vc2 = state_at(o2, res["nu2"], MU)
    assert abs(res["q"]) < q_max(np.linalg.norm(r1), np.linalg.norm(r2))
    u1, u2 = res["v1"] - vc1, res["v2"] - vc2
    rh = r1 / np.linalg.norm(r1)
    assert u1 @ rh / np.linalg.norm(u1) + u2 @ rh / np.linalg.norm(u2) == pytest.approx(0, abs=1e-6)


def test_delta_v_is_convex_in_q():
    o1, o2 = Orbit(9000, 7000, 0.2, 0.5, 0.7), Orbit(16000, 11000, 0.6, 1.2, 0.2)
    qs = np.linspace(-3, 0.9, 200)
    f = np.array([nodal_delta_v(o1, o2, 1, 1.0, q, MU) for q in qs])
    assert np.all(f[:-2] + f[2:] - 2 * f[1:-1] >= -1e-9)


def test_coplanar_orbits_have_no_nodal_line():
    assert nodal_direction(Orbit(7000, 7000), Orbit(9000, 8000)) is None
    with pytest.raises(ValueError):
        min_delta_v_nodal(Orbit(7000, 7000), Orbit(9000, 8000), MU)
