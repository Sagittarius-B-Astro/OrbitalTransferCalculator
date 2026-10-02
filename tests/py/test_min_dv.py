"""End-to-end checks of the free-time minimum delta-v search against closed-form answers."""
import numpy as np
import pytest
from planechange.frames import Orbit
from planechange.min_dv import min_delta_v, delta_v
from planechange.main import solve

MU = 3.986004418e5


def hohmann_dv(r1, r2):
    at = (r1 + r2) / 2
    v = lambda r, a: np.sqrt(MU * (2 / r - 1 / a))
    return abs(v(r1, at) - v(r1, r1)) + abs(v(r2, r2) - v(r2, at))


def test_coplanar_circular_recovers_hohmann():
    res = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000), MU, n_grid=16)
    assert res["dv"] == pytest.approx(hohmann_dv(7000, 14000), rel=1e-4)
    assert res["a"] == pytest.approx(10500, rel=1e-3) and res["e"] == pytest.approx(1 / 3, rel=1e-3)


def test_pure_plane_change_is_2v_sin_half_di():
    di = np.radians(20)
    res = min_delta_v(Orbit(7000, 7000), Orbit(7000, 7000, di), MU, n_grid=16)
    assert res["dv"] == pytest.approx(2 * np.sqrt(MU / 7000) * np.sin(di / 2), rel=1e-6)


def test_combined_burn_beats_plane_change_only_at_apogee():
    di = np.radians(20)
    res = min_delta_v(Orbit(7000, 7000), Orbit(14000, 14000, di), MU, n_grid=16)
    at = 10500
    vp, va = (np.sqrt(MU * (2 / r - 1 / at)) for r in (7000, 14000))
    v1, v2 = np.sqrt(MU / 7000), np.sqrt(MU / 14000)
    naive = abs(vp - v1) + np.sqrt(va ** 2 + v2 ** 2 - 2 * va * v2 * np.cos(di))
    assert res["dv"] < naive
    assert res["dv"] > hohmann_dv(7000, 14000)  # can never beat the coplanar cost


def test_result_is_a_local_minimum():
    o1, o2 = Orbit(7000, 7000), Orbit(14000, 14000, np.radians(20))
    res = min_delta_v(o1, o2, MU, n_grid=12)
    rng = np.random.default_rng(0)
    for _ in range(25):
        d = rng.normal(scale=0.01, size=3)
        p = res["p"] * (1 + d[2])
        assert delta_v(o1, o2, res["nu1"] + d[0], res["nu2"] + d[1], p, MU) >= res["dv"] - 1e-7


def test_main_solve_json_contract():
    out = solve(dict(r1a=7000, r1p=7000, i1=0, RAAN1=0, w1=0, r2a=14000, r2p=14000, i2=20, RAAN2=0, w2=0, mu=MU, n_grid=10))
    assert {"totalDeltaV", "transferTime", "nu1Deg", "nu2Deg", "p", "a", "e", "arc"} <= out.keys()
    assert len(out["arc"]) == 81 and np.linalg.norm(out["arc"][-1]) == pytest.approx(14000, rel=1e-5)
