"""Lambert solutions are verified by propagating them with the independent kepler.propagate."""
import numpy as np
import pytest
from planechange.lambert_izzo import lambert, x2tof, find_xy
from planechange.kepler import propagate

MU = 3.986004418e5


def pos(r, deg):
    t = np.radians(deg)
    return r * np.array([np.cos(t), np.sin(t), 0.0])


@pytest.mark.parametrize("deg", [60, 120, 175, 240, 300])
@pytest.mark.parametrize("tof", [1000.0, 3915.0, 8000.0, 30000.0, 60000.0])
def test_every_solution_hits_target(deg, tof):
    r1, r2 = pos(9000, 0), pos(15000, deg)
    sols = lambert(r1, r2, tof, MU)
    assert len(sols) >= 1
    for s in sols:
        r, v = propagate(r1, s["v1"], tof, MU)
        assert np.linalg.norm(r - r2) < 1e-6
        assert np.linalg.norm(v - s["v2"]) < 1e-9


def test_multirev_counts():
    r1, r2 = pos(9000, 0), pos(15000, 120)
    sols = lambert(r1, r2, 30000.0, MU)
    assert sorted({s["M"] for s in sols}) == [0, 1, 2]
    assert [s["M"] for s in sols].count(1) == 2          # a left and a right branch
    assert len(lambert(r1, r2, 30000.0, MU, max_revs=1)) == 3


def test_three_dimensional_geometry():
    r1 = np.array([7000.0, 1000, 500])
    r2 = np.array([-3000.0, 9000, 4000])
    for s in lambert(r1, r2, 4000.0, MU):
        r, _ = propagate(r1, s["v1"], 4000.0, MU)
        assert np.linalg.norm(r - r2) < 1e-6


def test_retrograde_flag_reverses_angular_momentum():
    r1, r2 = pos(9000, 0), pos(15000, 120)
    pro = lambert(r1, r2, 3915.0, MU)[0]
    ret = lambert(r1, r2, 3915.0, MU, retrograde=True)[0]
    assert np.cross(r1, pro["v1"])[2] > 0 > np.cross(r1, ret["v1"])[2]


def test_tof_curve_is_continuous_across_battin_lagrange_windows():
    lam = 0.5
    xs = np.array([0.7, 0.79, 0.8, 0.801, 0.99, 0.991, 1.0, 1.01, 1.19, 1.2, 1.21, 1.5])
    t = [x2tof(x, lam, 0) for x in xs]
    assert all(np.diff(t) < 0)  # strictly decreasing in x for M = 0


def test_invalid_inputs_raise():
    with pytest.raises(AssertionError):
        find_xy(1.0, 1.0)
    with pytest.raises(AssertionError):
        find_xy(0.5, -1.0)
