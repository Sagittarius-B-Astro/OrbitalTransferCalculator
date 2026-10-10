import numpy as np
import pytest
from planechange.optimizers import halley, householder3, brent_root, golden_min, nelder_mead

# Test functions carried over from the original test_function() in PlaneChange.py
f1 = lambda x: x ** 4 - 5 * x ** 3 + 5 * x ** 2 - 1
d1 = lambda x: 4 * x ** 3 - 15 * x ** 2 + 10 * x
d2 = lambda x: 12 * x ** 2 - 30 * x + 10
d3 = lambda x: 24 * x - 30


def test_halley_and_householder_find_same_root():
    xh, _ = halley(f1, d1, d2, 1.5)
    xr, _ = householder3(f1, d1, d2, d3, 1.5)
    assert abs(f1(xh)) < 1e-10 and abs(f1(xr)) < 1e-10
    assert xh == pytest.approx(xr, abs=1e-8)


def test_householder_converges_in_few_iterations():
    _, it = householder3(f1, d1, d2, d3, 1.5)
    assert it <= 8


def test_brent_root_polynomial_and_bracket_check():
    # NB: the original test used the bracket (0, 2), but f1(0) and f1(2) are both negative -> not a bracket
    for a, b in ((0, 0.9), (3, 4)):
        x, _ = brent_root(f1, a, b)
        assert abs(f1(x)) < 1e-9
    with pytest.raises(ValueError):
        brent_root(f1, 0, 2)
    with pytest.raises(ValueError):
        brent_root(lambda x: x ** 2 + 1, -1, 1)


def test_brent_root_sine():
    x, _ = brent_root(lambda x: np.sin(x), 3, 4)
    assert x == pytest.approx(np.pi, abs=1e-9)


def test_golden_min_parabola():
    x, fx = golden_min(lambda x: (x - 1.234) ** 2 + 5, -10, 10)
    assert x == pytest.approx(1.234, abs=1e-6) and fx == pytest.approx(5, abs=1e-10)


def test_nelder_mead_quadratic_2d_and_3d():
    x, fx, _ = nelder_mead(lambda v: (v[0] - 37) ** 2 + (v[1] - 121) ** 2, [(30, 100), (45, 100), (30, 140)])
    assert x == pytest.approx([37, 121], abs=1e-4)
    x3, _, _ = nelder_mead(lambda v: sum((v - np.array([1, -2, 3])) ** 2), [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert x3 == pytest.approx([1, -2, 3], abs=1e-4)


def test_nelder_mead_rosenbrock():
    rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
    x, fx, _ = nelder_mead(rosen, [(-1.2, 1), (-1.0, 1), (-1.2, 1.2)], max_iter=5000)
    assert x == pytest.approx([1, 1], abs=1e-3)
