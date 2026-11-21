import numpy as np
import pytest

from coil_flow import (critical_reynolds_ito, dean_number, friction_coiled_ito,
                       friction_straight, pressure_drop)


def test_blasius_value():
    assert friction_straight(1e4) == pytest.approx(0.0316, abs=1e-4)


def test_coil_exceeds_straight_pipe():
    d, big_d = 0.01, 0.2
    re = np.array([1e4, 3e4, 6e4])
    assert np.all(friction_coiled_ito(re, d, big_d) > friction_straight(re))


def test_tighter_coil_more_friction():
    re = 3e4
    assert friction_coiled_ito(re, 0.01, 0.2) > friction_coiled_ito(re, 0.01, 0.4)


def test_dean_number():
    assert dean_number(1e4, 0.01, 0.25) == pytest.approx(1e4 * np.sqrt(0.04))


def test_out_of_range_rejected():
    with pytest.raises(ValueError):
        friction_coiled_ito(10.0, 0.01, 0.2)


def test_critical_re_increases_with_curvature():
    assert critical_reynolds_ito(0.01, 0.2) > critical_reynolds_ito(0.01, 0.5)


def test_pressure_drop_scaling():
    assert pressure_drop(0.03, 2.0, 0.01, 1000, 2.0) == pytest.approx(0.03 * 200 * 0.5 * 1000 * 4)
