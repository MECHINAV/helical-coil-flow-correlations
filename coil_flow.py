"""Turbulent friction factor in helically coiled pipes vs. straight pipes.

Correlations (Darcy friction factor):
    Straight pipe, Blasius :  f_s = 0.3164 * Re**-0.25         (Re < ~1e5, smooth)
    Coiled pipe, Ito (1959):  f_c = [0.029 + 0.304 * (Re * (d/D)**2)**-0.25] * (d/D)**0.5
                              for 0.034 < Re*(d/D)**2 < 300
    Critical Reynolds, Ito :  Re_cr = 2e4 * (d/D)**0.32        for 15 < D/d < 860

d = tube inner diameter, D = coil (centreline) diameter.

NOTE: check these coefficients against the paper you validated in Fluent before
citing numbers. The point of this repo is a transparent baseline to compare CFD to.
"""
import numpy as np


def dean_number(re, d, big_d):
    return re * np.sqrt(d / big_d)


def friction_straight(re):
    return 0.3164 * np.asarray(re, dtype=float) ** -0.25


def friction_coiled_ito(re, d, big_d):
    re = np.asarray(re, dtype=float)
    ratio = d / big_d
    x = re * ratio**2
    if np.any((x < 0.034) | (x > 300)):
        raise ValueError("Re*(d/D)^2 outside Ito's validity range 0.034–300")
    return (0.029 + 0.304 * x**-0.25) * ratio**0.5


def critical_reynolds_ito(d, big_d):
    if not 15 <= big_d / d <= 860:
        raise ValueError("D/d outside Ito's validity range 15–860")
    return 2e4 * (d / big_d) ** 0.32


def pressure_drop(f, length, d, rho, velocity):
    """Darcy-Weisbach: dp = f * (L/d) * rho * v^2 / 2   [Pa]"""
    return f * (length / d) * 0.5 * rho * velocity**2
