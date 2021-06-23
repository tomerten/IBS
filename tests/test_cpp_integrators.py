#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module twiss.
"""

import os

import IBSLib as ibslib
import numpy as np
import pandas as pd
import pytest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
my_twiss_file = os.path.join(THIS_DIR, "b2_design_lattice_1996.twiss")


def test_cpp_integrator_simpson_decade():
    tau = np.empty(3)
    a = 1
    b = 2
    c = 3
    c1 = 4
    cx = 5
    cy = 6
    cprime = 7
    cyy = 8
    tl1 = 9
    tl2 = 10
    tx1 = 11
    tx2 = 12
    ty1 = 13
    ty2 = 14
    ibslib.integrator_simpson_decade(
        a, b, c, c1, cx, cy, cprime, cyy, tl1, tl2, tx1, tx2, ty1, ty2, tau
    )
    print("*** simpson-decade ***")
    print(tau)
    expected = [0.00283582, 0.00429457, 0.00605323]
    assert np.allclose(tau, expected)


def test_cpp_integrator_twsint():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    tau = np.zeros(3)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005
    sige = 7e-4
    gammas = twissheader["GAMMA"]
    len = twissheader["LENGTH"]
    bxavg = len / (2.0 * ibslib.pi * twissheader["Q1"])
    byavg = len / (2.0 * ibslib.pi * twissheader["Q2"])
    alx = 1
    aly = 2
    dx = 1
    dpx = 0.1
    dy = 0
    dpy = 0
    ibslib.integrator_twsint(
        pnumber, ex, ey, sigs, sige, gammas, bxavg, byavg, alx, aly, dx, dpx, dy, dpy, tau
    )
    print("*** simpson-twsint ***")
    print(tau)
    expected = [1.66573660e01, 1.88711448e03, 1.27508746e-02]
    assert np.allclose(tau, expected)


def test_cpp_integrand():
    l, ax, bx, a, b, c = (1, 2, 3, 4, 5, 6)
    integrand = ibslib.integrand(l, ax, bx, a, b, c)
    print("*** integrand ***")
    print(integrand)
    expected = 0.078125
    assert integrand - expected < 1e-6


def test_cpp_simpson():
    ax, bx, a, b, c, al, bl = (1, 2, 3, 4, 5, 6, 7)
    n = 10

    res = ibslib.integrator_simpson(ibslib.integrand, ax, bx, a, b, c, al, bl, n)
    print("*** simpson ***")
    print(res)

    assert res - 0.0026280561864451294 < 1e-9


def test_cpp_intSimpson():
    ax, bx, a, b, c, ay, by, ass, bs, ci = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    integral = np.zeros(3)
    ibslib.integrator_simpson_bjorken_mtingwa(
        ibslib.integrand, ax, bx, ay, by, ass, bs, a, b, ci, integral
    )
    print("*** intSimpson ***")
    print(integral)

    expected = np.array([0.72578111, 0.11911581, 0.55244817])
    assert np.allclose(integral, expected)


def test_cpp_integral_bjorken_mtingwa():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    tau = np.zeros(3)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005
    sige = 7e-4
    gammas = twissheader["GAMMA"]
    len = twissheader["LENGTH"]
    bxavg = len / (2.0 * ibslib.pi * twissheader["Q1"])
    byavg = len / (2.0 * ibslib.pi * twissheader["Q2"])
    alx = 1
    aly = 2
    dx = 1
    dpx = 0.1
    dy = 0
    dpy = 0
    ibslib.integral_bjorken_mtingwa(
        pnumber, ex, ey, sigs, sige, gammas, bxavg, byavg, alx, aly, dx, dpx, dy, dpy, tau
    )
    print("*** BM ***")
    print(tau)
    expected = [1.66576664e01, 1.88714806e03, 1.27503449e-02]
    assert np.allclose(tau, expected)


def test_cpp_integral_conte_martini():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    tau = np.zeros(3)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005
    sige = 7e-4
    gammas = twissheader["GAMMA"]
    len = twissheader["LENGTH"]
    bxavg = len / (2.0 * ibslib.pi * twissheader["Q1"])
    byavg = len / (2.0 * ibslib.pi * twissheader["Q2"])
    alx = 1
    aly = 2
    dx = 1
    dpx = 0.1
    dy = 0
    dpy = 0
    ibslib.integral_conte_martini(
        pnumber, ex, ey, sigs, sige, gammas, bxavg, byavg, alx, aly, dx, dpx, dy, dpy, tau
    )
    print("*** BM ***")
    print(tau)
    expected = [1.66573660e01, 1.88711448e03, 1.27508746e-02]
    assert np.allclose(tau, expected)


def test_cpp_integral_zimmerman():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    tau = np.zeros(3)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005
    sige = 7e-4
    gammas = twissheader["GAMMA"]
    len = twissheader["LENGTH"]
    bxavg = len / (2.0 * ibslib.pi * twissheader["Q1"])
    byavg = len / (2.0 * ibslib.pi * twissheader["Q2"])
    alx = 1
    aly = 2
    dx = 1
    dpx = 0.1
    dy = 0
    dpy = 0
    ibslib.integral_zimmerman(
        pnumber, ex, ey, sigs, sige, gammas, bxavg, byavg, alx, aly, dx, dpx, dy, dpy, tau
    )
    print("*** BM ***")
    print(tau)
    expected = [1.66573660e01, 1.88712328e03, 9.15404495e-03]
    assert np.allclose(tau, expected)
