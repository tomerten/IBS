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

ode_threshold = 1e-4


def test_cpp_RadiationDampingApprox():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)

    out = np.zeros(7)

    ibslib.RadiationDampingApprox(
        twissheader["LENGTH"], twissheader["GAMMA"], twissheader["GAMMATR"], 4.35, 7, 10, out
    )
    print(out)

    expected = [1.54199518, 1.44441042, 0.33204837, 0.0, 0.0, 0.23375962, 0.0]

    assert np.allclose(out, expected)


def test_cpp_RadiationDampingLattice():
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)

    out = np.zeros(7)

    ibslib.RadiationDampingLattice(twisstable, out)
    print(out)

    expected = [0.26299979, 1.44292462, 0.3313656, -0.00465372, 0.0, 0.00201665, 0.0]

    assert np.allclose(out, expected)


def test_cpp_get_rad_damp_equi():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)

    out = np.zeros(7)

    ibslib.RadiationDampingLattice(twisstable, out)

    equi = np.zeros(9)
    ibslib.RadiationDampingEquilibria(
        twissheader, out, ibslib.electron_mass / ibslib.proton_mass, 6.39e-3, equi
    )
    print(equi)

    expected = [
        7.96752677e-03,
        7.99322361e-03,
        3.99017725e-03,
        5.90832712e-09,
        9.41556346e-14,
        4.87765480e-07,
        3.05467987e-03,
        1.00322520e00,
        1.00000000e00,
    ]

    assert np.allclose(equi, expected)


def test_cpp_get_energy_loss_per_turn():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)

    out = np.zeros(7)

    ibslib.RadiationDampingLattice(twisstable, out)

    loss = ibslib.get_energy_loss_per_turn(
        twissheader, out[1], ibslib.electron_mass / ibslib.proton_mass
    )

    print(loss)

    assert loss - 170261.9087271691 < 1e-6


def test_cpp_ode1():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)

    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005

    t = [0.0]
    exa = [ex]
    eya = [ey]
    sigsa = [sigs]
    sigea = []
    harmon = [400.0]
    voltages = [-4.0 * 375e3]
    coupling = 0
    threshold = 1e-3

    res = ibslib.runODE(
        twissheader,
        twisstable,
        harmon,
        voltages,
        t,
        exa,
        eya,
        sigsa,
        sigea,
        1,
        pnumber,
        coupling,
        threshold,
        "der",
    )

    print(res)
    exfinal = 6.038066e-09
    eyfinal = 2.848387e-13
    sigsfinal = 3.929674e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < ode_threshold
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < ode_threshold
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < ode_threshold


def test_cpp_ode_rlx():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)

    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005

    t = [0.0]
    exa = [ex]
    eya = [ey]
    sigsa = [sigs]
    sigea = []
    harmon = [400.0]
    voltages = [-4.0 * 375e3]
    coupling = 0
    threshold = 1e-4

    res = ibslib.runODE(
        twissheader,
        twisstable,
        harmon,
        voltages,
        t,
        exa,
        eya,
        sigsa,
        sigea,
        1,
        pnumber,
        coupling,
        threshold,
        "rlx",
    )

    # print(res)
    exfinal = 6.030113e-09
    eyfinal = 2.862700e-13
    sigsfinal = 3.915411e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < ode_threshold
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < ode_threshold
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < ode_threshold


def test_cpp_ode_der_coupling():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005

    t = [0.0]
    exa = [ex]
    eya = [ey]
    sigsa = [sigs]
    sigea = []
    harmon = [400.0]
    voltages = [-4.0 * 375e3]
    coupling = 10
    threshold = 1e-3

    res = ibslib.runODE(
        twissheader,
        twisstable,
        harmon,
        voltages,
        t,
        exa,
        eya,
        sigsa,
        sigea,
        1,
        pnumber,
        coupling,
        threshold,
        "der",
    )

    print(res)
    exfinal = 5.911103e-09
    eyfinal = 5.903462e-10
    sigsfinal = 3.085612e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < ode_threshold
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < ode_threshold
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < ode_threshold


def test_cpp_ode2():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005

    t = [0.0]
    exa = [ex]
    eya = [ey]
    sigsa = [sigs]
    sigea = []
    harmon = [400.0]
    voltages = [-4.0 * 375e3]
    coupling = 0

    res = ibslib.runODE(
        twissheader,
        twisstable,
        harmon,
        voltages,
        t,
        exa,
        eya,
        sigsa,
        sigea,
        1,
        pnumber,
        10,
        0.0005,
        coupling,
        "der",
    )

    print(res)
    exfinal = 5.679857e-09
    eyfinal = 2.635102e-11
    sigsfinal = 3.624877e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < ode_threshold
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < ode_threshold
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < ode_threshold


def test_cpp_ode2_rlx():
    twissheader = ibslib.GetTwissHeader(my_twiss_file)
    twisstable = ibslib.GetTwissTable(my_twiss_file)
    twisstable = ibslib.updateTwiss(twisstable)
    pnumber = 1e10
    ex = 5e-9
    ey = 1e-10
    sigs = 0.005

    t = [0.0]
    exa = [ex]
    eya = [ey]
    sigsa = [sigs]
    sigea = []
    harmon = [400.0]
    voltages = [-4.0 * 375e3]
    coupling = 0

    res = ibslib.runODE(
        twissheader,
        twisstable,
        harmon,
        voltages,
        t,
        exa,
        eya,
        sigsa,
        sigea,
        1,
        pnumber,
        10,
        0.0005,
        coupling,
        "rlx",
    )

    print(res)
    exfinal = 5.004572e-09
    eyfinal = 9.950159e-11
    sigsfinal = 4.990405e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < ode_threshold
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < ode_threshold
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < ode_threshold


# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (otherwise all tests are normally run with pytest)
# Make sure that you run this code with the project directory as CWD, and
# that the source directory is on the path
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_cpp_ode_rlx

    print("__main__ running", the_test_you_want_to_debug)
    the_test_you_want_to_debug()
    print("-*# finished #*-")
