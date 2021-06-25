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
    ibslib.get_rad_damp_equi(
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

    res = ibslib.runODE(
        twissheader, twisstable, harmon, voltages, t, exa, eya, sigsa, sigea, 1, pnumber
    )

    print(res)
    exfinal = 5.908329e-09
    eyfinal = 2.617289e-13
    sigsfinal = 3.054559e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < 1e-6
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < 1e-6
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < 1e-6


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

    res = ibslib.runODE(
        twissheader, twisstable, harmon, voltages, t, exa, eya, sigsa, sigea, 1, pnumber, 10, 0.005
    )

    print(res)
    exfinal = 5.908329e-09
    eyfinal = 2.617127e-13
    sigsfinal = 3.054559e-03
    assert abs((res["ex"][-1] - exfinal) / exfinal) < 1e-6
    assert abs((res["ey"][-1] - eyfinal) / eyfinal) < 1e-6
    assert abs((res["sigs"][-1] - sigsfinal) / sigsfinal) < 1e-6
