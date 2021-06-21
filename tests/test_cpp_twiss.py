#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module twiss.
"""

import IBSLib as ibslib
import pandas as pd


def test_GetTwissHeader():
    expected = {
        "ALFA": 0.000731694275,
        "BCURRENT": 0.002001335116,
        "BETXMAX": 17.43722834,
        "BETYMAX": 21.09492526,
        "BV_FLAG": 1.0,
        "CHARGE": -1.0,
        "DELTAP": 0.0,
        "DQ1": 1.07774528,
        "DQ2": 1.371772343,
        "DXMAX": 0.4496533373,
        "DXRMS": 0.2243964618,
        "DYMAX": 0.0,
        "DYRMS": 0.0,
        "ENERGY": 1.7,
        "ET": 0.001,
        "EX": 5e-09,
        "EY": 5e-12,
        "GAMMA": 3326.817012,
        "GAMMATR": 36.96878449,
        "KBUNCH": 1.0,
        "LENGTH": 240.0,
        "MASS": 0.00051099895,
        "NPART": 10000000000.0,
        "ORBIT5": -0.0,
        "PC": 1.699999923,
        "Q1": 17.84996564,
        "Q2": 6.743031467,
        "SIGE": 0.0001,
        "SIGT": 0.005,
        "SYNCH_1": 0.175606626,
        "SYNCH_2": 1.44292462,
        "SYNCH_3": 0.3313655984,
        "SYNCH_4": -0.00465371553,
        "SYNCH_5": 0.001749753835,
        "XCOMAX": 0.0,
        "XCORMS": 0.0,
        "YCOMAX": 0.0,
        "YCORMS": 0.0,
    }
    actual = ibslib.GetTwissHeader("b2_design_lattice_1996.twiss")

    assert expected == actual


if __name__ == "__main__":
    twissheader = ibslib.GetTwissHeader("b2_design_lattice_1996.twiss")
    print(twissheader)
