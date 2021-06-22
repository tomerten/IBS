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


def test_cpp_twclog():
    clog = np.empty(2)
    ibslib.twclog(1e10, 7, 12, 1, 0, 5e-9, 1e-10, 1e-18, 3600, 1, 1700, 0.5, 0.005, 2e-9, clog)
    print(clog)
    expected = np.array([2.29372683e01, 3.25793624e-03])
    # are not exaclty due to rounding in printing out -> allclose
    assert np.allclose(clog, expected)


def test_cpp_twclogtailcut():
    clog = np.empty(2)
    ibslib.twclogtail(
        1e10,
        0.5,
        7,
        12,
        1,
        -0.1,
        0,
        0,
        0.1,
        0,
        0.0,
        0.0,
        0.0,
        5e-9,
        1e-10,
        1e-18,
        ibslib.electron_mass / ibslib.proton_mass,
        3600,
        1700,
        240,
        ibslib.electron_mass,
        -1.0,
        7e-4,
        0.005,
        clog,
    )
    print(clog)
    expected = np.array([2.31189787e01, 9.38213087e-09])
    # are not exaclty due to rounding in printing out -> allclose
    assert np.allclose(clog, expected)
