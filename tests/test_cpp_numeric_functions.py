#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module twiss.
"""

import IBSLib as ibslib
import numpy as np
import pandas as pd
import pytest

constants = [
    (ibslib.clight, 299792458.0),
    (ibslib.hbarGeV, 6.582119569e-25),
    (ibslib.electron_mass, 0.51099895000e-3),
    (ibslib.proton_mass, 0.93827208816),
    (ibslib.neutron_mass, 0.93956542052),
    (ibslib.mu_mass, 0.1056583755),
    (ibslib.atomic_mass_unit, 0.93149410242),
    (ibslib.pi, 3.141592653589793),
    (ibslib.electric_charge, 1.602176634e-19),
    (ibslib.euler, 0.577215664901533),
    (ibslib.electron_radius, 2.8179403262e-15),
    (ibslib.proton_radius, 1.5346982671888944e-18),
]


@pytest.mark.parametrize("name, value", constants)
def test_constants(name, value):
    assert name == value


def test_cpp_sigefromsigs():
    assert (ibslib.sige_from_sigs(ibslib.pi * 2 * 1.25e6, 0.005, 5e-4, 3326.0, 37.0)) < 1e-2


def test_cpp_sigsfromsige():
    val = ibslib.sigs_from_sige(8.96628617341675e-05, 3326.0, 37.0, 5e-4 * ibslib.pi * 2 * 1.25e6)
    assert (val < 0.0051) & (val > 0.004999)


def test_cpp_eta():
    assert ibslib.eta(3600.0, 37.0) - 0.00073046018996082 < 1e-9


def test_cpp_fmohl():
    a = 5.709563671168914e-04
    b = 2.329156389696222e-01
    q = 2.272866910079534e00
    npp = 1000

    actual = ibslib.fmohl(a, b, q, npp)
    expected = 6824.655537384558
    assert expected - actual < 1e-9


def test_cpp_particle_radius():
    charge = -1
    aatom = 1

    actual = ibslib.particle_radius(charge, aatom)
    expected = 1.5346982671888944e-18
    assert actual == expected


def test_cpp_BetaRelativisticFromGamma():
    gamma = 1
    expected = 0
    actual = ibslib.beta_relativistic_from_gamma(gamma)
    assert expected == actual


def test_cpp_rds():
    x, y, z = 1, 2, 3
    actual = ibslib.rds(x, y, z)
    expected = 0.29046028102188937
    assert actual == expected


hvphi = [
    ([1.0], [1.0], 90, -1.0),
    ([1.0, 17.6, 20.0], [400.0, 1200.0, 1400.0], 90, 30.742135),
]


@pytest.mark.parametrize("voltages, harmonics, phi, expected", hvphi)
def test_cpp_rfvoltages(voltages, harmonics, phi, expected):
    actual = ibslib.rf_voltage_in_ev(phi, -1.0, harmonics, voltages)
    assert actual - expected < 1.0e-6
    print(actual)


hvphip = [
    ([1.0], [1.0], 90, -6.123233995e-17),
    ([1.0, 17.6, 20.0], [400.0, 1200.0, 1400.0], 90, 30.742135),
]


@pytest.mark.parametrize("voltages, harmonics, phi, expected", hvphip)
def test_cpp_rfvoltages(voltages, harmonics, phi, expected):
    actual = ibslib.rf_voltage_in_ev_prime(phi, -1.0, harmonics, voltages)
    assert actual - expected < 1.0e-6
    print(actual)
