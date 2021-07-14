===========================
IBS - Intra Beam Scattering
===========================


C++ Library with PYBIND11 Python wrapper for IBS calculations, including an ODE model simulation.


*  Free software: MIT license
*  Documentation: https://ibs.readthedocs.io.


Features
--------

*  Read MADX Twiss module
*  Numeric Functions Lib necessary for the IBS calculations, including necessary constants.
*  Coulomb Log module
*  Radiation Damping module
*  Integration methods module (Simpson, Simpson Decade, Simpson Decade with scaling)
*  IBS models 
*  ODE simulation module with CLI interface

Current Supported Models
------------------------

-  Piwinski smooth lattice approximation
-  Piwinski Lattice element by element weighted
-  Piwinski Lattice Modified taking some vertical effects into account
-  Nagaitsev's high-energy approximation (standard Coulomb Log and with Tailcut)
-  Bjorken-Mtingwa (standard Coulomb Log and with Tailcut)
-  Conte-Martini (standard Coulomb Log and with Tailcut)
-  Madx (CERN note AB-2006-002) using `TWINT` and `SIMPSONDECADE` methods to perform the integration (standard Coulomb Log and with Tailcut).

Coublomb Log methods
--------------------

- twclog - uses element by element twiss data
- twclogtail - uses element by element twiss data
- CoublombLog - uses ring averages 
- TailCutCoulombLog - uses ring averages

Integration methods
-------------------

- Simpson (standard implementation)
- SimpsonDecade - Simspon per decade for covering large spread in integration ranges (ususally 50 orders of magnitude difference between low and high)
- TWINT, SimsponDecade with scaling method

Radiation Damping
-----------------

- Radiation Damping using smooth lattice approximation 
- Radiation Damping element by element
- Equilibrium from pure radiation damping and exitation (taux, tauy, taus, exinf, eyinf, sigeoe2, jx, jy)
- Radiation losses per turn
- Critical omega, theta, photon energy

Numeric Functions
-----------------

- sigefromsigs
- eta 
- fmohl
- particle_radius
- BetaRelativisticFromGamma
- rds 
- VeffRFeV
- VeffRFeVPrime
- synchronuousphase
- VeffRFeVPotentialWellDistortion
- VeffRFeVPotentialWellDistortionPrime
- synchronuousphasewithPWD
- synchrotronTune
- synchrotronTunewithPWD
- csige (calculates sige from RF settings, radiation losses and sigs)

ODE 
---

The tables below show examples of the output of the ODE method, with `threshold` 
set to `1e-3` written to a csv file.

.. csv-table:: ODE Model using Piwinski Smooth
    :file: cpp/tests/bin/ODE_test_output_piwinski_smooth.csv 
    :header-rows: 1

|

.. csv-table:: ODE Model using Piwinski Lattice
    :file: cpp/tests/bin/ODE_test_output_piwinski_lattice.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Piwinski Lattice Modified   
    :file: cpp/tests/bin/ODE_test_output_piwinski_latticemodified.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Nagaitsev   
    :file: cpp/tests/bin/ODE_test_output_nagaitsev.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Nagaitsev Tailcut
    :file: cpp/tests/bin/ODE_test_output_nagaitsevtailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman)
    :file: cpp/tests/bin/ODE_test_output_madx.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman) with Tailcut
    :file: cpp/tests/bin/ODE_test_output_madxtailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with standard Simpson integration (Fails for ey)
    :file: cpp/tests/bin/ODE_test_output_bjorken_mtingwa2.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with Simpson Decade Integration 
    :file: cpp/tests/bin/ODE_test_output_bjorken_mtingwa.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Bjorken-Mtingwa with Simpson Decade Integration and Tailcut
    :file: cpp/tests/bin/ODE_test_output_bjorken_mtingwatailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Conte-Martini using Simspon Decade Integration
    :file: cpp/tests/bin/ODE_test_output_conte_martini.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using Conte-Martini using Simspon Decade Integration and Tailcut
    :file: cpp/tests/bin/ODE_test_output_conte_martini_tailcut.csv 
    :header-rows: 1

|
    
.. csv-table::  ODE Model using MADX (Zimmerman) using Simpson Decade Integration 
    :file: cpp/tests/bin/ODE_test_output_madxibs.csv 
    :header-rows: 1