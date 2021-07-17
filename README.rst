===========================
IBS - Intra Beam Scattering
===========================

.. image:: https://readthedocs.org/projects/ibs/badge/
    :alt: Documentation Status
    :scale: 100%
    :target: http://ibslib.readthedocs.io

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/tomerten/ibs/main?filepath=IBS%2Fdocs%2Fnotebooks%2FTutorial.ipynb


C++ Library with PYBIND11 Python wrapper for IBS calculations, including an ODE model simulation.


*  Free software: MIT license
*  Documentation: https://ibslib.readthedocs.io.

Tutorial on Binder
------------------

Visit the link below for a live jupyter notebook with the turorial.

https://mybinder.org/v2/gh/tomerten/ibs/main?filepath=IBS%2Fdocs%2Fnotebooks%2FTutorial.ipynb

Requirements
------------

For the installation:

* CMake >= 3.10.2
* GCC compiler 
* Python >= 3.7.1
* Git - for the submodules

For the docs (see docs/requirements.txt - except doxygen):

* doxygen 
* breathe
* sphinx-click
* click
* cookiecutter
* semantic-version
* matplotlib
* ipython
* nbsphinx
* ipykernel
* pandoc

Quick start
-----------

If you do not want to go through the trouble of installing everything on your local system, 
you can also use either `Docker`_ or `Singularity`_ containers.

.. code-block:: console

    $ git clone https://github.com/tomerten/IBS.git
    $ cd IBS
    $ docker build -t ibslib:latest .
    $ sudo singularity build ibslib.sim ibslib.sdef


For quick access to the ``runode`` simulations (make sure input file and Twiss file are in the run directory):

.. code-block:: console

    $ ./ibslib.sim JSON_SIM_INPUT_FILE

For using docker have a look at the docker documentation on how to mount folders for using files 
from your local system. One interesting feature is running jupyter from within the container, giving
you direct access to the tutorial notebook.

.. code-block:: console

    $ docker run -it --user 1000 -v DIR_TO_MOUNT:/home -p 8889:8888 ibslib:latest jupyter notebook --ip=0.0.0.0 --no-browser --allow-root

After this, open a browser and go to ``localhost:8889``, copy the token from the terminal window and you
are good to go. For the turorial navigate to ``IBS/docs/notebooks/Tutorial.ipynb``. In the above case 
files from the ``DIR_TO_MOUNT`` directory will be in the home directory in the container. If you create 
a notebook and save it there it will also appear on your local (host) system.


.. _Docker: https://www.docker.com/
.. _Singularity: https://sylabs.io/guides/3.0/user-guide/quick_start.html

Source Installation
-------------------

The sources for ibs can be downloaded from the `Github repo`_.

We recommend to install everything in a dedicated `Conda` environment. 

.. code-block:: console

    $ conda create python=3.8 --name=ibslibenv
    $ conda activate ibslibenv
    $ conda install pip
    $ pip install poetry

Clone the public repository and run the `build_all.sh` script:

.. code-block:: console

    $ git clone git://github.com/tomerten/ibs
    $ cd ibs
    $ bash build_all.sh

For more details see the installation guidelines on `Installation documentation`_. 

.. _Installation documentation: https://ibslib.readthedocs.io/en/latest/installation.html
.. _Github repo: https://github.com/tomerten/ibs

Build Docs
----------

Building of the documentation is based on `this article`_ making use of Cmake, Doxygen, Breathe and Sphinx. 
It can be build locally by the following command:

.. code-block:: console

    $ bash build_docs.sh 

.. note:: 
    Make sure you have `Doxygen`_ and `Pandoc`_ installed.

To install `Doxygen` in a Conda environment:

.. code-block:: console

    $ conda install -c conda-forge doxygen

.. _Pandoc: https://pandoc.org/installing.html
.. _Doxygen: https://www.doxygen.nl/download.html
.. _this article: https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

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
-  Zimmerman (Madx - CERN note AB-2006-002) using `TWINT` and `SIMPSONDECADE` methods to perform the integration (standard Coulomb Log and with Tailcut).

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
- Equilibrium from pure radiation damping and exitation (taux, tauy, taus, exinf, eyinf, sigeoe2, sigsinf, jx, jy)
- Radiation losses per turn
- Critical omega, theta, photon energy

Numeric Functions
-----------------

- Method to calculate standard accelerator quantities.
- Various methods related to RF calculations to derive longitudinal quantities necessary for the IBS algorithms.
- Numeric functions used by specific IBS models (fmohl, rds). 

ODE 
---

The tables below show examples of the output of the ODE method, applied to the BESSY II design lattice. 
For the examples below the `stability threshold` set to `1e-3` written to a csv file.

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