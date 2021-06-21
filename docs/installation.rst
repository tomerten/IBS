.. highlight:: shell

************
Installation
************


Stable release
==============

To install ibs, run this command in your terminal:

.. code-block:: console

    $ pip install ibs

This is the preferred method to install ibs, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
============

The sources for ibs can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/tomerten/ibs

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/tomerten/ibs/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/tomerten/ibs
.. _tarball: https://github.com/tomerten/ibs/tarball/master

================================
Detailed info about Installation
================================

The generation of this package is based on `cmake_cpp_pybind11_tutorial <https://github.com/smrfeld/cmake_cpp_pybind11_tutorial>`_.

This IBS Library can be used as a independent ``C++`` Library and as ``Python`` Library.
Both are generated using ``CMake``.

The C++ Library
===============

We start by setting up the directory structure for the ``C++`` Library.


.. code-block:: console

    $ mkdir cpp 
    $ cd cpp
    $ mkdir include
    $ mkdir src 
    $ cd include 
    $ mkdir ibs_bits 


As a next step we add some code, the content of the first file will be explained later.


.. code-block:: console

    $ touch ibs 
    $ cd ibs_bits
    $ touch twiss.hpp
    $ cd ../../src 
    $ touch twiss.cpp 
    

We now add a `CMakeLists.txt` file.

.. code-block:: console

    $ cd ..
    $ touch CMakeLists.txt 

The final directory structure should look like:

.. code-block:: console

    cpp/CMakeLists.txt
    cpp/include/ibs
    cpp/include/ibs_bits/twiss.hpp
    cpp/src/twiss.cpp

The `CMakeLists.txt` file in the ``cpp`` folder will allow to build and install the independent ``C++`` Library.
Below we show the content of this file for convenience and later reference.

.. include:: ../cpp/CMakeLists.txt 
    :literal:
