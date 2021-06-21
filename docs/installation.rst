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


Detailed info about Installation
================================

The generation of this package is based on `cmake_cpp_pybind11_tutorial <https://github.com/smrfeld/cmake_cpp_pybind11_tutorial>`_.

This IBS Library can be used as a independent ``C++`` Library and as ``Python`` Library.
Both are generated using ``CMake``.

The C++ Library
---------------

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


We also need to have the main header file find all the actual libraries.

.. include:: ../cpp/include/ibs
    :literal:

Building the C++ Library
------------------------

.. code-block:: console

    $ cd cpp
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ make install


For generating an ``XCode`` project ``IBSLib.xcodeproject`` in the build director:

.. code-block:: console

    $ cd cpp/build
    $ cmake .. -GXcode

Testing the C++ Library
-----------------------

.. note:: The tests are not actual proper unit tests, but rather a place
    where some test code can be place to test the library functions.

We again start by creating the directory and adding files to them:

.. code-block:: console

    $ cd cpp 
    $ mkdir tests
    $ cd tests
    $ touch CMakeLists.txt
    $ mkdir src
    $ cd src
    $ touch test_cpp.cpp

The specific content for the a quick test of the Twiss library is shown below.

.. include:: ../cpp/tests/src/test_cpp.cpp
    :literal:


The content of the corresponding ``CMakeLists.txt`` file:


.. include:: ../cpp/tests/CMakeLists.txt
    :literal:

Now we are prepared to build and run the tests.

.. code-block:: console

    $ cd cpp/tests
    $ mkdir build
    $ cd build
    $ cmake ..
    $ cd ../bin 
    4 ./test_cpp