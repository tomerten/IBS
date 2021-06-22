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


Python wrapper
==============

For the ``Python`` wrapper we rely on  ``pybind11``. There are several
methods to this package (e.g. using conda install) but here we decided to
add it in third party folder using ``git submodules``. Therefor
to get it one needs to run:

.. code-block:: console

    $ git init 
    $ git add .
    $ git commit -am "initial commit"
    $ git submodule add https://github.com/pybind/pybind11.git third_party/pybind11-2.6.2
    $ git add .


.. note:: Note that the order of commits and immediately after adding the submodule doing 
    a new commit is important to not get errors. 

As git submodules are not updated automatically we need to add some extra code to the 
``CMakeLists.txt`` file:

.. literalinclude:: ../CMakeLists.txt
    :lines: 22-42

Wrapping the code
-----------------

In order to wrap the ``C++`` code we need to add the following
 extra code into the `Python` directory.

.. code-block:: console

    $ touch python/pytwiss.cpp

The content of this file is given by:

.. include:: ../python/pytwiss.cpp
    :literal:

creating the ``Python`` function wrappers around the ``C++`` functions.
As before we use ``CMake`` to build and install the code. Due to the fact we use 
a third party library and we need also `Python` the `CMakeLists.txt` is more complicated (for more info and a tutorial on how
to use this with ``CMake`` see `Cmake tutorial c++ <https://www.youtube.com/watch?v=ED-WUk440qc>`_.).

Using the method describe above to include ``pybind11`` one needs 
to make sure ``pybind11`` can be found by ``CMake`` which can be accomplished by adding
the following:

.. code-block:: cmake 

    SET(MY_PYBIND ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pybind11-2.5.0)
    add_subdirectory(${MY_PYBIND})

As we rely on ``Python`` their will also be extra code to set the appropriate 
directories for installing the wrapper function and making them available from
within python.

.. code-block:: cmake

    # find python
    find_program(
        PYTHON_EXECUTABLE
        NAMES python
    )
    message("Python Executable ${PYTHON_EXECUTABLE}")
    # if python found 
    # get the PYTHON_LIBRARY_DIR to allow to install files
    if(PYTHON_EXECUTABLE)
        message("Python used: ${PYTHON_EXECUTABLE}")
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c
            "from distutils import sysconfig;print(sysconfig.get_config_var('LIBDIR'))"
                OUTPUT_VARIABLE PYTHON_LIBRARY_DIR
                RESULT_VARIABLE PYTHON_LIBRARY_DIR_RESULT
                ERROR_QUIET)    
        string(STRIP ${PYTHON_LIBRARY_DIR} PYTHON_LIBRARY_DIR)
        if(NOT PYTHON_LIBRARY_DIR_RESULT MATCHES 0)
            message(SEND_ERROR "Failed to determine PYTHON_LIBRARY")
        endif()
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c "import site; print(site.getsitepackages()[0])"
            OUTPUT_VARIABLE _site_packages
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    else()
        message(FATAL_ERROR "python executable not found.")
    endif()

The full ``CMakeLists.txt`` is given below:

.. include:: ../CMakeLists.txt
    :literal:

