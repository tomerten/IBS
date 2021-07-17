#!/bin/bash

cd docs
cmake .
make clean
make Doxygen
make Sphinx