#!/bin/bash

# create conda python3.8 env and setting it up
conda create python=3.8 --name=ibslibenv
conda activate ibslibenv
conda install poetry

# create build dirs
mkdir -p build 
mkdir -p cpp/build
mkdir -p cpp/tests/build

# build cpp lib
cd cpp/build
cmake .. #-DCMAKE_INSTALL_PREFIX=~/.local
make
make install

# build cpp tests
cd ..
cd tests/build
cmake ..
make

# build main package
cd ../../../build
cmake ..
make
make install

# install python dep and package
cd ..
export CPLUS_INCLUDE_PATH=`pwd`/cpp/include/
poetry install
