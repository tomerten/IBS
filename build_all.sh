mkdir -p cpp/build
mkdir -p cpp/build/tests
cd cpp/build
cmake .. -DCMAKE_INSTALL_PREFIX=~/.local
make
make install
cd ..
cd tests/build
cmake ..
make
cd ../../../build
cmake ..
make
make install 
cd ..
export CPLUS_INCLUDE_PATH=`pwd`/cpp/include/
poetry install