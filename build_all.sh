mkdir -p cpp/build
mkdir -p cpp/build/tests
cd cpp/build
cmake ..
make
sudo make install
cd ..
cd tests/build
cmake ..
make
cd ../../../build
cmake ..
make
make install 