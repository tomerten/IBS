cd cpp/build
cmake ..
make
sudo make install
cd ..
cd tests/build
cmake ..
make
../bin/test_integrators
cd ../../../build
cmake ..
make
make install 