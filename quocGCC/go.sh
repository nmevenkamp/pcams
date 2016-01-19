#!/bin/bash
export CC="/Applications/MacPorts/bin/gcc-mp-5"
export CXX="/Applications/MacPorts/bin/g++-mp-5"
cmake -G"Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=1 -DDYNAMIC_LINKING=1 -DUSE_FFTW=1 -DUSE_PNG=1 -DUSE_BLAS=1 -DUSE_LAPACK=1 ../quocmesh