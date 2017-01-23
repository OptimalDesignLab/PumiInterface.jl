#!/bin/bash

mkdir -v ./install

cmake .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
  -DENABLE_THREADS=OFF \
  -DBUILD_SHARED_LIBS=True \
  -DCMAKE_INSTALL_PREFIX=`pwd`"/install" \
  -DIS_TESTING=True \
  -DMESHES=`pwd`"/../../meshes"
#  -DENABLE_ZOLTAN=True
