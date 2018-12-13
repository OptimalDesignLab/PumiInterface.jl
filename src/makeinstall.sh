#!/bin/bash

startdir=`pwd`
cd ./build
make -j 4
make install

cd $startdir

if [ -e ./use_julialib.sh ]; then
  rm ./use_julialib.sh
fi

# if scorec prefix not specified
if [ -z ${SCOREC_PREFIX} ]; then
  cp -vs ./use_julialib_scorec.sh ./use_julialib.sh
else
  cp -vs ./use_julialib_general.sh ./use_julialib.sh
fi

