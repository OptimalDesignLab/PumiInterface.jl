#!/bin/bash
#module load gcc/4.9.2
#module load pumi/core-sim
export LD_LIBRARY_PATH=`pwd`/../install/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/local/pumi/core-sim/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/../deps/core/build/install/lib:$LD_LIBRARY_PATH
export PATH=`pwd`/../deps/core/build/install/bin:$PATH
export SCOREC_PREFIX=`pwd`/../deps/core/build/install
