#!/bin/bash

module load gcc/4.9.2
#module load pumi/core-sim
module load pumi/core
#module load core-sim
export LD_LIBRARY_PATH=`pwd`/../install/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/local/pumi/core/lib:$LD_LIBRARY_PATH