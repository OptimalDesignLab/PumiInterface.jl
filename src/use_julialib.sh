#!/bin/bash
module load gcc/4.9.2
module load pumi/core-sim
export LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/pumi/core-sim/lib:$LD_LIBRARY_PATH
