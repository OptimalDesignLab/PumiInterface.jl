#!/bin/bash
#module load gcc/4.9.2
#module load pumi/core-sim
export LD_LIBRARY_PATH=`pwd`/../install/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/local/pumi/core-sim/lib:$LD_LIBRARY_PATH

echo "SCOREX_PREFIX = $SCOREC_PREFIX"
if [[ $SCOREC_PREFIX == "" ]]  # if someone else set SCOREC_PREFIX, they are
                               # responsible for the bin and lib paths
then
  echo "setting SCOREC_PREFIX"
  export LD_LIBRARY_PATH=`pwd`/../deps/core/build/install/lib:$LD_LIBRARY_PATH
  export PATH=`pwd`/../deps/core/build/install/bin:$PATH
  export SCOREC_PREFIX=`pwd`/../deps/core/build/install
fi
