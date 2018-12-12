#!/bin/bash

# sourcing this file sets the environment variables necessary to link to
# the Pumi that might have been built locally

if [[ $SCOREC_PREFIX == "" ]]  # if someone else set SCOREC_PREFIX, they are
                               # responsible for the bin and lib paths
then
  export LD_LIBRARY_PATH=`pwd`/../deps/core/build/install/lib:$LD_LIBRARY_PATH
  export PATH=`pwd`/../deps/core/build/install/bin:$PATH
  export SCOREC_PREFIX=`pwd`/../deps/core/build/install
fi
