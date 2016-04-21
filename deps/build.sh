#!/bin/bash

rm -vfr ./core

#git clone https://github.com/SCOREC/core.git core
if [ ! -d core ]; then
  ./download.sh
fi

cd ./core
git pull
mkdir -v ./build
cd ./build
mkdir -v ./install

# configure
../../config.sh
make -j 4
make install


# export environmental variables
ldpath="export LD_LIBRARY_PATH=`pwd`/install/lib:$LD_LIBRARY_PATH"
pkgpath="export PKG_CONFIG_PATH=`pwd`/install/lib/pkgconfig:$PKG_CONFIG_PATH"

echo $ldpath > ../../evars.sh
echo $pkgpath >> ../../evars.sh

export LD_LIBRARY_PATH=`pwd`/install/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=`pwd`/install/lib/pkgconfig:$PKG_CONFIG_PATH


