#!/bin/bash

rm -rf ./core
rm -r ./meshes
rm -rf ./pumi-meshes

start_dir=`pwd`
cd ../src
./cleanup.sh

cd $start_dir

exit 0
