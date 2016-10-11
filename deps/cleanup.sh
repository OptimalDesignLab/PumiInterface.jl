#!/bin/bash

rm -rf ./core
rm -r ./meshes
rm -v ./pumi_test_meshes.tar.gz*

start_dir=`pwd`
cd ../src
./cleanup.sh

cd $start_dir

exit 0
