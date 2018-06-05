#!/bin/bash

if [ ! -f "core.tar.gz" ]; then 
  git clone https://github.com/SCOREC/core.git core
  tar cfvz core.tar.gz core
fi

if [ ! -f "pumi-meshes.tar.gz" ]; then
  git clone https://github.com/SCOREC/pumi-meshes.git
  tar cfvz pumi-meshes.tar.gz pumi-meshes
fi
#wget https://scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz
#tar xfz ./pumi_test_meshes.tar.gz
