#!/bin/bash

export LD_LIBRARY_PATH=`pwd`/core/build/install/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=`pwd`/core/build/install/lib/pkgconfig:$PKG_CONFIG_PATH
