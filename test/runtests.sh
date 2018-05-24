#!/bin/bash

jj=julia
jflags="--optimize=3"

err=0
$jj $jflags ./runtests.jl
tmp=$?
err=$((err + tmp))

mpirun -np 2 $jj $jflags ./runtests_parallel.jl
tmp=$?
err=$((err + tmp))

mpirun -np 4 $jj $jflags ./runtests_parallel4.jl
tmp=$?
err=$((err + tmp))


exit $err
