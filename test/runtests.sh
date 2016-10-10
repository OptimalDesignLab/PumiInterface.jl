#!/bin/bash

jj=julia
#err=0
#$jj ./runtests.jl
#tmp=$?
#err=$((err + tmp))

mpirun -np 2 $jj ./runtests_parallel.jl
tmp=$?
err=$((err + tmp))

mpirun -np 4 $jj ./runtests_parallel4.jl
tmp=$?
err=$((err + tmp))


exit $err
