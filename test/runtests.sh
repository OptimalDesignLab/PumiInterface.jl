#!/bin/bash

jj=julia
err=0
$jj ./runtests.jl
tmp=$?
err=$((err + tmp))
echo "\nafter first tests, err = $err"
mpirun -np 2 $jj ./runtests_parallel.jl
tmp=$?
err=$((err + tmp))

echo "\nafter second tests, err = $err"

exit $err
