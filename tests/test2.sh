#! /bin/sh

ulimit -s unlimited
export OMP_STACKARRAYS=4g

./main.x 1784 89
