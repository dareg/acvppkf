#! /bin/sh

ulimit -s unlimited
export OMP_STACKARRAYS=4G
export OMP_STACKSIZE=4G

./main.x 1 89
