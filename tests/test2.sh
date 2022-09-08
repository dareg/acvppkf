#! /bin/sh

ulimit -s unlimited
export OMP_STACKARRAYS=4g

[ -f build/main.x ] || exit 1
ln -s build/main.x main.x
./main.x 1784
