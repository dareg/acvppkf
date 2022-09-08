#! /bin/sh
set -x
cd build
make
ln -s main.x ../main.x
cd ..
