#! /bin/sh

if [ ! -d build ]; then
	mkdir build
	cd build
	cmake ..
	cd ..
fi

cd build
make -j
ln -fs $PWD/main.x $PWD/../main.x
cd ..
