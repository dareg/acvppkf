#! /bin/sh
if [ -d acvppkf ];then
	echo "acvppkf directory already exists"
	exit 1
fi
if [ -f acpluiz.tar.xz ];then
	echo "acvppkf.tar.gz already exists"
	exit 1
fi

mkdir acvppkf
cp -r DATA acvppkf/
cp CMakeLists.txt acvppkf/
cp -r src acvppkf/
cp main.F90 acvppkf/
tar -cJvf acvppkf.tar.xz acvppkf
rm -r acvppkf
