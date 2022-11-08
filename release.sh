#! /bin/sh
if [ -d acvppkf ];then
	echo "acvppkf directory already exists"
	exit 1
fi
if [ -f acpluiz.tar.gz ];then
	echo "acvppkf.tar.gz already exists"
	exit 1
fi

mkdir acvppkf
mkdir acvppkf/DATA
cp DATA/ACVPPKF.IN.DAT acvppkf/DATA/
cp DATA/ACVPPKF.CONST.DAT acvppkf/DATA/
cp DATA/ACVPPKF.OUT.DAT acvppkf/DATA/
cp -r data_medium acvppkf/
cp -r tests acvppkf/
cp -r src acvppkf/
cp CMakeLists.txt acvppkf/
cp main.F90 acvppkf/
tar -czvf acvppkf.tar.gz acvppkf
rm -r acvppkf
