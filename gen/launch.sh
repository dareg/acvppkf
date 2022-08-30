#! /bin/sh
target=$1
util_target=$2
[ -f $target ] || exit 1

fxtran "$target" -o gen/out.xml
#python3 gen/simple_module.py gen/out.xml > "$util_target"
./gen/simple_module.pl gen/out.xml > test.F90 > "$util_target"
rm gen/out.xml
