#!/bin/bash

name=`basename $1 .test`
settings=$2
testset=testset/{name}.test
newtest=testset/${name}dec.test

pushd check/decs
rm -f *.dec
mkdir -p $name.$settings/gp
popd

make TEST=$name MODE=detectall SETTINGS=$settings test

pushd check
pushd decs
mv *.dec $name.$settings/
popd
mv *.gp decs/$name.$settings/gp/

echo > $newtest
for file in `cat testset/${name}.test`;
do
  f=`basename ${file} .gz`;
  f=${f%.mps}
  f=${f%.lp}
  for dec in `ls -1 decs/$name.$settings/ |grep $f`
  do
    echo "$file;decs/$name.$settings/$dec" >> $newtest
  done
done
popd
