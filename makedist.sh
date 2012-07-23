#!/bin/sh

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="1.0.0"
NAME="gcg-$VERSION"
rm -f $NAME
ln -s . $NAME
if test ! -e release
then
    mkdir release
fi
rm -f release/$NAME.tgz
tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/COPYING $NAME/README $NAME/LICENSE $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk \
$NAME/check/testset/short.test $NAME/check/testset/short.solu $NAME/check/cmpres.awk \
$NAME/settings/*.set \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.h \
$NAME/check/instances/cpmp/*.gz \
$NAME/check/instances/bpp/*.gz \
$NAME/check/instances/gap/*.gz \
$NAME/check/instances/cs/*.gz \
$NAME/check/instances/miplib/*.gz \
$NAME/doc/inc/*.inc
rm -f $NAME
echo ""
echo "finished"
