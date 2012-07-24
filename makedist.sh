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
$NAME/check/instances/cpmp/*.lp \
$NAME/check/instances/bpp/*.lp \
$NAME/check/instances/gap/*.lp \
$NAME/check/instances/cs/*.lp \
$NAME/check/instances/miplib/*.mps \
$NAME/check/instances/cpmp/*.dec \
$NAME/check/instances/bpp/*.dec \
$NAME/check/instances/gap/*.dec \
$NAME/check/instances/cs/*.dec \
$NAME/check/instances/miplib/*.dec \
$NAME/doc/inc/*.inc
rm -f $NAME
echo ""
echo "finished"
