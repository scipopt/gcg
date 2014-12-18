#!/bin/bash

# For release versions, only use VERSION="x.x.x".
# For development versions, use VERSION="x.x.x.x" with subversion number.
VERSION="2.0.1"
NAME="gcg-$VERSION"
rm -f $NAME
ln -s . $NAME
if test ! -e release
then
    mkdir release
fi

# run git status to clean the dirty git hash
git status

echo generating default setting files
make LPS=none OPT=opt-gccold READLINE=false ZLIB=false ZIMPL=false scip -j4
make LPS=none OPT=opt-gccold READLINE=false ZLIB=false ZIMPL=false -j4
bin/gcg -c "set default set save doc/inc/parameters.set quit"

sed -i '$ d' doc/inc/parameters.set

# Before we create a tarball change the director and file rights in a command way
echo adjust file modes
find ./ -name lib -prune -o -type d -exec chmod 750 {} \;
find ./ -name lib -prune -o -type f -exec chmod 640 {} \;
find ./ -name lib -prune -o -name "*.sh" -exec chmod 750 {} \;
find ./ -name lib -prune -o -name "*.py" -exec chmod 750 {} \;
find ./ -name lib -prune -o -name "*.prl" -exec chmod 750 {} \;
find ./ -name lib -prune -o -name "hmetis" -exec chmod 750 {} \;
chmod 750 bin/*

rm -f release/$NAME.tgz
tar --no-recursion --ignore-failed-read -cvzhf release/$NAME.tgz \
--exclude="*CVS*" \
--exclude="*cvs*" \
--exclude="*~" \
--exclude=".*" \
$NAME/README $NAME/LICENSE $NAME/INSTALL $NAME/CHANGELOG $NAME/Makefile $NAME/doc/* \
$NAME/release-notes/release-notes* \
$NAME/check/check.sh $NAME/check/evalcheck.sh $NAME/check/check.awk $NAME/check/eval.sh \
$NAME/check/testset/short.test $NAME/check/testset/short.solu $NAME/check/cmpres.awk \
$NAME/settings/earlybranching.set \
$NAME/settings/heurpricing.set \
$NAME/src/depend.* \
$NAME/src/*.c $NAME/src/*.cpp $NAME/src/*.h \
$NAME/src/graph/*.cpp $NAME/src/graph/*.c $NAME/src/graph/*.h \
$NAME/check/instances/cpmp/*.lp \
$NAME/check/instances/bpp/*.lp \
$NAME/check/instances/gap/*.lp \
$NAME/check/instances/cs/*.lp \
$NAME/check/instances/miplib/*.mps \
$NAME/check/instances/cpmp/*.dec $NAME/check/instances/cpmp/*.blk \
$NAME/check/instances/bpp/*.dec $NAME/check/instances/bpp/*.blk \
$NAME/check/instances/gap/*.dec \
$NAME/check/instances/cs/*.dec \
$NAME/check/instances/miplib/*.dec \
$NAME/doc/inc/*.inc
rm -f $NAME
echo ""
echo "Building documentation and webpage"
make doc

echo ""
echo "check version numbers in main.c, doc/xternal.c, Makefile and makedist.sh ($VERSION):"
grep "VERSION" src/main.c
grep "@version" doc/xternal.c
grep "^VERSION" Makefile
tail src/githash.c


echo "Collecting webpage"
mkdir temp-webpage
pushd temp-webpage
rm -f web-$NAME
ln -s ../doc/html doc
ln -sf ../doc/index.html index.html
tar jchf ../release/web-$NAME.tbz2 doc index.html
rm -f doc
rm -f index.html
popd
rm -r temp-webpage
echo "finished"
