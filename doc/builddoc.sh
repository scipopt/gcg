#!/bin/bash

# Generate doxygen documentation for GCG.
#
# Optionally, a custom .dxy file can be passed for the doxygen configuration.

# Stop on error.
set -e

if [[ -z ${BINDIR} ]]; then export BINDIR="$PWD/../bin"; fi

./resources/devs/howtoadd/createindexes.sh
./resources/devs/howtouse/createindexes.sh
./resources/devs/detection/classifiers/createindexes.sh
./resources/devs/detection/detectors/createindexes.sh
./resources/users/features/interactive-menu/generateMenu.sh
cd $(dirname $0)

mkdir -p html
cp -r resources/misc/scripts html

if [ "$1" == "--mathjax" ]
then
   DOXYGEN_USE_MATHJAX="YES"
   if [ -d html/MathJax ]
   then
      echo "Updating git repository for MathJax."
      cd html/MathJax
      git checkout 2.7.7
      rm *.md
      cd ../..
   else
      echo "Cloning git repository for MathJax."
      cd html
      git clone https://github.com/mathjax/MathJax.git --branch=2.7.7 --single-branch --depth 1
      rm MathJax/*.md
      cd ..
   fi
else
   DOXYGEN_USE_MATHJAX="NO"
fi

if [ "$HTML_FILE_EXTENSION" = "" ]
then
    HTML_FILE_EXTENSION=shtml
fi

# Find relevant documentation versions.

CURRENT_VERSION=`grep '@version' resources/main.md | awk '{ printf("%s", $2); }'`

echo "Building documentation in html/doc-${CURRENT_VERSION}."
echo "  Generating FAQ..."
echo "    Please ensure that you have php installed."
cd resources/misc/faq
python parser.py --linkext $HTML_FILE_EXTENSION  && php localfaq.php > faq.inc
cd ../../../
echo "<li><a href='../doc-${CURRENT_VERSION}/index.html'>GCG ${CURRENT_VERSION}</a></li>" > docversions.html

echo "  Moving visualization images..."
mkdir -p html/doc/img/visu
cp -r ./resources/devs/howtouse/visualizations/img/* ./html/doc/img/visu

echo "Generating parameters file."
echo "  Please ensure that GCG was installed correctly."
cd ..
"$BINDIR"/gcg -c "set default set save doc/resources/misc/parameters.set quit" > /dev/null 2>&1
cd doc

echo "Generating Doxygen docu..."
echo "  Please ensure that graphviz is installed on your system."
# Create index.html and gcgheader.html.
SCIPOPTSUITEHEADER=`sed 's/\//\\\\\//g' scipoptsuiteheader.html.in | tr -d '\n'`
DOCVERSIONS=`sed 's/\//\\\\\//g' docversions.html | tr -d '\n'`

sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" -e "s/..\/doc/doc/g" < index.html.in > html/index.html
sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" < gcgheader.html.in > gcgheader.html

# Build the gcg documentation.
DOXYGEN_USE_MATHJAX=${DOXYGEN_USE_MATHJAX} doxygen gcg.dxy


echo "Cleaning up..."
rm -rf html/doc-${CURRENT_VERSION} docversions.html gcgheader.html
mv html/doc html/doc-${CURRENT_VERSION}

# Remove citelist.html (the Bibliography) manually from the menu (but still reachable via link)
cd html/doc-${CURRENT_VERSION}
sed -i "/citelist/d" pages.html
sed -i "/citelist/d" navtreedata.js
sed -i "s/\:\[5/\:\[4/g" navtreeindex*.js # citelist is the third item in the navigation (after Users Guide and Devs Guide,
sed -i "s/\:\[6/\:\[5/g" navtreeindex*.js # since Installation counts as homepage and thus 0)
echo "Done."
