#!/bin/bash

# Generate doxygen documentation for GCG.
#
# Optionally, a custom .dxy file can be passed for the doxygen configuration.

# Stop on error.
set -e

if [ "$1" == "--mathjax" ]
then
   DOXYGEN_USE_MATHJAX="YES"
   if [ -d html/MathJax ]
   then
      echo "Updating git repository for MathJax."
      cd html/MathJax
      git pull
      cd ../..
   else
      echo "Cloning git repository for MathJax."
      cd html
      git clone https://github.com/mathjax/MathJax.git
      cd ..
   fi
else
   DOXYGEN_USE_MATHJAX="NO"
fi

# Find relevant documentation versions.

CURRENT_VERSION=`grep '@version' xternal.c | awk '{ printf("%s", $3); }'`

echo "Building documentation in html/doc-${CURRENT_VERSION}."
echo "<li><a href='../doc-${CURRENT_VERSION}/index.html'>GCG ${CURRENT_VERSION}</a></li>" > docversions.html

# Create index.html and gcgheader.html.

SCIPOPTSUITEHEADER=`sed 's/\//\\\\\//g' scipoptsuiteheader.html.in | tr -d '\n'`
DOCVERSIONS=`sed 's/\//\\\\\//g' docversions.html | tr -d '\n'`

sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" -e "s/..\/doc/doc/g" < index.html.in > html/index.html
sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" < gcgheader.html.in > gcgheader.html

# Build the gcg documentation.
DOXYGEN_USE_MATHJAX=${DOXYGEN_USE_MATHJAX} doxygen gcg.dxy

echo "Cleaning up."
rm -rf html/doc-${CURRENT_VERSION} docversions.html gcgheader.html
mv html/doc html/doc-${CURRENT_VERSION}

