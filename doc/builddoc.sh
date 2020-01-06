#!/bin/bash

# Generate doxygen documentation for GCG.
#
# Optionally, a custom .dxy file can be passed for the doxygen configuration.

# Stop on error.
set -e

if [[ -z ${BINDIR} ]]; then export BINDIR="$PWD/../bin"; fi

makeSubpageIndexing () {
  # Adds new .md pages to the table of contents in the folder
  # with a file named as the folder (.md).
  DIR=$1
  TITLE=$2
  cd $DIR
  OUT=$(basename $PWD).md
  echo "# $TITLE {#$(echo $OUT | sed 's/.md//')}" > $OUT
  # Get index list and append to .md
  ls | egrep '\.md$' | sed "/$OUT/d" | sed 's/.md//' | sed 's/^/- \@subpage /' >> $OUT

  #echo "Subpage indexing for ${DIR} built sucessfully."
  cd - > /dev/null 2>&1
}

makeSubpageIndexing "resources/devs/howtoadd/" "How to add"
makeSubpageIndexing "resources/devs/howtouse/" "How to use"
makeSubpageIndexing "resources/devs/detection/classifiers/clsvar/" "Variable Classifiers"
makeSubpageIndexing "resources/devs/detection/classifiers/clscons/" "Constraint Classifiers"
makeSubpageIndexing "resources/devs/detection/classifiers/clsindex/" "Index Classifiers"

./resources/users/features/interactive-menu/generateMenu.sh

mkdir -p html

if [ "$1" == "--mathjax" ]
then
   DOXYGEN_USE_MATHJAX="YES"
   if [ -d html/MathJax ]
   then
      echo "Updating git repository for MathJax."
      cd html/MathJax
      git checkout 2.7.7
      rm -f *.md
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
echo " Getting required CSS, .js and picture files..."
mkdir -p html/bootstrap/css
mkdir -p html/bootstrap/js
mkdir -p html/css
mkdir -p html/js
mkdir -p html/img

# Getting Bootstrap stuff
wget https://scip.zib.de/bootstrap/css/bootstrap.min.css --output-document html/bootstrap/css/bootstrap.min.css --quiet
wget https://scip.zib.de/bootstrap/css/custom.css --output-document html/bootstrap/css/custom.css --quiet
sed -i.bak 's/https:\/\/scip.zib.de\/images/..\/..\/img/g' html/bootstrap/css/custom.css && rm html/bootstrap/css/custom.css.bak
# Getting fonts and css
wget https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css --output-document html/bootstrap/css/font-awesome.min.css --quiet
wget https://fonts.googleapis.com/css?family=Open+Sans --output-document html/bootstrap/css/font-googleapis.css --quiet
wget https://fonts.gstatic.com/s/opensans/v17/mem8YaGs126MiZpBA-UFW50bbck.woff2 --output-document html/bootstrap/css/font-googleapis.woff2 --quiet
# Getting js
wget https://scip.zib.de/bootstrap/js/custom.js --output-document html/bootstrap/js/custom.js --quiet
wget https://scip.zib.de/bootstrap/js/bootstrap.min.js --output-document html/bootstrap/js/bootstrap.min.js --quiet
wget https://code.jquery.com/jquery.min.js --output-document html/js/jquery.min.js --quiet

echo "Generating FAQ..."
echo " Requirement: Correctly installed php."
cd resources/misc/faq
python parser.py --linkext $HTML_FILE_EXTENSION  && php localfaq.php > faq.inc
cd ../../../
echo "<li><a href='../doc-${CURRENT_VERSION}/index.html'>GCG ${CURRENT_VERSION}</a></li>" > docversions.html

echo " Moving resources to doc folder..."
cp -r resources/misc/scripts html
mkdir -p html/doc/img/visu
mkdir -p html/doc/bootstrap
mkdir -p html/doc/js
mkdir -p html/doc/css
cp -r resources/devs/howtouse/visualizations/img/* html/doc/img/visu
cp -r html/bootstrap html/doc
cp -r html/js html/doc
cp -r html/css html/doc
cp -r html/img/newscippy.png html/doc/img/
cp -r html/img/scribble_light_@2X.png html/doc/img/

echo "Generating parameters file."
echo " Requirement: Correctly installed GCG"
cd ..
"$BINDIR"/gcg -c "set default set save doc/resources/misc/parameters.set quit" > /dev/null 2>&1
cd doc

echo "Generating Doxygen docu..."
echo " Requirement: Correctly installed graphviz"
# Create index.html and gcgheader.html.
SCIPOPTSUITEHEADER=`sed 's/\//\\\\\//g' scipoptsuiteheader.html.in | tr -d '\n'`
DOCVERSIONS=`sed 's/\//\\\\\//g' docversions.html | tr -d '\n'`

sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" -e "s/..\/doc/doc/g" < index.html.in > html/index.html
sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" < gcgheader.html.in > gcgheader.html

# Build the gcg documentation.
DOXYGEN_USE_MATHJAX=${DOXYGEN_USE_MATHJAX} doxygen gcg.dxy

#pwd
#ls html/doc
#sed -i 's/https:\/\/scip.zib.de\/images//g' html/doc/bootstrap/css/custom.css

printf "\nCleaning up...\n"
rm -rf html/doc-${CURRENT_VERSION} docversions.html gcgheader.html
mv html/doc html/doc-${CURRENT_VERSION}

# Remove citelist.html (the Bibliography) manually from the menu (but still reachable via link)
cd html/doc-${CURRENT_VERSION}
sed -i.bak "/citelist/d" pages.html && rm pages.html.bak
sed -i.bak "/citelist/d" navtreedata.js && rm navtreedata.js.bak
sed -i.bak "s/\:\[5/\:\[4/g" navtreeindex*.js && rm navtreeindex*.js.bak # citelist is the third item in the navigation (after Users Guide and Devs Guide,
sed -i.bak "s/\:\[6/\:\[5/g" navtreeindex*.js && rm navtreeindex*.js.bak # since Installation counts as homepage and thus 0)
echo "Done."
