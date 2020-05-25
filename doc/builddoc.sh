#!/bin/bash

# Generate doxygen documentation for GCG.

# Stop on error.
#set -e

# Text colors/style (bold, white, underlined, red)
W='\e[0m'
B='\e[1m'
U='\e[4m'
R='\e[91m'

# Needed for cmake compatibility
if [[ -z ${BINDIR} ]]; then export BINDIR="$PWD/../bin"; fi

# Find relevant documentation versions.
CURRENT_VERSION=`grep '@version' resources/main.md | awk '{ printf("%s", $2); }'`

# Adds new .md pages to the table of contents in the folder
# with a file named as the folder (.md).
makeSubpageIndexing () {( set -e
  DIR=$1
  TITLE=$2
  cd $DIR
  OUT=$(basename $PWD).md
  echo "# $TITLE {#$(echo $OUT | sed 's/.md//')}" > $OUT
  # Get index list and append to .md
  ls | egrep '\.md$' | sed "/$OUT/d" | sed 's/.md//' | sed 's/^/- \@subpage /' >> $OUT

  #echo "Subpage indexing for ${DIR} built sucessfully."
  cd -  > /dev/null 2>&1
)}

# Get GCG Menu for interactive menu page
makeInteractiveMenuDocu () {( set -e
  cd resources/users/features/interactive-menu
  rm -f menu.html
  
  python3 getMenu.py

  cat menu_start.html.in  > menu.html
  cat menu.txt            >> menu.html
  cat menu_end.html.in    >> menu.html

  # Remove the text file that contains all menu entries (except for submenus, e.g. master/explore)
  rm menu.txt
  cd - 
)}

# Check if mathjax is wanted and clone repository on a fixed working version
checkMathjax () {( set -e
  if [ "$1" == "--mathjax" ]
  then
     DOXYGEN_USE_MATHJAX="YES"
     if [ -d html/MathJax ]
     then
        printf ": updating repository\n"
        cd html/MathJax
        git checkout 2.7.7  > /dev/null 2>&1
        rm -f *.md
        cd ../..
     else
        printf ": cloning repository\n"
        cd html
        git clone https://github.com/mathjax/MathJax.git --branch=2.7.7 --single-branch --depth 1  > /dev/null 2>&1
        rm MathJax/*.md
        cd ..
     fi
  else
    printf ": compiling without mathjax\n"
    DOXYGEN_USE_MATHJAX="NO"
  fi
)}

# Download SCIP css files, fonts etc. such that no accesses to American sites etc. are performed
getAdditionalResources () {( set -e
  mkdir -p html/bootstrap/css
  mkdir -p html/bootstrap/js
  mkdir -p html/css
  mkdir -p html/js
  mkdir -p html/img
  # Getting Bootstrap stuff
  wget https://scip.zib.de/bootstrap/css/bootstrap.min.css --output-document html/bootstrap/css/bootstrap.min.css --no-check-certificate
  wget https://scip.zib.de/bootstrap/css/custom.css --output-document html/bootstrap/css/custom.css --no-check-certificate
  sed -i.bak 's/https:\/\/scip.zib.de\/images/..\/..\/img/g' html/bootstrap/css/custom.css && rm html/bootstrap/css/custom.css.bak
  # Getting fonts and css
  wget https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css --output-document html/bootstrap/fonts/font-awesome.min.css
  wget https://fonts.googleapis.com/css?family=Open+Sans --output-document html/bootstrap/fonts/font-googleapis.css
  wget https://fonts.gstatic.com/s/opensans/v17/mem8YaGs126MiZpBA-UFW50bbck.woff2 --output-document html/bootstrap/fonts/font-googleapis.woff2
  wget https://fonts.gstatic.com/s/opensans/v17/mem8YaGs126MiZpBA-UFW50bbck.woff2 --output-document html/bootstrap/fonts/font-googleapis.woff2
  wget https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/fonts/fontawesome-webfont.woff2 --output-document html/bootstrap/fonts/fontawesome-webfont.woff2
  # Getting js
  wget https://scip.zib.de/bootstrap/js/custom.js --output-document html/bootstrap/js/custom.js --no-check-certificate
  wget https://scip.zib.de/bootstrap/js/bootstrap.min.js --output-document html/bootstrap/js/bootstrap.min.js --no-check-certificate
  wget https://code.jquery.com/jquery.min.js --output-document html/js/jquery.min.js
  # move additional resources to html folder 
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
  echo $?
)}

# Generate interactively created FAQ (txt -> php)
generateFAQ () {( set -e
  cd resources/misc/faq
  python3 parser.py --linkext shtml  && php localfaq.php > faq.inc
  cd -
)}

# Generate parameter file (includes SCIP params)
generateParamsFile () {( set -e
  cd ..
  "$BINDIR"/gcg -c "set default set save doc/resources/misc/parameters.set quit"
  cd -
)}

# Remove citelist.html (the Bibliography) manually from the menu (but still reachable via link)
removeBibliography () {( set -e
  cd html/doc
  sed -i.bak "/citelist/d" pages.html && rm pages.html.bak
  sed -i.bak "/citelist/d" navtreedata.js && rm navtreedata.js.bak
  sed -i.bak "s/\:\[5/\:\[4/g" navtreeindex*.js && rm navtreeindex*.js.bak # citelist is the third item in the navigation (after Users Guide and Devs Guide,
  sed -i.bak "s/\:\[6/\:\[5/g" navtreeindex*.js && rm navtreeindex*.js.bak # since Installation counts as homepage and thus 0)
  cd -
)}

# Create Doxygen documentation for pages and source code
generateDoxy () {( set -e
  # add version to the dropdown
  echo "<li><a href='../doc-${CURRENT_VERSION}/index.html'>GCG ${CURRENT_VERSION}</a></li>" >> docversions.html

  # remove duplicates
  sort -u docversions.html -o docversions.html

  # Create index.html and gcgheader.html.
  SCIPOPTSUITEHEADER=`sed 's/\//\\\\\//g' scipoptsuiteheader.html.in | tr -d '\n'`
  DOCVERSIONS=`sed 's/\//\\\\\//g' docversions.html | tr -d '\n'`
  YEAR=`date +"%Y"`
  
  # Replace year by current and version by installed one
  sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" -e "s/..\/doc/doc/g" -e "s/<YEAR\/>/${YEAR}/g" -e "s/<CURRGCG\/>/${CURRENT_VERSION}/g" < index.html.in > html/index.html
  sed -e "s/<SCIPOPTSUITEHEADER\/>/${SCIPOPTSUITEHEADER}/g" -e "s/<DOCVERSIONS\/>/${DOCVERSIONS}/g" < gcgheader.html.in > gcgheader.html

  # Set mathjax flag to export it
  DOXYGEN_USE_MATHJAX="USE_MATHJAX=${DOXYGEN_USE_MATHJAX}"
  
  # Build the gcg documentation.
  printf "${R}" # make doxygen errors red
  export ${DOXYGEN_USE_MATHJAX}; doxygen gcg.dxy
  printf "${W}" 
)}

main () {
  # The $? conditions below are a try-catch method to alert the user 
  # about the origin of the issue

  n=8 # number of steps to be performed
  i=1 # step counter
  
  printf "${B}${U}Building GCG HTML Documentation in html/doc-${CURRENT_VERSION}${W}\n"
  mkdir -p html

  # Requirement: Correctly installed git
  printf "[${i}/${n}] Checking mathjax status"; let "i++"
    checkMathjax $1
    
  # Requirement: Internet connection
  echo "[${i}/${n}] Downloading additional resources"; let "i++"
    getAdditionalResources  > /dev/null 2>&1
    if [ $? -ne 0 ]; then printf " ${R}Error:${W} Please check your internet connection.\n"; fi

  # Requirement: none
  echo "[${i}/${n}] Generating subpage indexing"; let "i++"
    makeSubpageIndexing "resources/devs/howtoadd/" "How to add"
    makeSubpageIndexing "resources/devs/howtouse/" "How to use"
  
  # Requirement: Correctly installed php
  echo "[${i}/${n}] Generating FAQ"; let "i++"
    generateFAQ  > /dev/null 2>&1
    if [ $? -ne 0 ]; then printf " ${R}Error:${W} Have you installed PHP and python3 correctly?\n"; fi
    
  # Requirement: Correctly installed GCG
  echo "[${i}/${n}] Generating GCG parameters file"; let "i++"
    generateParamsFile  > /dev/null 2>&1
    if [ $? -ne 0 ]; then printf " ${R}Error:${W} Have you installed GCG correctly?\n"; fi
    
  # Requirement: Correctly installed GCG
  echo "[${i}/${n}] Generating GCG interactive menu documentation"; let "i++"
    makeInteractiveMenuDocu  > /dev/null 2>&1
    if [ $? -ne 0 ]; then printf " ${R}Error:${W} Have you installed GCG and python3 correctly?\n"; fi
  
  # Requirement: Doxygen, graphviz
  echo "[${i}/${n}] Generating Doxygen documentation"; let "i++"
    generateDoxy #> /dev/null 2>&1 # Show doxygen warnings!
    if [ $? -ne 0 ]; then printf " ${R}Error:${W} Have you installed Doxygen and graphviz correctly?\n"; fi
  
  # Requirement: none
  echo "[${i}/${n}] Finalizing"
    # remove bibliography page (used for the use case pages)
    removeBibliography  > /dev/null 2>&1
    # remove old documentation under this name
    rm -rf html/doc-${CURRENT_VERSION} gcgheader.html
    # move freshly generated docu into the desired (versionized) folder
    mv html/doc html/doc-${CURRENT_VERSION}

  printf "${B}Done!${W}\n"
}

main $@
