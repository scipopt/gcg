#!/bin/bash
# This is an automated installer Script for a GCG installation from the RWTH Gitlab.

RED='\033[0;31m'
W='\e[0m\n'
B='\n\e[1m'

FOLDER='gcg_myrepo'
BRANCH='master'
SYS='make'
MODE='default'

function usage(){
  echo "Usage: ./installGCGFromGit.sh [options]"
  echo "Options:"
  echo "  -f, --folder    folder to install GCG to (default: 'gcg')"
  echo "  -s, --system    buildsystem (make/cmake, default: 'make')"
  echo "  -b, --branch    branch to clone (default: 'master')"
  echo "  -m, --mode      "
  echo "                  default: clone everything (recommended)"
  echo "                  fast:    only clone given branch"
  echo "                  fastest: only clone given branch and newest submodules"
  exit
}

function getArgs(){
  echo $1

  echo "OK"
}

function start(){
  while true; do
    printf "${B}Do you want to start GCG?${W}"
    read -p "Please answer. [Y/n] " yn
    case $yn in
      [Yy]* )
              if [ '$SYS' == "cmake" ]; then
                ./build/bin/gcg
              else
                ./bin/gcg
              fi
              break;;
      [Nn]* ) exit;;
          * ) printf "Please answer yes or no.";;
    esac
  done
}

function test(){
  while true; do
    printf "${B}Do you want to test GCG? (this should not take long).${W}"
    read -p "Please answer. [Y/n] " yn
    case $yn in
      [Yy]* )
              printf "  Running GCG on some sample instances...";
              make test
              printf "${B}Done!${W}";
              break;;
      [Nn]* ) break;;
          * ) printf "Please answer yes or no.";;
    esac
  done
}

function installPrerequisites(){
  while true; do
    printf "${B}Do you want to check and install the prerequisites? (not required nor possible on the chair's computers).${W}"
    read -p "Please answer. [Y/n] " yn
    case $yn in
      [Yy]* )
              echo "Installing prerequisites..."
              if [ '$SYS' == "cmake" ]; then
                sudo apt-get install cmake
              else
                sudo apt-get install make
              fi
              sudo apt-get install gcc gpp git build-essential libgmp-dev libreadline-dev zlib1g-dev bison flex libncurses-dev
              break;;
      [Nn]* ) echo "Not checking prerequisites."
              break;;
          * ) printf "Please answer yes or no.";;
    esac
  done
}

function download(){
  printf "${B}Downloading everything from the git...${W}"
  printf "${B}Cloning GCG...${W}"
  if [ $MODE = 'fast' ] || [ $MODE = 'fastest' ]; then
    git clone git@git.or.rwth-aachen.de:gcg/gcg.git $FOLDER --branch $BRANCH --depth=1
  else
    git clone git@git.or.rwth-aachen.de:gcg/gcg.git $FOLDER --branch $BRANCH
  fi
  cd $FOLDER

  # Get SCIP and SoPlex from GCG git
  printf "${B}Cloning submodules (SCIP and SoPlex)...${W}"
  git submodule init
  git submodule sync
  if [ $MODE = 'fastest' ]; then
    git submodule update --depth=1
  else
    git submodule update
  fi
}

function setLinks(){
  printf "${B}Setting links...${W}"
  # SCIP link in GCG folder
  ln -sfn ../lib/scip-git lib/scip
  # SoPlex links inside SCIP lib folder
  mkdir -p lib/scip/lib/include/spxinc
  ln -sfn $PWD/lib/soplex-git/src/* $PWD/lib/scip/lib/include/spxinc/
  mkdir -p lib/scip/lib/static/
  ln -sfn $PWD/lib/soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a $PWD/lib/scip/lib/static/libsoplex.linux.x86_64.gnu.opt.a
}

function compile(){
  # Make SCIP and SoPlex
  printf "${B}Compiling dependencies...${W}"
  make deps

  # Make GCG
  printf "${B}Compiling GCG...${W}"
  make
}

function test(){
  while true; do
    read -p "Do you want to test GCG? (this should not take long) [Y/n] " yn
    case $yn in
      [Yy]* )
              printf "Running GCG on some sample instances...";
                if [ $SYS == "cmake" ]
                then
                  make check
                elif [ $SYS == "make" ]
                then
                  make test
                fi
              printf "${B}Done!${W}";
              break;;
      [Nn]* ) exit;;
          * ) printf "Please answer yes or no.";;
    esac
  done
}

function fcmake(){
  printf "${B}Compiling GCG${W}"
  mkdir build
  cd build
  cmake .. | tail -n +91
  make
  test
}

function fmake(){
  setLinks;

  # Make SCIP and SoPlex
  printf "${B}Compiling dependencies...${W}"
  make deps

  # Make GCG
  printf "${B}Compiling GCG...${W}"sudo
  make
  test
}

function install(){
  installPrerequisites;
  download;
  if [ $SYS == 'make' ]; then
    fmake;
  else
    fcmake;
  fi
  start;
}

while true; do
  # get arguments
  while [ "$1" != "" ]; do
    case $1 in
        -b | --branch )         shift
                                BRANCH=$1
                                ;;
        -f | --folder )         shift
                                FOLDER=$1
                                ;;
        -m | --mode )           shift
                                MODE=$1
                                ;;
        -s | --system )         shift
                                SYS=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
  done
  printf "${B}GCG, SCIP and SoPlex will be installed from Git.${W}"
  echo "  Installation Mode: ${MODE}"
  if [ $MODE == 'fastest' ]; then
    echo "    WARNING: Installation Mode 'fastest' is experimental."
  fi
  echo "  Buildsystem:       ${SYS}"
  echo "  Folder:            ${FOLDER}"
  echo "  Branch:            ${BRANCH}"

  read -p "Do you wish to continue? [Y/n] " yn
  case $yn in
    [Yy]* )   install;
              break;;
    [Nn]* )   exit;;
        * )   printf "Please answer yes or no.";;
  esac
done
