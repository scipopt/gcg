#!/bin/bash
# This is an automated installer Script for the SCIP Optimization Suite including GCG.

RED='\033[0;31m'
W='\e[0m\n'
B='\n\e[1m'

function test(){
  while true; do
    read -p "Do you want to test the suite? (this should not take long) [Y/n] " yn
    case $yn in
      [Yy]* )
              printf "Running GCG on some sample instances...";
              cd gcg/
                if [ $1 = 1 ]
                then
                  make gcg_check
                elif [ $1 = 2 ]
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
  printf "${B}Updating CMake${W}"
  sudo apt-get install cmake
  printf "${B}Compiling the SCIP Optimization Suite${W}"
  mkdir build
  cd build
  cmake .. | tail -n +91
  make
  printf "${B}Creating Executable${W}"
  make install
  test "1"
}

function fmake(){
  printf "${B}Updating Make${W}"
  sudo apt-get install make
  printf "${B}Compiling SCIP, SoPlex and ZIMPL${W}"
  make scipoptlib
  printf "${B}Compiling GCG${W}"
  make gcg
  printf "${B}Creating Executable${W}"
  make install
  test "2"
}

function installPrerequisites(){
  while true; do
    printf "${B}Do you want to check and install the prerequisites?${W}"
    read -p "Please answer. [Y/n] " yn
    case $yn in
      [Yy]* )
              echo "Updating Package Manager..." 
              sudo apt-get update
              echo "Installing prerequisites..."
              sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev lib32ncurses5-dev libboost-program-options-dev
              break;;
      [Nn]* ) echo "Not checking prerequisites."
              break;;
          * ) printf "Please answer yes or no.";;
    esac
  done
}

function install(){
  installPrerequisites;
  printf "${B}Initializing Installation...${W}"
  while [[ ! -f $SCIPtar || -z $SCIPtar ]]; do
    read -e -p " Please enter path to SCIP Optimization Suite tarball: " SCIPtar
  done
  SCIPtar=$(realpath $SCIPtar)
  printf " SCIP .tar found."
  
  if [[ -d ${SCIPtar%.tgz} ]]; then
    printf " Removing old folder: '${SCIPtar%.tgz}'\n"
    rm -r ${SCIPtar%.tgz}
  fi

  printf "${B}Unpacking SCIP Optimization Suite${W}"
  tar xvzf $SCIPtar > /dev/null 2>&1
  rm $SCIPtar
  cd ${SCIPtar%.tgz}
  printf "Done.\n"

  while true; do
      printf "${B}Install the SCIP Optimization Suite using...${W}"
      printf "(1) CMake [recommended for beginners] or\n"
      printf "(2) Makefile [for advanced testing]\n"
      read -p "Please enter how you want to install it. [1/2] " yn
      case $yn in
          [1]* ) fcmake;
                  break;;
          [2]* ) fmake;
                  break;;
          * ) printf "Please answer 1 or 2.";;
      esac
  done
}


while true; do
    printf "${B}The SCIP Optimization Suite will be installed into the current folder.${W}"
    read -p "Do you wish to continue? [Y/n] " yn
    case $yn in
        [Yy]* ) install;
                break;;
        [Nn]* ) exit;;
        * ) printf "Please answer yes or no.";;
    esac
done
