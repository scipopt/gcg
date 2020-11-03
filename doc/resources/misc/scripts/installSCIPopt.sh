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
              printf "Running SCIP on some sample instances...";
                if [ $1 = 1 ]
                then
                  sudo make check
                elif [ $1 = 2 ]
                then
                  sudo make test
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
  printf "${B}Compiling SCIP Optimization Suite${W}"
  mkdir build
  cd build
  sudo cmake .. | tail -n +91
  sudo make
  printf "${B}Creating Executable${W}"
  sudo make install
  test "1"
}

function fmake(){
  printf "${B}Updating Make${W}"
  sudo apt-get install make
  printf "${B}Creating a single library containing SCIP, SoPlex and ZIMPL${W}"
  sudo make scipoptlib
  printf "${B}Compiling GCG${W}"
  sudo make gcg
  printf "${B}Creating Executable${W}"
  sudo make install
  test "2"
}

function install(){
  printf "${B}Updating Package Manager...${W}"
  sudo apt-get update
  printf "${B}Installing Prerequisites...${W}"
  sudo apt-get install build-essential libreadline-dev libz-dev libgmp3-dev lib32ncurses5-dev libboost-program-options-dev
  printf "${B}Initializing Installation...${W}"
  while [[ ! -f $SCIPtar || -z $SCIPtar ]]; do
    read -e -p " Please enter path to SCIP Optimization Suite tarball: " SCIPtar
  done
  SCIPtar=$(realpath $SCIPtar)
  printf " SCIP tar found.\n"
  
  if [[ -d ${SCIPtar%.tgz} ]]; then
    printf " Removing old folder: '${SCIPtar%.tgz}'"
    rm -r ${SCIPtar%.tgz}
  fi

  printf "\n${B}Unpacking SCIP Optimization Suite${W}"
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
