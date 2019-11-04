RED='\033[0;31m'
W='\e[0m\n'
B='\n\e[1m'

if [ ! -z $1 ];
  then
    NAME=$1
  else
    echo "Usage: ./installGCGFromGit.sh <folder name> <branch>"
    exit
fi

if [ ! -z $2 ];
  then
    BRANCH=$2
  else
    BRANCH=master
fi

function start(){
  while true; do
    printf "${B}Do you want to start GCG?${W}"
    read -p "Please answer. [Y/n] " yn
    case $yn in
      [Yy]* )
              ./bin/gcg
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
              sudo apt-get install gcc gpp git cmake build-essential libgmp-dev libreadline-dev zlib1g-dev bison flex libncurses-dev
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
  git clone git@git.or.rwth-aachen.de:gcg/gcg.git $NAME
  cd $NAME
  git checkout $BRANCH

  # Get SCIP and SoPlex from GCG git
  printf "${B}Cloning submodules (SCIP and SoPlex)...${W}"
  git submodule init
  git submodule sync
  git submodule update
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

function install(){
  installPrerequisites;
  download;
  setLinks;
  compile;
  test;
  start;
}

while true; do
  printf "${B}GCG, SCIP and SoPlex will be installed from Git into the folder '${1}'.${W}"
  read -p "Do you wish to continue? [Y/n] " yn
  case $yn in
    [Yy]* ) install;
    break;;
    [Nn]* ) exit;;
    * ) printf "Please answer yes or no.";;
  esac
done
