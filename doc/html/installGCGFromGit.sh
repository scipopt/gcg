if [ ! -z $1 ];
  then
    echo $1
    NAME=$1
  else
    echo "Usage: ./installGCGFromGit.sh <folder name>"
    exit
fi

echo "Installing prerequisites..."
sudo apt-get install gcc gpp git cmake build-essential libgmp-dev libreadline-dev zlib1g-dev bison flex libncurses-dev

echo "Cloning repository..."
git clone git@git.or.rwth-aachen.de:gcg/gcg.git $NAME
cd $NAME

echo "Getting submodules..."
git submodule init
git submodule sync
git submodule update

echo "Setting links..."
sudo ln -sfn ../lib/scip-git lib/scip
# SoPlex links inside SCIP lib folder
sudo ln -sfn ../../../soplex-git/src/ lib/scip/lib/include/spxinc
sudo ln -sfn ../../../soplex-git/lib/libsoplex.linux.x86_64.gnu.opt.a lib/scip/lib/static/libsoplex.linux.x86_64.gnu.opt.a

echo "Compiling dependencies..."
sudo make deps

echo "Compiling GCG..."
sudo make
