sudo apt-get install gcc gpp git cmake build-essential libgmp-dev libreadline-dev zlib1g-dev bison flex libncurses-dev
git clone git@git.or.rwth-aachen.de:gcg/gcg.git $1
cd $1
git submodule init
git submodule sync
git submodule update
sudo make deps
sudo make
