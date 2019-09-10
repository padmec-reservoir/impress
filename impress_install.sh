echo "####################################"
echo "STEP 1: INSTALLING MOAB DEPENDENCIES"
echo "####################################"
sudo apt-get -y update &&
sudo apt-get -y install wget pkg-config git libopenblas-dev liblapack-dev make cmake autoconf automake libtool clang gcc g++ gfortran libhdf5-mpich-dev libnetcdf-c++4 libeigen3-dev libmetis-dev doxygen liboce-foundation-dev liboce-modeling-dev liboce-ocaf-dev liboce-visualization-dev oce-draw build-essential bzip2 tar m4 file swig tcl tk libssl-dev python3-pip python-pip cython &&
sudo apt-get clean &&
sudo pip3 install mpi4py cython numpy pytest colorlog configobj pytest-runner && sudo pip install cython

echo "#######################"
echo "STEP 2: INSTALLING MOAB"
echo "#######################"

export BASE_INSTALL_DIR=$PWD &&
sudo git clone https://bitbucket.org/fathomteam/moab.git &&
    cd moab/ &&
    sudo autoreconf -fi &&
    sudo mkdir build &&
    cd $BASE_INSTALL_DIR/moab/build &&
    sudo ../configure \
            --prefix=$BASE_INSTALL_DIR/MOAB \
            --with-mpi=/usr \
            --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/mpich \
            --with-netcdf=/usr \
            --with-metis=/usr \
            --enable-optimize \
            --enable-debug \
            --enable-tools \
            --enable-pymoab \
            --enable-shared \
            CFLAGS='-O2 -fPIC -DPIC' \
            CXXFLAGS='-O2 -fPIC -DPIC' \
            FCFLAGS='-O2 -fPIC' \
            FFLAGS='-O2 -fPIC' \
            PYTHON=python3 &&
    sudo make -j4 &&
    sudo make install &&
    cd $BASE_INSTALL_DIR/moab/build/pymoab &&
    sudo python3 setup.py build
    sudo python3 setup.py install &&
    cd $BASE_INSTALL_DIR &&
    sudo echo "export PATH=$PATH:$HOME/MOAB/bin" >> $HOME/.bashrc &&
    sudo echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/MOAB/lib" >> $HOME/.bashrc &&
    sudo echo "export CPATH=$CPATH:$HOME/MOAB/include" >> $HOME/.bashrc &&
    sudo echo "export PYTHONPATH=$PYTHONPATH:$HOME/MOAB/lib/python3.6/" >> $HOME/.bashrc

echo "###########################################"
echo "STEP 3: UPDATING PYMOAB TO IMPRESS VERSION"
echo "###########################################"

cd moab/pymoab &&
sudo git clone https://github.com/gabrielmmats/PymoabIMPRESS &&
sudo rm -rf pymoab &&
sudo mv PymoabIMPRESS pymoab &&
cd ../build/pymoab &&
sudo python3 setup.py build &&
sudo python3 setup.py install &&
cd ../../.. &&
sudo rm -rf moab

echo "##########################"
echo "STEP 4: INSTALLING IMPRESS"
echo "##########################"

sudo pip3 install impress-padmec
