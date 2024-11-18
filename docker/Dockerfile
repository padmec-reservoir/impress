FROM phusion/baseimage:focal-1.2.0
LABEL author="Filipe Cumaru <facsa@cin.ufpe.br>"

ENV HOME /home/scientific

ENV MOAB_DIR /usr
ENV MPI_HOME /usr
ENV PYTHONPATH /usr/lib/python3.8/site-packages
ENV PATH /usr/local/bin:$PATH
ENV LANG C.UTF-8

RUN set -x \
    && apt-get update \ 
    && apt-get -qy install \
        wget pkg-config git libopenblas-dev liblapack-dev \
        make cmake autoconf automake libtool \
        clang gcc g++ gfortran \
        libhdf5-mpich-dev libnetcdf-c++4 \
        libeigen3-dev libmetis-dev doxygen \
        netgen libnglib-dev  \
        build-essential bzip2 tar m4 file swig unzip libbz2-dev liblzma-dev \
		tcl tk libssl-dev \
		python3.8 python3-pip \
    && apt-get clean

# make some useful symlinks that are expected to exist
RUN cd /usr/bin \
	&& ln -s pydoc3 pydoc \
	&& ln -s python3 python \
	&& ln -s python3-config python-config


RUN pip install cython==0.29.32 numpy==1.23.1 \
	pytest==7.1.2 colorlog==6.6.0 configobj==5.0.6 \
	pytest-runner==6.0.0 ipython==8.4.0 pyyaml==5.3.1 \
	scipy==1.9.0 jedi==0.17.2

# # Install modified version of MOAB and PyMOAB.
WORKDIR $HOME
ADD ./docker/moab.zip $HOME
RUN cd $HOME \
	&& unzip moab.zip \
    && cd $HOME/moab \
    && autoreconf -fi \
    && mkdir build \
    && cd $HOME/moab/build \
    && ../configure \
            --prefix=$MOAB_DIR \
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
    && make -j4\
    && make install \
	&& cd ${HOME}/moab \
	&& git clone https://github.com/gabrielmmats/PymoabIMPRESS.git \
    && mv PymoabIMPRESS ${HOME}/moab/build/pymoab \
    && cd ${HOME}/moab/build/pymoab \
    && rm -rf pymoab \
    && mv PymoabIMPRESS ./pymoab \
    && python setup.py build \
    && python setup.py install \
    && cd ${HOME} \
    && rm -rf moab/ \
    && echo "export MOAB_DIR=$MOAB_DIR" >> $HOME/.bashrc \
    && echo "export PATH=$PATH:$MOAB_DIR/bin" >> $HOME/.bashrc \
    && echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MOAB_DIR/lib" >> $HOME/.bashrc

RUN git clone https://github.com/padmec-reservoir/impress.git \
	&& echo "export PYTHONPATH=${PYTHONPATH}:${HOME}/impress" >> ${HOME}/.bashrc
