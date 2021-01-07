FROM phusion/baseimage:master
LABEL author="Filipe Cumaru <facsa@cin.ufpe.br>"

ENV HOME /home/scientific

ENV MOAB_DIR /usr
ENV MPI_HOME /usr
ENV PYTHONPATH /usr/lib/python3.6/site-packages
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
        liboce-foundation-dev liboce-modeling-dev liboce-ocaf-dev liboce-visualization-dev oce-draw \
        netgen libnglib-dev  \
        build-essential bzip2 tar m4 file swig unzip libbz2-dev liblzma-dev \
	tcl tk libssl-dev \
    && apt-get clean


ENV GPG_KEY 0D96DF4D4110E5C43FBFB17F2D347EA6AA65421D
ENV PYTHON_VERSION 3.6.3


# Install Python3
RUN set -ex \
	&& buildDeps=' \
		dpkg-dev \
		tcl-dev \
		tk-dev \
	' \
	&& apt-get update && apt-get install -y $buildDeps --no-install-recommends && rm -rf /var/lib/apt/lists/* \
	\
	&& wget -O python.tar.xz "https://www.python.org/ftp/python/${PYTHON_VERSION%%[a-z]*}/Python-$PYTHON_VERSION.tar.xz" \
	&& wget -O python.tar.xz.asc "https://www.python.org/ftp/python/${PYTHON_VERSION%%[a-z]*}/Python-$PYTHON_VERSION.tar.xz.asc" \
	&& export GNUPGHOME="$(mktemp -d)" \
	&& gpg --keyserver ha.pool.sks-keyservers.net --recv-keys "$GPG_KEY" \
	&& gpg --batch --verify python.tar.xz.asc python.tar.xz \
	&& rm -rf "$GNUPGHOME" python.tar.xz.asc \
	&& mkdir -p /usr/src/python \
	&& tar -xJC /usr/src/python --strip-components=1 -f python.tar.xz \
	&& rm python.tar.xz \
	\
	&& cd /usr/src/python \
	&& gnuArch="$(dpkg-architecture --query DEB_BUILD_GNU_TYPE)" \
	&& ./configure \
		--build="$gnuArch" \
		--enable-loadable-sqlite-extensions \
		--enable-shared \
		--with-system-expat \
		--with-system-ffi \
		--without-ensurepip \
	&& make -j "$(nproc)" \
	&& make install \
	&& ldconfig \
	\
	&& apt-get purge -y --auto-remove $buildDeps \
	\
	&& find /usr/local -depth \
		\( \
			\( -type d -a \( -name test -o -name tests \) \) \
			-o \
			\( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
		\) -exec rm -rf '{}' + \
	&& rm -rf /usr/src/python


# make some useful symlinks that are expected to exist
RUN cd /usr/local/bin \
	&& ln -s idle3 idle \
	&& ln -s pydoc3 pydoc \
	&& ln -s python3 python \
	&& ln -s python3-config python-config


# if this is called "PIP_VERSION", pip explodes with "ValueError: invalid truth value '<VERSION>'"

ENV PYTHON_PIP_VERSION 9.0.1


RUN set -ex; \
	\
	wget -O get-pip.py 'https://bootstrap.pypa.io/get-pip.py'; \
	\
	python get-pip.py \
		--disable-pip-version-check \
		--no-cache-dir \
		"pip==$PYTHON_PIP_VERSION" \
	; \
	pip --version; \
	\
	find /usr/local -depth \
		\( \
			\( -type d -a \( -name test -o -name tests \) \) \
			-o \
			\( -type f -a \( -name '*.pyc' -o -name '*.pyo' \) \) \
		\) -exec rm -rf '{}' +; \
	rm -f get-pip.py; \
	pip install -U pip



CMD ["/sbin/my_init"]

RUN pip install cython numpy pytest colorlog configobj pytest-runner ipython pyyaml numba scipy pandas jedi==0.17.2

# Install mpi4py
WORKDIR $HOME
RUN cd $HOME \
    && git clone https://bitbucket.org/mpi4py/mpi4py \
    && cd $HOME/mpi4py \
    && python setup.py build \
    && python setup.py install \
    && cd $HOME \
    && rm -rf mpi4py


# Install modified version of MOAB and PyMOAB.
WORKDIR $HOME
ADD ./moab.zip $HOME
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
