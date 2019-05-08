#!/bin/bash

env | sort



# GB: install xdrfile library
wget http://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
tar xzf xdrfile-1.1.4.tar.gz
(
    cd xdrfile-1.1.4
    ./configure --prefix=$PREFIX --enable-shared
    make
    make install
)



# Thanks to https://github.com/intbio/plumed-conda/blob/master/plumed2_v2.5.0/build.sh
condaldflags="$BUILD_PREFIX/lib/libz.a $BUILD_PREFIX/lib/libz$SHLIB_EXT
	      $BUILD_PREFIX/lib/libxdrfile.a $BUILD_PREFIX/lib/libxdrfile$SHLIB_EXT"

./configure --prefix=$PREFIX --enable-shared --disable-python --disable-external-lapack --disable-external-blas LDFLAGS="$condaldflags"
make -j4
make install

cd python
make pip
export plumed_default_kernel=$PREFIX/lib/libplumedKernel$SHLIB_EXT
$PYTHON -m pip install .

