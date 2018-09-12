#!/bin/bash

if [ "$TRAVIS_OS_NAME" != "osx" ]; then #
  pushd ..
  wget http://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.12/src/hdf5-1.8.12.tar.gz
  tar -xzf hdf5-1.8.12.tar.gz
  cd hdf5-1.8.12
  ./configure --prefix=/usr/local
  sudo make install
  popd
fi