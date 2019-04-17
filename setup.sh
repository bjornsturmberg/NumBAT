#!/bin/bash

function makeBAT() {
## install dependencies
  apt-get install -y python3-numpy python3-dev python3-scipy python3-nose python3-pip gfortran make gmsh libatlas-dev libblas-dev liblapack-dev libsuitesparse-dev
  pip3 install matplotlib # safer option than apt-get'ing as will install mpl2.0 without conflicting older versions.
## compile Fortran routines
  cd backend/fortran/
  make
## run tests
  echo ''
  echo '#####################################################################'
  echo '        Executing NumBAT test suite'
  echo '#####################################################################'

  cd ../../tests/
  nosetests3

##
  echo ''
  echo '#####################################################################'
  echo '        NumBAT and its dependencies have been installed        '
  echo '               and its test cases calculated.  '
  echo '        Tests outputs are: . passed; E errors; F failed. '
  echo ''
  echo 'NumBAT is brought to you by Bjorn Sturmberg, Kokou Dossou,'
  echo 'Chris Poulton and Michael Steel, with support from CUDOS'
  echo '#####################################################################'
  echo ''

}

makeBAT

