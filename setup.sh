#!/bin/bash

function makeBAT() {
## compile Fortran routines
  cd backend/fortran/
  make
## run tests
  cd ../../tests/
  nosetests

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

