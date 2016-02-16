#!/bin/bash

source /grid/fermiapp/products/common/etc/setup
export PRODUCTS=${PRODUCTS}:/grid/fermiapp/products/uboone
setup root 
source /annie/app/externals/geant4.9.6.p04-build/share/Geant4-9.6.4/geant4make/geant4make.sh
export G4WORKDIR=$PWD
export LD_LIBRARY_PATH=${G4WORKDIR}:$LD_LIBRARY_PATH
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1




