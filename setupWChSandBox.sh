#!/bin/bash

# setup root
source /grid/fermiapp/products/uboone/setup_uboone.sh
setup root

# setup G4parameters
source /annie/app/externals/geant4.9.6.p04-build/share/Geant4-9.6.4/geant4make/geant4make.sh
export G4BIN=/annie/app/externals/geant4.9.6.p04-build/bin/
# for shared libraries
export G4WORKDIR=$PWD
export LD_LIBRARY_PATH=$G4WORKDIR:$LD_LIBRARY_PATH
export WCSIM=$G4WORKDIR/src

export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1

