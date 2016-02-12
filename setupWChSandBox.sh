#!/bin/bash

# setup root
source /grid/fermiapp/products/uboone/setup_uboone.sh
setup root

# setup G4parameters
source /annie/app/users/wetstein/geant4.9.6.p04-build/share/Geant4-9.6.4/geant4make/geant4make.sh

# for shared libraries
export LD_LIBRARY_PATH=/annie/app/users/txin/WChSandBox:$LD_LIBRARY_PATH

export G4WORKDIR=/annie/app/users/txin/WChSandBox

export WCSIM=/annie/app/users/txin/WChSandBox/src/

export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1

