#!/bin/bash

# root version 5.30
source /myROOTpath/root/bin/thisroot.sh

# GEANT 4.9.4.p02
source /myGEANT4path/geant4.9.4.p02/env.sh

# for shared libraries
export LD_LIBRARY_PATH=/pathofmyworkingarea/WCSimLite:$LD_LIBRARY_PATH

export G4WORKDIR=/pathofmyworkingarea/WCSimLite

export WCSIM=/pathofmyworkingarea/WCSimLite/src/
