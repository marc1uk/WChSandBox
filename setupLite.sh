#!/bin/bash

source ../hk-hyperk/Source_At_Start_T2HK.sh

# for shared libraries
export LD_LIBRARY_PATH=${parentdir}/annie-WCSimLite:$LD_LIBRARY_PATH

export G4WORKDIR=${parentdir}/annie-WCSimLite

export WCSIM=${parentdir}/annie-WCSimLite/src/
