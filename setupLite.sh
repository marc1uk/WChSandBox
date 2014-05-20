#!/bin/bash

source ../hk-hyperk/Source_At_Start_T2HK.sh

# for shared libraries
export LD_LIBRARY_PATH=${parentdir}/hk-wchsandbox:$LD_LIBRARY_PATH

export G4WORKDIR=${parentdir}/hk-wchsandbox

export WCSIM=${parentdir}/hk-wchsandbox/src/
