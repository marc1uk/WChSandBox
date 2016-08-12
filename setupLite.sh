#!/bin/bash

#source ./Source_At_Start_T2HK.sh
#** Almost entirely redundant: replace with **
curdir=${PWD}
packagedir=$(dirname ${BASH_SOURCE})
parentdir=/home/marc/LinuxSystemFiles/
export parentdir
#** End replacement **


# for shared libraries
export LD_LIBRARY_PATH=${parentdir}/WChSandBox/WChSandBox_v0:$LD_LIBRARY_PATH

export G4WORKDIR=${parentdir}/WChSandBox/WChSandBox_v0

export WCSIM=${parentdir}/WChSandBox/WChSandBox_v0/src/

