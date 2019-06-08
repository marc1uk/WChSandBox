#!/bin/bash
curdir=${PWD}
packagedir=$(dirname ${BASH_SOURCE})
parentdir=$(cd "${packagedir}/../.."; pwd)
export parentdir

unset GEANT4_BASE_DIR
#export GEANT4_BASE_DIR=${parentdir}/GEANT4/
export GEANT4_BASE_DIR=/usr/local/

#if [[ -e "${GEANT4_BASE_DIR}/install/bin/geant4.sh" ]]; then
#    source "${GEANT4_BASE_DIR}/install/bin/geant4.sh"
#fi

if [[ ":${LD_LIBRARY_PATH}:" != "" ]]; then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}
fi

#unset G4NEUTRONHPDATA
#export G4NEUTRONHPDATA=${GEANT4_BASE_DIR}/data/G4NDL3.14
unset G4WORKDIR
export G4WORKDIR="${parentdir}/WChSandBox/WChSandBox_v0/tmp/Linux-g++/WChSandBox/exe"
unset ROOTSYS
#export ROOTSYS=${parentdir}/ROOT/build/
export ROOTSYS=/usr/local/root/

if [[ -e "${ROOTSYS}/bin/thisroot.sh" ]]; then
    source ${ROOTSYS}/bin/thisroot.sh
else
    echo "Cannot source bin/thisroot.sh"
fi

#unset CLHEP_BASE_DIR
#export CLHEP_BASE_DIR=${parentdir}/CLHEP/install

#if [[ ":${LD_LIBRARY_PATH}:" != ":${CLHEP_BASE_DIR}/lib:" ]]; then
#    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
#fi


#unset UTILSDIR
#export UTILSDIR=${parentdir}/hk-utils
#unset NEUTDIR
#export NEUTDIR=${parentdir}/hk-neut
unset WCSIMDIR
export WCSIMDIR=${parentdir}/WChSandBox/WChSandBox_v0
#unset FITQUN_ROOT
#export FITQUN_ROOT=${parentdir}/hk-fitqun

