#!/bin/bash
curdir=${PWD}
packagedir=$(dirname ${BASH_SOURCE})
#parentdir=$(cd "${packagedir}/../.."; pwd)
parentdir=/home/marc/LinuxSystemFiles/GEANT4
export parentdir

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#unset GEANT4_BASE_DIR
#export GEANT4_BASE_DIR=${parentdir}/GEANT4/

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#if [[ -e "${GEANT4_BASE_DIR}/install/bin/geant4.sh" ]]; then
#    source "${GEANT4_BASE_DIR}/install/bin/geant4.sh"
#fi

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#if [[ ":${LD_LIBRARY_PATH}:" != "" ]]; then
#    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${G4LIB}/${G4SYSTEM}
#fi

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#unset G4NEUTRONHPDATA
#export G4NEUTRONHPDATA=${GEANT4_BASE_DIR}/data/G4NDL3.14

# DONE IN BASH.RC
#unset ROOTSYS
#export ROOTSYS=~/LinuxSystemFiles/ROOT/build/

# DONE IN BASH.RC
#if [[ -e "${ROOTSYS}/bin/thisroot.sh" ]]; then
#    source ${ROOTSYS}/bin/thisroot.sh
#else
#    echo "Cannot source bin/thisroot.sh"
#fi

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#unset CLHEP_BASE_DIR
#export CLHEP_BASE_DIR=${parentdir}/CLHEP/install

# DONE IN GEANT4 ENV.SH (CALLED IN BASH.RC)
#if [[ ":${LD_LIBRARY_PATH}:" != ":${CLHEP_BASE_DIR}/lib:" ]]; then
#    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CLHEP_BASE_DIR}/lib
#fi

# DON'T HAVE THEM OR KNOW WHAT THEY'RE FOR
#unset UTILSDIR
#export UTILSDIR=${parentdir}/hk-utils
#unset NEUTDIR
#export NEUTDIR=${parentdir}/hk-neut
#unset FITQUN_ROOT
#export FITQUN_ROOT=${parentdir}/hk-fitqun

