# ts-WChSandBox
WChSandBox package (core in the ANNIE repository)

WChSanBox is a GEANT based Simulation package for the TITUS code. 

****************************************************************
ts-WChSanBox Usage:
****************************************************************
The package is run by executing the binary file bin/Linux-g++/WChSandBox and passing a configuration (.mac) file as an argument.

eg.

./bin/Linux-g++/WChSandBox MainRun.mac

*****************************************************************
*****************************************************************
Configuration file (.mac): 
*****************************************************************

The .mac files contain configuration settings for the generator run.

eg:"

/control/verbose 0

/mygen/vecfile fluxesandtables/numu_center.txt

/tracking/verbose 0
/tracking/verbose 0
/control/verbose 0
/particle/process/verbose 0
/process/verbose 0
/material/verbose 0
/vis/verbose 0
/run/verbose 0

/run/beamOn 200

"

The verbose statements are bool flags for printouts of each step in the simulation. 0 = no print out, 1 = print out.


The type of particle an location of interaction to be simulated are determined by a vecfile at the line /mygen/vecfile. This points to a vectfile that is stored in the fluxesandtables directory.


The number of events to be simulated are determined by the int argument /run/beamon
 
*************************************************************************
*************************************************************************
Output files:
*************************************************************************

FullEvent.root
generatorcardfile.root
TextOut.txt

descriptions of output to follow....






************************************************************************
Written by Benjamin Richards b.richards@qmul.ac.uk