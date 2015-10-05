# ts-WChSandBox
WChSandBox package (core in the ANNIE repository)

WChSanBox is a GEANT based Simulation package for the TITUS code. 

****************************************************************
Building ts-WChSandBox
****************************************************************
The package can be built using the using the automated build scripts in ts-titus package.

It can also be built seperately by first sourcing the "Source_At_Start.sh" script in the ts-titus package and then running make.

source ../ts-titus/Source_At_Start.sh
make
   
****************************************************************
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

descriptions of output files (WORK IN PROGRESS!)

FullEvent.root: Tree "EventTree"

   Declaration of leaf types
   Int_t           evt; //No. of event simulated
   Int_t           nphot; //Number of photos per event
   Int_t           npart; //Number of particles per event
   Int_t           ncapturecount; //Number of captured neutrons
   Int_t           neutroncount; // Number of neutrons 
   Double_t        phot_xStart[623422];   //[nphot] //photon start location X
   Double_t        phot_yStart[623422];   //[nphot]  //photon start location Y
   Double_t        phot_zStart[623422];   //[nphot] //photon start location Z
   Double_t        phot_tStart[623422];   //[nphot] //photon start location time
   Double_t        phot_xEnd[623422];   //[nphot] //photon end location X
   Double_t        phot_yEnd[623422];   //[nphot] //photon end location Y
   Double_t        phot_zEnd[623422];   //[nphot] //photon end location Z
   Double_t        phot_tEnd[623422];   //[nphot] //photon end location time
   Double_t        phot_wavelength[623422];   //[nphot] // Wavelength of photon
   Int_t           phot_processStart[623422];   //[nphot] // 
   Int_t           phot_isScat[623422];   //[nphot] // Bool if scattered scatter
   Int_t           phot_parentid[623422];   //[nphot] //
   Int_t           phot_trackid[623422];   //[nphot] // 
   Int_t           phot_hit[623422];   //[nphot] // Bool hit pmt
   Int_t           phot_capnum[623422];   //[nphot] 
   Double_t        part_xStart[3623];   //[npart]
   Double_t        part_yStart[3623];   //[npart]
   Double_t        part_zStart[3623];   //[npart]
   Double_t        part_tStart[3623];   //[npart]
   Double_t        part_xEnd[3623];   //[npart]
   Double_t        part_yEnd[3623];   //[npart]
   Double_t        part_zEnd[3623];   //[npart]
   Double_t        part_tEnd[3623];   //[npart]
   Double_t        part_pxStart[3623];   //[npart] //Momentum start
   Double_t        part_pyStart[3623];   //[npart]
   Double_t        part_pzStart[3623];   //[npart]
   Double_t        part_pxEnd[3623];   //[npart] // momentum end
   Double_t        part_pyEnd[3623];   //[npart]
   Double_t        part_pzEnd[3623];   //[npart]
   Double_t        part_KEstart[3623];   //[npart] //kenetic energy start
   Double_t        part_KEend[3623];   //[npart]
   Int_t           part_processStart[3623];   //[npart]//
   Int_t           part_processEnd[3623];   //[npart]
   Int_t           part_parentid[3623];   //[npart] /// Parent particle
   Int_t           part_trackid[3623];   //[npart] //
   Int_t           part_pid[3623];   //[npart] //
   Double_t        capt_x[14];   //[ncapturecount]
   Double_t        capt_y[14];   //[ncapturecount]
   Double_t        capt_z[14];   //[ncapturecount]




************************************************************************
Written by Benjamin Richards b.richards@qmul.ac.uk