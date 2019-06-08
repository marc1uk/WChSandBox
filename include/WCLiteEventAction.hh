//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: WCLiteEventAction.hh,v 1.2 2006/06/29 17:53:57 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef WCLiteEventAction_h
#define WCLiteEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "TTree.h"
#include "TRandom3.h"


// Modified 2014-06-13 by MSM ?!? 
//  - Added declaration of knphotmax here for general usage
//  - Increased photon buffer by factor 10 (1e6 -> 1e7)
// Remodified 2014-06-16 by MJW ##
// Initialization moved outside of class

const int knphotmax = 10000000;
const int knpartmax = 10000;
const int kcapmax = 100;

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WCLiteEventAction : public G4UserEventAction
{
  public:
    WCLiteEventAction();
   ~WCLiteEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
 


  // out text file
  fstream* textout;

  TFile *no;
  TTree *nt,*ht,*trk,*gtrk,*mul,*capt,*gcapt;
  TTree *evttree;
  
  TRandom3* nR;

  G4int eventcount;
  G4int neutroncount;
  G4int ncapturecount;
  G4int nphot;
  G4int npart;

  G4double *phot_xStart;
  G4double *phot_yStart;
  G4double *phot_zStart;
  G4double *phot_tStart;
  G4double *phot_xEnd;
  G4double *phot_yEnd;
  G4double *phot_zEnd;
  G4double *phot_tEnd;
  G4double *phot_pxStart;
  G4double *phot_pyStart;
  G4double *phot_pzStart;
  G4double *phot_pxEnd;
  G4double *phot_pyEnd;
  G4double *phot_pzEnd;
  G4double *phot_wavelength;
  G4int *phot_capnum;
  G4int *phot_processStart;
  G4int *phot_isScat;
  G4int *phot_parentid;
  G4int *phot_trackid;
  G4int *phot_hit;

  G4double *part_xStart;
  G4double *part_yStart;
  G4double *part_zStart;
  G4double *part_tStart;
  G4double *part_xEnd;
  G4double *part_yEnd;
  G4double *part_zEnd;
  G4double *part_tEnd;
  G4double *part_pxStart;
  G4double *part_pyStart;
  G4double *part_pzStart;
  G4double *part_pxEnd;
  G4double *part_pyEnd;
  G4double *part_pzEnd;
  G4double *part_KEstart;
  G4double *part_KEend;
  G4int *part_processStart;
  G4int *part_processEnd;
  G4int *part_parentid;
  G4int *part_trackid;
  G4int *part_pid;

  G4int *capt_num;
  G4int *capt_nucleus;
  G4int *capt_pid;
  G4int *capt_nphot;
  G4int *capt_ngamma;
  G4double *capt_x;
  G4double *capt_y;
  G4double *capt_z;
  G4double *capt_t0;
  G4double *capt_E;


  int totalcapturecount;

  double partcode,xStart,yStart,zStart,tStart,xStartDir,yStartDir,zStartDir;
  double xEnd,yEnd,zEnd,tEnd,xEndDir,yEndDir,zEndDir,phhit,processEnd;
  double trackid,parentid,processStart,KEStart,TEStart,TEEnd,KEEnd;
  double wavelengthStart,isScatPhot;
  
  double xStep,yStep,zStep,tStep,iStep,processStep,TEStep,KEStep;
  
  double nphotspercapt;
  double cpartcode,cxStart,cyStart,czStart,ctStart;
  double cxEnd,cyEnd,czEnd,ctEnd,cphhit,cprocessEnd;
  double ctrackid,cparentid,cprocessStart,cKEStart,cTEStart,cTEEnd,cKEEnd;
  double cwavelengthStart,cisScatPhot;
  double captenergy,ncaptgammas,captT,captx,capty,captz;

 
  double mhxEnd,mhyEnd,mhzEnd,mhtEnd,mhprocessStart,mhwavelength;
  double totalmucount;
  double mhisScatPhot;
  double mhparentid; 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
