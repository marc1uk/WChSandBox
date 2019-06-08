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
const int knpartmax = 100000;
const int kmrdhitnmax = 1000; 
const int kcapmax = 100;
const int kpmthitnmax = 100000;

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
    G4int ConvertProcessNameToCode(TString processName);
    G4int ConvertParticleNameToCode(TString particleName);
 


  // out text file
  std::fstream* textout;
  std::ofstream* unaccountedparticlesandprocesses;

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
//  G4int *phot_mrdhit;
//  G4int *phot_mrdprimary;
//  G4int *phot_mrdnumsecs;
//  G4int *phot_mrdpritrackid;
//  G4double *phot_mrdedep;
//  G4double *phot_mrdstartx;
//  G4double *phot_mrdstarty;
//  G4double *phot_mrdstartz;
//  G4int *phot_mrddetected;

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
  G4int *part_GdOriginalTrackID;
  G4int *part_passThruNCV;
  G4double *part_NCVentryX;
  G4double *part_NCVentryY;
  G4double *part_NCVentryZ;
  G4double *part_NCVentryT;
  G4double *part_NCVentryE;
  G4double *part_NCVexitX;
  G4double *part_NCVexitY;
  G4double *part_NCVexitZ;
  G4double *part_NCVexitT;
  G4double *part_NCVexitE;
//  G4int *part_mrdhit;
//  G4int *part_mrdprimary;
//  G4int *part_mrdnumsecs;
//  G4int *part_mrdpritrackid;
//  G4double *part_mrdedep;
//  G4double *part_mrdstartx;
//  G4double *part_mrdstarty;
//  G4double *part_mrdstartz;
//  G4int *part_mrddetected;

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
  
  G4double *mrdhit_x;
  G4double *mrdhit_y;
  G4double *mrdhit_z;
  G4double *mrdhit_t;
  G4int *mrdhit_process;
  G4int *mrdhit_particleID;
  G4int *mrdhit_trackID;
  G4double *mrdhit_edep;
  G4int *mrdhit_copynum;
  G4int *mrdhit_objnum;
  
  G4double *facchit_x;
  G4double *facchit_y;
  G4double *facchit_z;
  G4double *facchit_t;
  G4int *facchit_process;
  G4int *facchit_particleID;
  G4int *facchit_trackID;
  G4double *facchit_edep;
  G4int *facchit_copynum;
  G4int *facchit_objnum;
 
  TFile *mrdfile;
  TTree *mrdtree;
  TTree *facctree;
  TTree *ncvtree;
  G4int mrd_numhits;
  G4int facc_numhits;
  G4int ncv_numhits;
 
  TTree *mrdpmttree;
  TTree *faccpmttree;
  G4int mrdpmt_numhits;
  G4int faccpmt_numhits;
  
  G4double *mrdpmthit_x;
  G4double *mrdpmthit_y;
  G4double *mrdpmthit_z;
  G4double *mrdpmthit_t;
  G4int *mrdpmthit_process;
  G4int *mrdpmthit_trackID;
  G4int *mrdpmthit_parentID;
  G4double *mrdpmthit_wavelength;
  G4int *mrdpmthit_copynum;
  
  G4double *faccpmthit_x;
  G4double *faccpmthit_y;
  G4double *faccpmthit_z;
  G4double *faccpmthit_t;
  G4int *faccpmthit_process;
  G4int *faccpmthit_trackID;
  G4int *faccpmthit_parentID;
  G4double *faccpmthit_wavelength;
  G4int *faccpmthit_copynum;
  
  G4double *ncvhit_x;
  G4double *ncvhit_y;
  G4double *ncvhit_z;
  G4double *ncvhit_t;
  G4int *ncvhit_process;
  G4int *ncvhit_trackID;
  G4int *ncvhit_parentID;
  G4int *ncvhit_particleID;
  G4double *ncvhit_edep;

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
  
  int mrdhit, ismrdprimary, mrdnumsecs, mrdpritrackid, hitTrackID, objnum, copynum, mrddetected, hitPMTnumber, hitParentID;
  double mrdedep, mrdstartx, mrdstarty, mrdstartz, hitWavelength;
  G4double hitPosx, hitPosy, hitPosz, hitTime, hitEdep;
  G4int hitPartCode, hitProcessCode, hitCopyNum;
  G4String hitParticleName, hitProcessName; 
  std::string hitPhysical;
  
  int GdOriginalTrackID, passThruNCV;
  double ncvEntryX, ncvEntryY, ncvEntryZ, ncvEntryT, ncvEntryE, ncvExitX, ncvExitY, ncvExitZ, ncvExitT, ncvExitE;
  

  
  //G4int ConvertParticleNameToCode(G4String particleName);
  //G4int ConvertProcessNameToCode(G4String processName);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
