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

#ifndef WCLiteTrackingAction_h
#define WCLiteTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"

class RecorderBase;

class WCLiteTrackingAction : public G4UserTrackingAction {

public:  
  WCLiteTrackingAction();
  ~WCLiteTrackingAction();
  
  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);
  
private:
  RecorderBase* recorder;

  // Histogram Declarations

  TTree* HitTreeEnd;

  // Variable Decalarations
  int line;
  int numhits;
  int numphots;
  int nummuplus;
  int nummuminus;
  int numeplus;
  int numeminus;
  int numgamma;

  int back_count;
  int callcount;

  double xStart,yStart,zStart,xDirStart,yDirStart,zDirStart,tStart,KEStart,TEStart,xEnd,yEnd,zEnd,xDirEnd,yDirEnd,zDirEnd,tEnd,KEEnd,TEEnd,partcode,TrackID_S,TrackID_E,parentID_S,parentID_E,phhit,isScatPhot;
  TString partname, detectorStart, detectorEnd; 
  double processEnd,processStart;
  int numsteps;
  double Etot;
  double wavelength;
  double R;

  int iter;
  int iter2;
};

#endif
