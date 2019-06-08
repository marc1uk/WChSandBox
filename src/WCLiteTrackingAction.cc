
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
#include "WCLiteTrajectory.hh"
#include "WCLiteTrackingAction.hh"
//#include "WCLiteUserTrackInformation.hh"
#include "WCLiteDetectorConstruction.hh"
//#include "RecorderBase.hh" 
#include <iostream>
#include <fstream>
#include "G4TrackingManager.hh"
#include "G4Trajectory.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::WCLiteTrackingAction()
{

  callcount==0;

  // Tree Declarations

  HitTreeEnd = new TTree("hittreeend","hittreeend");

  HitTreeEnd->Branch("trackID0",&TrackID_S);
  HitTreeEnd->Branch("pdgid",&partcode);
  HitTreeEnd->Branch("parentID0",&parentID_S);
  HitTreeEnd->Branch("Et0",&TEStart);
  HitTreeEnd->Branch("E0",&KEStart);
  HitTreeEnd->Branch("x0",&xStart);
  HitTreeEnd->Branch("y0",&yStart);
  HitTreeEnd->Branch("z0",&zStart);
  HitTreeEnd->Branch("t0",&tStart);
  HitTreeEnd->Branch("processStart",&processStart);
  HitTreeEnd->Branch("trackIDf",&TrackID_E);
  HitTreeEnd->Branch("parentIDf",&parentID_E);
  HitTreeEnd->Branch("Etf",&TEEnd);
  HitTreeEnd->Branch("Ef",&KEEnd);
  HitTreeEnd->Branch("xf",&xEnd);
  HitTreeEnd->Branch("yf",&yEnd);
  HitTreeEnd->Branch("zf",&zEnd);
  HitTreeEnd->Branch("tf",&tEnd);
  HitTreeEnd->Branch("wavelength",&wavelength);
  HitTreeEnd->Branch("processEnd",&processEnd);
  HitTreeEnd->Branch("phhit",&phhit);
  HitTreeEnd->Branch("isScatPhot",&isScatPhot);


  // Variable Decalarations
  line=0;
  numhits=0;
  numphots=0;
  nummuplus=0;
  nummuminus=0;
  numeplus=0;
  numeminus=0;
  numgamma=0;
  back_count=0;
  Etot=0;
  iter=1;
  iter2=1;

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void WCLiteTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);
  //Use custom trajectory class
  //  fpTrackingManager->SetTrajectory(new WCLiteTrajectory(aTrack));
  
  G4String theProcess;
  G4String theLocation;

  isScatPhot=0.0;

  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
  G4int stepno = (aTrack->GetCurrentStepNumber());
  G4Step* aStep = (G4Step*) (aTrack->GetStep());
 
  theProcess+="UserLimit";

  processStart = -1;
 
  detectorStart = theLocation;

  
  //  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() ){
      WCLiteTrajectory* thisTrajectory = new WCLiteTrajectory(aTrack);
      fpTrackingManager->SetTrajectory(thisTrajectory);
      fpTrackingManager->SetStoreTrajectory(true);
      //    }
      //  else 
      //    fpTrackingManager->SetStoreTrajectory(false);
  
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void WCLiteTrackingAction::PostUserTrackingAction(const G4Track* aTrack){ 

  callcount++;

  phhit=0;

  G4String theProcess;
  G4String theLocation;
  G4String sProcess;
  
  // G4Step* aStep = aStep->GetTrack();


  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
  G4int stepno = (aTrack->GetCurrentStepNumber());
  G4Step* aStep = (G4Step*) (aTrack->GetStep());

  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
    theProcess = ((G4VProcess*) (aStep->GetPostStepPoint()->GetProcessDefinedStep()))->GetProcessName();
  } else {
    theProcess+="UserLimit";
  }

  if(line%100000==0) G4cout<<"  Post tracking call number: "<<line<<" "<<detectorEnd<<G4endl;
  
  if(partname=="opticalphoton") numphots++;
  line++;  

  //  G4cout<<"post-end"<<G4endl;
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::~WCLiteTrackingAction()
{



}
















