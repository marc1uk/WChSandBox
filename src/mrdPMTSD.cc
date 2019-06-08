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
// $Id: mrdPMTSD.cc  $
//
//
//
#include "mrdPMTSD.hh"
#include "mrdPMThit.hh"
#include "WCLiteDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mrdPMTSD::mrdPMTSD(G4String name)
  : G4VSensitiveDetector(name),fPMTHitCollection(0)
{
  collectionName.insert("mrdpmtHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mrdPMTSD::~mrdPMTSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMTSD::Initialize(G4HCofThisEvent* hitsCE){
  fPMTHitCollection = new mrdPMThitsCollection(SensitiveDetectorName,collectionName[0]);
  //Store collection with event and keep ID
  static G4int hitCID = -1;
  if(hitCID<0){
    hitCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitCID, fPMTHitCollection );
  nhit=0;
  G4cout<<"Initializing LG/PMT sensitive detectors."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool mrdPMTSD::ProcessHits(G4Step* ,G4TouchableHistory* ){
  return false;	
  G4cout<<"rejected hit on LG."<<G4endl;
  // don't handle hits automatically - we'll call our own processhits manually for optical photons on boundary steps.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode

G4bool mrdPMTSD::ProcessHits_PhotonHit(const G4Step* aStep, G4TouchableHistory* ){
  //G4cout<<"*************Manual call to process hits!*********"<<G4endl;

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  G4int pmtNumber = aStep->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber();
  G4VPhysicalVolume* physVol = aStep->GetPostStepPoint()->GetPhysicalVolume();

// Default behaviour from LXe example is to have just one hit per PMT, 
// each of which just stores the TOTAL number of photons on that PMT. 
// This gives no information about individual hits (e.g. position).
// -------------------------------------------------------- 
/*  //Find the correct hit collection
  G4int n=fPMTHitCollection->entries();
  mrdPMThit* hit=NULL;
  for(G4int i=0;i<n;i++){
    if((*fPMTHitCollection)[i]->GetPMTNumber()==pmtNumber){
      hit=(*fPMTHitCollection)[i];
      break;
    }
  }
 
  if(hit==NULL){//this pmt wasnt previously hit in this event
    hit = new mrdPMThit(); //so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetPMTPhysVol(physVol);
    fPMTHitCollection->insert(hit);
    hit->SetHitPos(aStep->GetPostStepPoint()->GetPosition());
  }

  hit->IncPhotonCount(); //increment hit for the selected pmt
*/
// -------------------------------------------------------- 
// Gonna instead use the same mechanism as for previous sensitive detector, 
// so we can record full info for each hit.
// --------------------------------------------------------
     mrdPMThit* hit = new mrdPMThit(aStep); // create new hit
     hit->SetHitID(nhit);
     nhit++;
     hit->SetPMTNumber(pmtNumber);
     hit->SetPMTPhysVol(physVol);
     fPMTHitCollection->insert(hit);
     
/*  if(!LXeDetectorConstruction::GetSphereOn()){
    hit->SetDrawit(true);
    //If the sphere is disabled then this hit is automaticaly drawn
  }
  else{//sphere enabled
    LXeUserTrackInformation* trackInfo=
      (LXeUserTrackInformation*)aStep->GetTrack()->GetUserInformation();
    if(trackInfo->GetTrackStatus()&hitSphere)
      //only draw this hit if the photon has hit the sphere first
      hit->SetDrawit(true);
  }
*/

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMTSD::EndOfEvent(G4HCofThisEvent* ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMTSD::clear() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMTSD::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMTSD::PrintAll() {}
