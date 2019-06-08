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
// $Id: mrdPMThit.cc $
//
//
//
#include "mrdPMThit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"

G4Allocator<mrdPMThit> mrdPMThitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//mrdPMThit::mrdPMThit() : fPmtNumber(-1),fPhysVol(0),fDrawit(false) {} // not used. fPhotons(0),

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mrdPMThit::mrdPMThit(const G4Step* aStep){
	//G4cout<<"MRD Photon Hit Step constructor"<<G4endl;
	G4StepPoint* postStepPoint= aStep->GetPostStepPoint();
	fPos = postStepPoint-> GetPosition();
	fTime = postStepPoint-> GetGlobalTime();
	G4StepPoint* preStepPoint= aStep->GetPreStepPoint();
	G4Track* aTrack = aStep->GetTrack();
	fTrackID = aTrack->GetTrackID();
	fParentID = aTrack->GetParentID();
	const G4VProcess* creationproc = aTrack->GetCreatorProcess();
	if(creationproc){fCreationProcess = (const G4String)(creationproc->GetProcessTypeName(creationproc->GetProcessType()));}
	else {fCreationProcess = "Primary";}
	//G4cout<<"creator process is "<<fCreationProcess<<G4endl;
	G4double KE = preStepPoint->GetKineticEnergy();
	fhitWavelength = (4.13566733*2.99792458 *0.0001 )/KE;
	//fPhotons=0;
	fDrawit=false;
	//G4cout<<"step constructor end."<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mrdPMThit::~mrdPMThit() {} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

mrdPMThit::mrdPMThit(const mrdPMThit &right) : G4VHit()
{
  fHitID=right.fHitID;
  fTrackID=right.fTrackID;
  fParentID=right.fParentID;
  fPmtNumber=right.fPmtNumber;
  //fPhotons=right.fPhotons;
  fPos=right.fPos;
  fTime=right.fTime;
  fPhysVol=right.fPhysVol;
  fDrawit=right.fDrawit;
  fCreationProcess=right.fCreationProcess;
  fhitWavelength=right.fhitWavelength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const mrdPMThit& mrdPMThit::operator=(const mrdPMThit &right){
  fHitID=right.fHitID;
  fTrackID=right.fTrackID;
  fParentID=right.fParentID;
  fPmtNumber=right.fPmtNumber;
  //fPhotons=right.fPhotons;
  fPos=right.fPos;
  fTime=right.fTime;
  fPhysVol=right.fPhysVol;
  fDrawit=right.fDrawit;
  fCreationProcess=right.fCreationProcess;
  fhitWavelength=right.fhitWavelength;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int mrdPMThit::operator==(const mrdPMThit &right) const{
  return (fHitID==right.fHitID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMThit::Draw(){
  if(fDrawit&&fPhysVol){ //ReDraw only the PMTs that have hit counts > 0
    //Also need a physical volume to be able to draw anything
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){//Make sure that the VisManager exists
      G4VisAttributes attribs(G4Colour(1.,0.,0.));
      attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(fPhysVol->GetRotation())//If a rotation is defined use it
        rot=*(fPhysVol->GetRotation());
      G4Transform3D trans(rot,fPhysVol->GetTranslation());//Create transform
      pVVisManager->Draw(*fPhysVol,attribs,trans);//Draw it
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void mrdPMThit::Print() {}
