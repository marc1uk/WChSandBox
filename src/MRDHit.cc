// ====================================================================
//   MRDHit.cc
//
//   26/11/15 M. O'Flaherty (mostly copied from SBsimMRDHit; 2006/03/03 K. Hiraide)
// ====================================================================
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "MRDHit.hh"

#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"

#include <iostream>

// allocator
G4Allocator<MRDHit> MRDHitAllocator;

//////////////////////////
MRDHit::MRDHit(){}

//////////////////////////
MRDHit::MRDHit(G4Step* aStep) {

	G4StepPoint* preStepPoint= aStep->GetPreStepPoint();
	G4Track* aTrack = aStep->GetTrack();
	G4TouchableHistory* touchable= (G4TouchableHistory*)(preStepPoint-> GetTouchable());    
  	G4VPhysicalVolume* thePhysical = touchable->GetVolume();

	hitTrackID = aTrack->GetTrackID();
	hitPos = preStepPoint-> GetPosition();
	hitTime = preStepPoint-> GetGlobalTime();
	hitParticleName = aTrack-> GetDefinition()->GetParticleName();
	hitParticleID = aTrack-> GetDefinition()->GetPDGEncoding(); 
	hitProcessName = ((G4VProcess*) (aStep->GetPostStepPoint()->GetProcessDefinedStep()))->GetProcessName();
	eDeposited = aStep-> GetTotalEnergyDeposit();
	hitDetectTime = 0;
	hitCopyNum = thePhysical->GetCopyNo();
	//hitPanelNum;
	//hitPaddleNum;

	/*
	//  G4double              kineticE = preStepPoint->GetKineticEnergy();
 	G4double              kineticE = aTrack->GetKineticEnergy();
  	G4Material*           material = aTrack->GetMaterial();
  	*/
}
  	
///////////////////////////
MRDHit::~MRDHit(){}

////////////////////////
void MRDHit::Draw()
////////////////////////
{
  /*
  G4VVisManager* pVVisManager= G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(hitpos);
    circle.SetScreenSize(5.0);
    circle.SetFillStyle(G4Circle::filled);

    G4Color color, goodColor(1.,0.,0.), badColor(0.,0.,1.);

    color=goodColor;

    G4VisAttributes attribs(color);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
  */
  /*
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Transform3D trans(rot.inverse(),hitpos);
      G4VisAttributes attribs;
      const G4VisAttributes* pVA = pLogV->GetVisAttributes();
      if(pVA) attribs = *pVA;
      G4Colour colour(1.,0.,0.);
      attribs.SetColour(colour);
      attribs.SetForceWireframe(false);
      attribs.SetForceSolid(true);
      pVVisManager->Draw(*pLogV,attribs,trans);
    }
  */
}

/////////////////////////
void MRDHit::Print() {
  G4cout.precision(2);	//?
  G4cout << "hitID= "    << hitID << ", "
	 << "trackID= " << hitTrackID   << "("
	 << hitParticleName << ", " << hitProcessName << "), "
         << "t= "    << hitTime/ns     << " ns, ";
  G4cout.precision(0);
  G4cout << "pos = "  << hitPos*(1./cm) << " cm,";
  G4cout << "hitcopynum = " << hitCopyNum << G4endl;

}



