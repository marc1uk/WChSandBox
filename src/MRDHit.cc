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
#include "G4SystemOfUnits.hh"

#include <iostream>

// allocator
G4Allocator<MRDHit> MRDHitAllocator;

//////////////////////////
MRDHit::MRDHit(){}

//////////////////////////
MRDHit::MRDHit(G4Step* aStep) {

	G4StepPoint* preStepPoint= aStep->GetPreStepPoint();
	G4Track* aTrack = aStep->GetTrack();
	G4TouchableHandle touchable = preStepPoint->GetTouchableHandle();	//not called a pointer but USED as one!
	//G4TouchableHistory* touchablehist= (G4TouchableHistory*)(preStepPoint-> GetTouchable());    
  	G4VPhysicalVolume* thePhysical = touchable->GetVolume();
  	//G4VPhysicalVolume* thePhysical2 = touchablehist->GetVolume();
	//G4cout<<"touchable history says volume is:"<<(thePhysical->GetName())<<" while touchablehandle says "<<(thePhysical->GetName())<<G4endl; <-gives same result
	hitTrackID = aTrack->GetTrackID();
	hitPos = preStepPoint-> GetPosition();
	hitTime = preStepPoint-> GetGlobalTime();
	hitParticleName = aTrack-> GetDefinition()->GetParticleName();
	hitParticleID = aTrack-> GetDefinition()->GetPDGEncoding(); 
	//G4cout<<"getting process"<<G4endl;
	const G4VProcess* theProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep();
	if(theProcess){hitProcessName = theProcess->GetProcessName();} else {hitProcessName="UserLimit";}
	// above must use POST step point for process; otherwise segv. Sometimes blank though? 
	//G4cout<<"process is "<<hitProcessName<<G4endl;
	hitEdeposited = aStep-> GetTotalEnergyDeposit();
	hitDetectTime = 0;
	hitCopyNum = thePhysical->GetCopyNo();
	hitPhysical = thePhysical->GetName();
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
  
  G4VVisManager* pVVisManager= G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(hitPos);
    circle.SetScreenSize(5.0);
    circle.SetFillStyle(G4Circle::filled);

    G4Color color, goodColor(1.,0.,0.), badColor(0.,0.,1.);
    
     //G4Colour  white   ()              ;  // white
     G4Colour  white   (1.0, 1.0, 1.0) ;  // white
     G4Colour  grey    (0.5, 0.5, 0.5) ;  // grey
     G4Colour  black   (0.0, 0.0, 0.0) ;  // black
     G4Colour  red     (1.0, 0.0, 0.0) ;  // red
     G4Colour  green   (0.0, 1.0, 0.0) ;  // green
     G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
     G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
     G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta 
     G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow
    

    color=green;

    G4VisAttributes attribs(color);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
  
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
  G4cout << "touchable = " << hitPhysical << G4endl;

}



