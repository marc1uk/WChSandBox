// ====================================================================
//   NCVSD.cc
//
//   26/11/15 M. O'Flaherty (based on SBsimMRDSD by 2006/03/03 K. Hiraide)
// ====================================================================
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4RunManager.hh"

#include "NCVSD.hh"
#include "NCVHit.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"

//#include <CLHEP/Random/Randomize.h>

//////////////////////////////////////////////////
NCVSD::NCVSD(const G4String& name) : G4VSensitiveDetector(name) {

  collectionName.insert("NCVHitsCollection");
  G4cout<<"Constructing NCV Sensitive Detector."<<G4endl;
  
  hcid = -1;

  //G4RunManager *runMgr = G4RunManager::GetRunManager();
  //detector = (WCLiteDetectorConstruction*)runMgr->GetUserDetectorConstruction();
  //NCVresp = new SBsimNCVResponse(detector);
  //NCVdb   = detector->GetNCVDB();
}

///////////////////////////////
NCVSD::~NCVSD(){}

//////////////////////////////////////////////////
void NCVSD::Initialize(G4HCofThisEvent* HCTE){

  // create hit collection(s)
  G4cout<<"Initializing NCV hit collection"<<G4endl;
  hitsCollection = new NCVHitsCollection(SensitiveDetectorName, collectionName[0]); 
  
  // push H.C. to "Hit Collection of This Event"
  
  hcid= GetCollectionID(0); 
  //if(hcid < 0) { hcid = G4SDManager::GetSDMpointer()->GetCollectionID(0); }
  HCTE-> AddHitsCollection(hcid, hitsCollection);

  // initialization
  nhit=0;

}

//////////////////////////////////////////////////////////////
G4bool NCVSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) { // weird syntax; type declared but no variable??

    G4Track* aTrack = aStep->GetTrack();
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    // don't make hits from optical absorption...
    if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){return true;}
    // seem to have StepStatus = fGeomBoundary, ProcessDefinedStep = Transportation. not absorption..??
    	
    // don't make hits if no energy deposited...
    G4double eDeposited = aStep->GetTotalEnergyDeposit();
    if (eDeposited == 0) return true;

    ++nhit;    
    newHit= new NCVHit(aStep);
    newHit->SetHitID(nhit);      
    hitsCollection->insert(newHit);
    //G4cout<<G4endl<<"An NCV hit!"<<G4endl;
    //newHit->Print();

  // Scintillator response (quenching, attenuation)
  // This calls in code from SBSimNCVResponse which applies attenuation and Birk's constant to modify the deposited enegy. 
  // Probably good to implement, but depends on the geometry.
  //NCVresp->ApplyScintiResponse(&edepbuf,aTrack,preStepPoint);	//pass in ref to edepbuf (eDeposited from track) which is modified.

  /*
  if(particle->GetPDGCharge()!=0 && kineticE>0) {
    G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/mm);
    steplength[id] += (aStep->GetStepLength())/mm;
    dedxtimesstep[id] += dedx*((aStep->GetStepLength())/mm);
  }
  */
  return true;
}

//////////////////////////////////////////////////////
void NCVSD::EndOfEvent(G4HCofThisEvent* /*HCTE*/)
//////////////////////////////////////////////////////
{
	//SciBooNE code adds noise at this point. Good idea?
	G4cout << "**NCVSD end of event action**" << G4endl;
	G4cout << "Number of NCV hits: " << nhit << G4endl;
	/*for (G4int i=0;i<5;i++){
		NCVHit* ahit = (NCVHit*)hitsCollection->GetHit(i);
		ahit->Print();
	}*/
	//hitsCollection->PrintAllHits();
}

/////////////////////////////
void NCVSD::DrawAll() {
  hitsCollection-> DrawAllHits();
}

//////////////////////////////
void NCVSD::PrintAll() {
  if(nhit>0){
    G4cout << "##### NCV Hit #####" << G4endl;
    hitsCollection-> PrintAllHits();
    G4cout << G4endl;
  }
}
