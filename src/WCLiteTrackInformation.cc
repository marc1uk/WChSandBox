#include "WCLiteTrackInformation.hh"
#include "WCLiteTrajectory.hh"
#include "G4ios.hh"

G4Allocator<WCLiteTrackInformation> aTrackInformationAllocator;
    
WCLiteTrackInformation::WCLiteTrackInformation() {  
	//G4cout << "construction() of: " << this << G4endl;
	// -----------------
	// neutron counting
	// -----------------
	scatteredn=false;
	primaryn=true;
	parenttraj=0;		// no parent trajectory
	trackID =0;			// set in function
	trajectoryID=0;
	sameasparenttrackid=0;
	numSecondaries=0;
	originalTrackID=0;
	originalPosition=G4ThreeVector(0.,0.,0.);
	originalTime=0.;
//	particleDefinition=0;
//	originalMomentum=G4ThreeVector(0.,0.,0.);
//	originalEnergy=0.;
// -----------------------
// NCV efficiency checking
// -----------------------
	GdOriginalTrackID=-1;
	passesthroughNCV=0;
	ncvEntryPos=G4ThreeVector();
	ncvEntryTime=0.;
	ncvExitPos=G4ThreeVector();
	ncvExitTime=0.;
}

WCLiteTrackInformation::WCLiteTrackInformation(const G4Track* aTrack) {
	//G4cout << "construction(G4Track*) of: " << this << G4endl;
	// -----------------
	// neutron counting
	// -----------------
	scatteredn = true;
	primaryn = false;
	parenttraj = 0;	// set in function
	trackID = 0;		// set in function
	sameasparenttrackid=0;
	trajectoryID=0;
	numSecondaries=0;
	originalTrackID=aTrack->GetTrackID();
	originalPosition=aTrack->GetPosition();
	originalTime=aTrack->GetGlobalTime();
//	particleDefinition=aTrack->GetDefinition();
//	originalMomentum=aTrack->GetMomentum();
//	originalEnergy=aTrack->GetTotalEnergy();
// -----------------------
// NCV efficiency checking
// -----------------------
	GdOriginalTrackID=-1;
	passesthroughNCV=0;
	ncvEntryPos=G4ThreeVector();
	ncvEntryTime=0.;
	ncvExitPos=G4ThreeVector();
	ncvExitTime=0.;
}

//WCLiteTrackInformation::WCLiteTrackInformation(const WCLiteTrackInformation* aTrackInfo) { 
//    scatteredn = aTrackInfo->scatteredn;
//    originalTrackID = aTrackInfo->originalTrackID;
//    particleDefinition = aTrackInfo->particleDefinition;
//    originalPosition = aTrackInfo->originalPosition;
//    originalMomentum = aTrackInfo->originalMomentum;
//    originalEnergy = aTrackInfo->originalEnergy;
//    originalTime = aTrackInfo->originalTime;
//}

WCLiteTrackInformation::~WCLiteTrackInformation(){
	//G4cout << "destruction of: " << this << G4endl;
}

void WCLiteTrackInformation::Print() const
{
    G4cout<< "scatteredn = " << scatteredn << G4endl;	
    G4cout<< "primaryn = " << primaryn << G4endl;
    G4cout<< "parenttraj = " << parenttraj << G4endl;
    G4cout<< "trackID = " << trackID << G4endl;
    G4cout<< "sameasparenttrackid = " << sameasparenttrackid << G4endl;
    G4cout<< "trajectoryID = " << trajectoryID << G4endl;
    G4cout<< "numSecondaries = " << numSecondaries << G4endl;
    G4cout<< "ncvOriginalTrackID = " << GdOriginalTrackID << G4endl;
    G4cout<< "originalTrackID = " << originalTrackID << G4endl;
    G4cout<< "originalPosition = " << originalPosition << G4endl;
    G4cout<< "originalTime = " << originalTime << G4endl;
}


