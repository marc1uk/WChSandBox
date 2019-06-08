#include "WCLiteTrackInformation.hh"
#include "WCLiteTrajectory.hh"
#include "G4ios.hh"

G4Allocator<WCLiteTrackInformation> aTrackInformationAllocator;

WCLiteTrackInformation::WCLiteTrackInformation() {
	scatteredn = false;
	primaryn = true;
	parenttraj = 0; 	// no parent trajectory
	trackID = 0;		// set in function
	 
    originalTrackID = 0;
//    particleDefinition = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
//    originalMomentum = G4ThreeVector(0.,0.,0.);
//    originalEnergy = 0.;
    originalTime = 0.;
}
WCLiteTrackInformation::WCLiteTrackInformation(const G4Track* aTrack) {
	scatteredn = true;	
	primaryn = false;	
	parenttraj = 0;	// set in function
	trackID = 0;	// set in function

    originalTrackID = aTrack->GetTrackID();
//    particleDefinition = aTrack->GetDefinition();
    originalPosition = aTrack->GetPosition();
//    originalMomentum = aTrack->GetMomentum();
//    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();
}
//WCLiteTrackInformation::WCLiteTrackInformation(const WCLiteTrackInformation* aTrackInfo) { 
//	scatteredn = aTrackInfo->scatteredn;
//    originalTrackID = aTrackInfo->originalTrackID;
//    particleDefinition = aTrackInfo->particleDefinition;
//    originalPosition = aTrackInfo->originalPosition;
//    originalMomentum = aTrackInfo->originalMomentum;
//    originalEnergy = aTrackInfo->originalEnergy;
//    originalTime = aTrackInfo->originalTime;
//}

WCLiteTrackInformation::~WCLiteTrackInformation(){;}

void WCLiteTrackInformation::Print() const
{
    if(scatteredn){ G4cout << "Particle is a scattered neutron." << G4endl;}
    if(primaryn) { G4cout << "Particle is a primary neutron." << G4endl;}
}


