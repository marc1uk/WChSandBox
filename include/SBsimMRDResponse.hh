
// ====================================================================
//   SBsimMRDResponse.hh
//
//   2006/09/04 K. Hiraide
// ====================================================================
#ifndef SBSIM_MRD_RESPONSE
#define SBSIM_MRD_RESPONSE

#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
#include "G4Track.hh"

//#include "SBsimDetectorConstruction.hh"
#include "WCLiteDetectorConstruction.hh"
#include "SBsimMRDDB.hh" //added

class SBsimMRDResponse {
public:
  //SBsimMRDResponse(SBsimDetectorConstruction* detector);
  SBsimMRDResponse(WCLiteDetectorConstruction* detector);
  ~SBsimMRDResponse();

  void ApplyScintiResponse(G4double* edeposit, G4Track* aTrack, G4StepPoint* preStepPoint);

private:
  //SciBarCard*    sbcard;
  //SBsimSciBarDB* sbdb;
  SBsimMRDDB* mrddb;
  G4EmCalculator emcal;

  void BirksSaturation(G4double* edeposit, G4Track* aTrack);
  void AttenuationLength(G4double* edeposit, G4Track* aTrack, G4StepPoint* preStepPoint);

};

#endif
