// ====================================================================
//   NCVSD.hh
//
//   
// ====================================================================
#ifndef NCV_SD_H
#define NCV_SD_H
 
#include "G4VSensitiveDetector.hh"
//#include "G4EmCalculator.hh"
#include "NCVHit.hh"
//#include "NCVResponse.hh"
#include "MRDDetectorConstruction.hh"  

class G4HCofThisEvent;
class G4Step;								
class G4TouchableHistory;		

class NCVSD : public G4VSensitiveDetector {

private:

  NCVHitsCollection* hitsCollection;
  NCVHit* newHit;
  G4int hcid;
  //WCLiteDetectorConstruction* detector;
  
  //NCVResponse* NCVresp;
  //SBsimNCVDB* NCVdb;

  G4int         nhit;

  
  /*
  G4double      dedxtimesstep[NBUF_NCV];
  G4double      steplength[NBUF_NCV];
  G4double      ADCcount[NBUF_NCV];

  G4EmCalculator emcal;
  */
  
public:

  NCVSD(const G4String& name);
  virtual ~NCVSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};
#endif
