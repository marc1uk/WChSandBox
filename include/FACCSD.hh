// ====================================================================
//   FACCSD.hh
//
//   
// ====================================================================
#ifndef FACC_SD_H
#define FACC_SD_H
 
#include "G4VSensitiveDetector.hh"
//#include "G4EmCalculator.hh"
#include "FACCHit.hh"
//#include "MRDResponse.hh"
#include "MRDDetectorConstruction.hh"  

class G4HCofThisEvent;
class G4Step;				// ?? why do we need these?
class G4TouchableHistory;		// ??

class FACCSD : public G4VSensitiveDetector {

private:

  FACCHitsCollection* hitsCollection;
  FACCHit* newHit;
  G4int hcid;
  G4int nhit;

  /*
  G4double      dedxtimesstep[NBUF_MRD];
  G4double      steplength[NBUF_MRD];
  G4double      ADCcount[NBUF_MRD];

  G4EmCalculator emcal;
  */
  
public:

  FACCSD(const G4String& name);
  virtual ~FACCSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};
#endif
