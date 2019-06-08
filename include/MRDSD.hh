// ====================================================================
//   MRDSD.hh
//
//   
// ====================================================================
#ifndef MRD_SD_H
#define MRD_SD_H
 
#include "G4VSensitiveDetector.hh"
//#include "G4EmCalculator.hh"
#include "MRDHit.hh"
//#include "MRDResponse.hh"
#include "MRDDetectorConstruction.hh"  

class G4HCofThisEvent;
class G4Step;				// ?? why do we need these?
class G4TouchableHistory;		// ??

class MRDSD : public G4VSensitiveDetector {

private:

  MRDHitsCollection* hitsCollection;
  MRDHit* newHit;
  G4int hcid;
  //WCLiteDetectorConstruction* detector;
  
  //MRDResponse* mrdresp;
  //SBsimMRDDB* mrddb;

  G4int         nhit;

  
  /*
  G4double      dedxtimesstep[NBUF_MRD];
  G4double      steplength[NBUF_MRD];
  G4double      ADCcount[NBUF_MRD];

  G4EmCalculator emcal;
  */
  
public:

  MRDSD(const G4String& name);
  virtual ~MRDSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};
#endif
