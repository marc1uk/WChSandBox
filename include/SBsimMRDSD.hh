// ====================================================================
//   SBsimMRDSD.hh
//
//   2006/03/03 K. Hiraide
// ====================================================================
#ifndef SBSIM_MRD_SD_H
#define SBSIM_MRD_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "G4EmCalculator.hh"
#include "SBsimMRDHit.hh"
#include "SBsimMRDResponse.hh"

#include "SBsimMRDDB.hh" //added

#define NBUF_MRD 512

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class SBsimMRDSD : public G4VSensitiveDetector {
private:
  SBsimMRDHitsCollection* hitsCollection;
  SBsimMRDResponse* mrdresp;
  SBsimMRDHit* newHit;
  //SBsimDetectorConstruction* detector;
  WCLiteDetectorConstruction* detector;
  SBsimMRDDB* mrddb;

  G4int         nhit;
  G4int         id;
  G4int         initflag[NBUF_MRD];

  //  G4double      edepbuf[NBUF_MRD];
  G4double      edepbuf;
  G4double      edeposit[NBUF_MRD];
  G4double      edeposit_resp[NBUF_MRD];
  G4int         itrack[NBUF_MRD];
  G4int         ipart[NBUF_MRD];
  G4String      parname[NBUF_MRD];
  G4double      hittime[NBUF_MRD];
  G4double      detecttime[NBUF_MRD];
  G4ThreeVector hitpos[NBUF_MRD];
  
  G4double      dedxtimesstep[NBUF_MRD];
  G4double      steplength[NBUF_MRD];
  G4double      ADCcount[NBUF_MRD];

  G4EmCalculator emcal;

public:
  SBsimMRDSD(const G4String& name);
  virtual ~SBsimMRDSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif
