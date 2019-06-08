// ====================================================================
//   NCVHit.hh
//
//   26/11/15
// ====================================================================
#ifndef NCVHit_h
#define NCVHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"

class NCVHit : public G4VHit
{
  public:

      NCVHit();
      NCVHit(G4Step* aStep);
      ~NCVHit();
      //NCVHit(const NCVHit &right);
      //const NCVHit& operator=(const NCVHit &right); //same as implicit definition
      int operator==(const NCVHit &right) const;

      inline void * operator new(size_t);
      inline void operator delete(void *aHit);

      virtual void Draw();
      virtual void Print();

  private:
      G4int hitID;
      G4int hitTrackID;
      G4ThreeVector hitPos;
      G4double hitTime;
      G4String hitParticleName;
      G4int hitParticleID; 	// from GetPDGEncoding
      G4String hitProcessName;
      G4double hitEdeposited;
      G4double hitDetectTime;
      G4int hitCopyNum;
      G4String hitPhysical;
      //G4int hitPanelNum;
      //G4int hitPaddleNum;

  public:
  void SetHitID (G4int hitIDin) {hitID = hitIDin;}
  G4int GetHitID() { return hitID; }
  
  void SetHitTrackID (G4int trackIDin) {hitTrackID = trackIDin;}
  G4int GetHitTrackID() { return hitTrackID; }
  
  void SetHitPos (G4ThreeVector hitPosin) { hitPos = hitPosin; }
  G4ThreeVector GetHitPos() { return hitPos; }
  
  void SetHitTime (G4double hitTimein) { hitTime = hitTimein; }
  G4double GetHitTime() { return hitTime; }
  
  void SetHitParticleName (G4String partNamein) { hitParticleName = partNamein; }
  G4String GetHitParticleName() { return hitParticleName; }
  
  void SetHitParticleID (G4int hitParticleIDin) { hitParticleID = hitParticleIDin; }
  G4int GetHitParticleID() { return hitParticleID; }
  
  void SetHitProcessName (G4String processNamein) { hitProcessName = processNamein; }
  G4String GetHitProcessName() { return hitProcessName; }  
  
  void SetHitEdeposit (G4double hitEdepositedin) {hitEdeposited = hitEdepositedin;}
  G4double GetHitEdeposit() { return hitEdeposited; }
    
  void SetHitDetectTime (G4double hitDetectTimein) { hitDetectTime = hitDetectTimein; }
  G4double GetHitDetectTime() { return hitDetectTime; }
  
  void SetHitCopyNum (G4int hitCopyNumin) {hitCopyNum = hitCopyNumin;}
  G4int GetHitCopyNum() { return hitCopyNum; }
  
  void SetHitPhysical (G4String hitPhysicalin) {hitPhysical = hitPhysicalin;}
  G4String GetHitPhysical() { return hitPhysical; }
  
/*  void SetHitPanelNum (G4int hitPanelNumin) {hitPanelNum = hitPanelNumin;}
  G4int GetHitPanelNum() { return hitPanelNum; }
  
  void SetHitPaddleNum (G4int hitPaddleNumin) {hitPaddleNum = hitPaddleNumin;}
  G4int GetHitPaddleNum() { return hitPaddleNum; }
*/
};

typedef G4THitsCollection<NCVHit> NCVHitsCollection;

extern G4Allocator<NCVHit> NCVHitAllocator;

inline void* NCVHit::operator new(size_t)
{
  return (void *) NCVHitAllocator.MallocSingle();
}

inline void NCVHit::operator delete(void *aHit)
{
  NCVHitAllocator.FreeSingle((NCVHit*) aHit);
}

inline G4int NCVHit::operator==
       (const NCVHit& right) const 
{
   return (this==&right) ? 1 : 0; 
}

/*inline const SBsimNCVHit& SBsimNCVHit::operator= (const SBsimNCVHit& right) {
// Think this is the default implicitly defined = operator so needn't be defined explicitly
	hitID = right.hitID
	hitTrackID = right.hitTrackID
	hitPos = right.hitPos
	hitTime = right.hitTime
	hitParticleName = right.hitParticleName
	hitParticleID = right.hitParticleID
	hitProcessName = right.hitProcessName
	eDeposited = right.eDeposited
	hitDetectTime = right.hitDetectTime
	hitPanelNum = right.hitPanelNum
	hitPaddleNum = right.hitPaddleNum
  return *this;
}
*/

#endif

