#ifndef WCLiteTrackInformation_h
#define WCLiteTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include "WCLiteTrajectory.hh"

class WCLiteTrackInformation : public G4VUserTrackInformation 
{
  public:
    WCLiteTrackInformation();
    WCLiteTrackInformation(const G4Track* aTrack);
    WCLiteTrackInformation(const WCLiteTrackInformation* aTrackInfo);
    virtual ~WCLiteTrackInformation();
   
    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);
    inline int operator ==(const WCLiteTrackInformation& right) const
    {return (this==&right);}

    void Print() const;

  private:
  	G4bool scatteredn;
  	G4bool primaryn;
  	WCLiteTrajectory* parenttraj;
  	G4int trackID;
  	
    G4int                 originalTrackID;
//    G4ParticleDefinition* particleDefinition;
    G4ThreeVector         originalPosition;
//    G4ThreeVector         originalMomentum;
//    G4double              originalEnergy;
    G4double              originalTime;


  public:
  	inline G4bool GetIsScatteredN() const {return scatteredn;}
  	inline void SetIsScatteredN(G4bool isscattdn) { scatteredn = isscattdn;}
  	
  	inline G4bool GetIsPrimaryN() const {return primaryn;}
  	inline void SetIsPrimaryN(G4bool isprimary) { primaryn = isprimary;}
  	  	
  	inline WCLiteTrajectory* GetParentTrajectory() const {return parenttraj;}
  	inline void SetParentTrajectory(WCLiteTrajectory* intraj) { parenttraj = intraj;}
  	
  	inline G4int GetTrackID() const {return trackID;}
  	inline void SetTrackID(G4int trackin){trackID = trackin;}
  	
    inline G4int GetOriginalTrackID() const {return originalTrackID;}
//    inline G4ParticleDefinition* GetOriginalParticle() const {return particleDefinition;}
    inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
//    inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
//    inline G4double GetOriginalEnergy() const {return originalEnergy;}
    inline G4double GetOriginalTime() const {return originalTime;}
};

extern G4Allocator<WCLiteTrackInformation> aTrackInformationAllocator;

inline void* WCLiteTrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void WCLiteTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((WCLiteTrackInformation*)aTrackInfo);}

#endif
