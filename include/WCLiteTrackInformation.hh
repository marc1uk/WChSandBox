#ifndef WCLiteTrackInformation_h
#define WCLiteTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
//#include "WCLiteTrajectory.hh"	// circular dependancy; can use forward declaration.
class WCLiteTrajectory;

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
    G4bool                scatteredn;			// main reason for creating user info class: keep track of primary...
    G4bool                primaryn;				// vs inelastically scattered neutrons for a proper count of total neutron prodn.
    WCLiteTrajectory*     parenttraj;
    G4int                 trajectoryID;
    G4int                 trackID;
    G4int                 sameasparenttrackid;			// for inelastic neutron scattering daughters, call one the same as parent.
    G4int                 originalTrackID;					// stores the trackID of the ancestor track without a parent
    G4int                 numSecondaries;						// needed to pass on ancestry info
    G4int                 GdOriginalTrackID;				// stores the trackID of the ancestor Gd capture if applicable
    G4int                 passesthroughNCV;					// record if particle has passed through NCV for efficiency testing
    G4ThreeVector         ncvEntryPos;							// position the particle entered the NCV (if applicable)
    G4double              ncvEntryTime;
    G4ThreeVector         ncvExitPos;								// position the particle left the NCV (if applicable)
    G4double              ncvExitTime;
    G4ThreeVector         originalPosition;
    G4double              originalTime;
//  G4double              totalEdep;								// not sure the validity of this
//  G4String              exitProcess;
//  G4ParticleDefinition* particleDefinition;
//  G4ThreeVector         originalMomentum;
//  G4double              originalEnergy;

  public:
    inline G4bool GetIsScatteredN() const {return scatteredn;}
    inline void SetIsScatteredN(G4bool isscattdn) { scatteredn = isscattdn;}
    
    inline G4bool GetIsPrimaryN() const {return primaryn;}
    inline void SetIsPrimaryN(G4bool isprimary) { primaryn = isprimary;}
    
    inline WCLiteTrajectory* GetParentTrajectory() const {return parenttraj;}
    inline void SetParentTrajectory(WCLiteTrajectory* intraj) { parenttraj = intraj;}
    
    inline G4int GetTrackID() const {return trackID;}
    inline void SetTrackID(G4int trackin){trackID = trackin;}
    
    inline G4int GetSameAsParentTrackID() const {return sameasparenttrackid;}
    inline void SetSameAsParentTrackID(G4int sameasparenttrackidin){sameasparenttrackid = sameasparenttrackidin;}
    
    inline G4int GetNumSecs() const {return numSecondaries;}
    inline void AddNumSecs(G4int numSecondariesin){numSecondaries += numSecondariesin;}
    
    inline G4int GetGdOriginalTrackID() const {return GdOriginalTrackID;}
    inline void SetGdOriginalTrackID(G4int GdOriginalTrackIDin){GdOriginalTrackID = GdOriginalTrackIDin;}  
    
    inline G4int GetPassThruNCV() const {return passesthroughNCV;}
    inline void SetPassThruNCV(G4int passesthroughNCVin){passesthroughNCV = passesthroughNCVin;} 
    
    inline G4ThreeVector GetNCVentryPos() const {return ncvEntryPos;}
    inline void SetNCVentryPos(G4ThreeVector ncvEntryPosin){ncvEntryPos = ncvEntryPosin;}
    
    inline G4double GetNCVentryTime() const {return ncvEntryTime;}
    inline void SetNCVentryTime(G4double ncvEntryTimein){ncvEntryTime = ncvEntryTimein;}
    
    inline G4ThreeVector GetNCVexitPos() const {return ncvExitPos;}
    inline void SetNCVexitPos(G4ThreeVector ncvExitPosin){ncvExitPos = ncvExitPosin;}
    
    inline G4double GetNCVexitTime() const {return ncvExitTime;}
    inline void SetNCVexitTime(G4double ncvExitTimein){ncvExitTime = ncvExitTimein;}
    
    inline G4double GetOriginalTime() const {return originalTime;}
    inline G4int GetOriginalTrackID() const {return originalTrackID;}
    inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
    
/*
    inline G4double GetTotalEdep() const {return totalEdep;}
    inline void AddTotalEdep(G4int totalEdepin){totalEdep += totalEdepin;}    
    
    inline G4String GetExitProcess() const {return exitProcess;}
    inline void SetExitProcess(G4String exitProcessin){exitProcess = exitProcessin;}    
    
    inline G4ParticleDefinition* GetOriginalParticle() const {return particleDefinition;}
    inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
    inline G4double GetOriginalEnergy() const {return originalEnergy;}
*/

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
