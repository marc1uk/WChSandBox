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
  	G4bool scatteredn;
  	G4bool primaryn;
  	WCLiteTrajectory* parenttraj;
  	G4int trajectoryID;
  	G4int trackID;
  	G4int sameasparenttrackid;	// for inelastic neutron scattering daughters, call one the same as parent. 
  	G4int passesthroughMRD;
  	G4int isMRDprimary;
  	G4int numSecondaries;
  	G4int mrdOriginalTrackID;	// stores the trackID of the first ancestor track to enter the MRD
  	G4ThreeVector mrdStartPos;	// position the particle entered the MRD (if applicable)
  	G4int mrdDetected;
  	//G4double totalEdep;		// not sure the validity of this
  	//G4String exitProcess;
    G4int                 originalTrackID;	// stores the trackID of the ancestor track without a parent
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
  	
  	inline G4int GetSameAsParentTrackID() const {return sameasparenttrackid;}
  	inline void SetSameAsParentTrackID(G4int sameasparenttrackidin){sameasparenttrackid = sameasparenttrackidin;}
  	
  	inline G4int GetPassThruMRD() const {return passesthroughMRD;}
  	inline void SetPassThruMRD(G4int passesthroughMRDin){passesthroughMRD = passesthroughMRDin;} 
  	// set in trajectory class AppendStep() method, since we don't have a stepping action yet.
  	
  	inline G4int GetIsMRDprimary() const {return isMRDprimary;}
  	inline void SetIsMRDprimary(G4int isMRDprimaryin){isMRDprimary = isMRDprimaryin;} 
  	
  	inline G4int GetNumSecs() const {return numSecondaries;}
  	inline void AddNumSecs(G4int numSecondariesin){numSecondaries += numSecondariesin;}  	
  	
  	inline G4int GetmrdOriginalTrackID() const {return mrdOriginalTrackID;}
  	inline void SetmrdOriginalTrackID(G4int mrdOriginalTrackIDin){mrdOriginalTrackID = mrdOriginalTrackIDin;}  
  	
  	inline G4ThreeVector GetMRDstartPos() const {return mrdStartPos;}
	inline void SetMRDstartPos(G4ThreeVector mrdStartPosin){mrdStartPos = mrdStartPosin;} 	
	
	inline G4int GetMRDdetected() const {return mrdDetected;}
	inline void SetMRDdetected(G4int mrdDetectedin){mrdDetected = mrdDetectedin;}
  	
  	//inline G4double GetTotalEdep() const {return totalEdep;}
  	//inline void AddTotalEdep(G4int totalEdepin){totalEdep += totalEdepin;}  	
  	
//  	inline G4String GetExitProcess() const {return exitProcess;}
//  	inline void SetExitProcess(G4String exitProcessin){exitProcess = exitProcessin;}  	
  	
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
