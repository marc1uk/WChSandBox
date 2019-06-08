class WCLiteTrajectory;

#ifndef WCLiteTrajectory_h
#define WCLiteTrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                  // Include from 'system'
#include "G4ios.hh"                  // Include from 'system'
#include <vector>                    // G4RWTValOrderedVector
#include "globals.hh"                // Include from 'global'
#include "G4ParticleDefinition.hh"   // Include from 'particle+matter'
#include "G4TrajectoryPoint.hh"      // Include from 'tracking'
#include "G4Track.hh"                // 
#include "G4Step.hh"                 // 
#include "WCLiteTrackInformation.hh" // 
#include <stdexcept>                 // std::out_of_range

class G4Polyline;                    // Forward declaration.

typedef std::vector<G4VTrajectoryPoint*>  TrajectoryPointContainer;
///////////////////
class WCLiteTrajectory : public G4VTrajectory
///////////////////
{

//--------
public: // with description
//--------

// Constructor/Destrcutor

   WCLiteTrajectory();

   WCLiteTrajectory(const G4Track* aTrack);
   WCLiteTrajectory(WCLiteTrajectory &);
   virtual ~WCLiteTrajectory();

// Operators
   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const WCLiteTrajectory& right) const
   {return (this==&right);} 

// Get/Set functions 
   inline G4int GetTrackID() const
   { return fTrackID; }
   inline G4int GetParentID() const
   { return fParentID; }
   inline G4String GetParticleName() const
   { return ParticleName; }
   inline G4double GetCharge() const
   { return PDGCharge; }
   inline G4int GetPDGEncoding() const
   { return PDGEncoding; }
   inline G4ThreeVector GetInitialMomentum() const
   { return initialMomentum; }
   inline G4String GetCreatorProcessName() const
   { return creatorProcess; }
   inline G4double GetGlobalTime() const
   { return globalTime; }
   inline G4bool GetSaveFlag() const 
   { return SaveIt; }
   
   inline void SetSaveFlag(G4bool value) { SaveIt = value; }

// New function we have added
   inline G4ThreeVector GetStoppingPoint() const
   { return stoppingPoint; }
   inline G4VPhysicalVolume* GetStoppingVolume() const
   { return stoppingVolume;}
   inline void SetStoppingPoint(G4ThreeVector& currentPosition) 
   { stoppingPoint = currentPosition;}
   inline void SetStoppingVolume(G4VPhysicalVolume* currentVolume)
   { stoppingVolume = currentVolume;}


// Other member functions
   virtual void ShowTrajectory(std::ostream& os=G4cout) const;
   using G4VTrajectory::DrawTrajectory;	//import all overloads
   virtual void DrawTrajectory(G4int i_mode=0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual int GetPointEntries() const { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const {
		 try{ return positionRecord->at(i); } 
		 catch (const std::out_of_range& oor) {
		  	std::cerr << "Out of Range error: " << oor.what() << '\n';
		 }
   }
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;
   
   void CopyInfo(WCLiteTrackInformation* infoin);	 //************ copy stuff from track information
   
	 inline G4int GetSameAsParentTrackID() const {return sameasparenttrackid;}
	 inline void SetSameAsParentTrackID(G4int sameasparenttrackidin){sameasparenttrackid = sameasparenttrackidin;}
  
	 inline G4int GetNumSecs() const {return numSecondaries;}
	 inline void AddNumSecs(G4int numSecondariesin){numSecondaries += numSecondariesin;}
	 
	 inline G4int GetGdOriginalTrackID() const { return GdOriginalTrackID; }
	 inline void SetGdOriginalTrackID(G4int GdOriginalTrackIDin){GdOriginalTrackID = GdOriginalTrackIDin;}
	 
	 inline G4int GetPassThruNCV() const {return passesthroughNCV;}
	 inline void SetPassThruNCV(G4int passesthroughNCVin){passesthroughNCV = passesthroughNCVin;} 
	 
	 inline G4ThreeVector GetNCVentryPos() const {return ncvEntryPos;}
	 inline void SetNCVentryPos(G4ThreeVector ncvEntryPosin){ncvEntryPos = ncvEntryPosin;} 
	 
	 inline G4double GetNCVentryTime() const { return ncvEntryTime; }
	 inline void SetNCVentryTime(G4double ncvEntryTimein){ncvEntryTime = ncvEntryTimein;}
	 
	 inline G4ThreeVector GetNCVexitPos() const {return ncvExitPos;}
	 inline void SetNCVexitPos(G4ThreeVector ncvExitPosin){ncvExitPos = ncvExitPosin;} 
	 
	 inline G4double GetNCVexitTime() const { return ncvExitTime; }
	 inline void SetNCVexitTime(G4double ncvExitTimein){ncvExitTime = ncvExitTimein;}
	 //inline G4double GetTotalEdep() const {return totalEdep;}
	 //inline void AddTotalEdep(G4int totalEdepin){totalEdep += totalEdepin;}  
	 
	 void Print() const;
   
   // ************
//---------
   private:
//---------

  TrajectoryPointContainer* positionRecord;
  G4int                     fTrackID;
  G4int                     fParentID;
  G4int                     PDGEncoding;
  G4double                  PDGCharge;
  G4String                  ParticleName;
  G4ThreeVector             initialMomentum;
  // These are new variables
  G4ThreeVector             stoppingPoint;
  G4VPhysicalVolume         *stoppingVolume;
  // ----------------
  // neutron counting
  // ----------------
  G4int                     sameasparenttrackid;
  G4int                     numSecondaries;
  //G4double                totalEdep;						// not sure the validity of this
  //G4String                exitProcess;
  // M Fechner : new saving mechanism
  G4bool                    SaveIt;
  G4String                  creatorProcess;
  G4double                  globalTime;
  // ----------------------
  // NCV efficiency testing
  // ----------------------
  G4int                     GdOriginalTrackID;
  G4int                     passesthroughNCV;
  G4ThreeVector             ncvEntryPos;
  G4double                  ncvEntryTime;
  G4ThreeVector             ncvExitPos;
  G4double                  ncvExitTime;
};

/***            TEMP  : M FECHNER ***********
** modification by Chris Walter that works for geant4 >= 4.6.2p01
** does not compile with 4.6.1
#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;
#else
  extern G4DLLIMPORT G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;
#endif
*/

extern G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;

inline void* WCLiteTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void WCLiteTrajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((WCLiteTrajectory*)aTrajectory);
}

#endif

