//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#ifndef WCLiteTrajectoryPoint_h
#define WCLiteTrajectoryPoint_h 1

#include "G4TrajectoryPoint.hh"
#include "G4Allocator.hh"
#include "G4ios.hh" 
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"

class WCLiteTrajectoryPoint : public G4TrajectoryPoint
{
public:
  
  //Constructor/Destructor
  WCLiteTrajectoryPoint();
  WCLiteTrajectoryPoint(G4ThreeVector pos, G4ThreeVector dir, G4double theGlobaltime, 
			G4double theTE, G4double theKE, G4String theDetLocation, 
			G4String theProcess);
  WCLiteTrajectoryPoint(const WCLiteTrajectoryPoint &right);

  virtual ~WCLiteTrajectoryPoint();
  
  //Operators
  G4bool operator==(const WCLiteTrajectoryPoint& right) const;
  inline void *operator new(size_t);
  inline void operator delete(void *aTrajectoryPoint);


  //Get/Set Functions
  //  virtual const G4ThreeVector GetPosition() const = 0;
  // Get/Set functions
  inline const G4ThreeVector GetPosition() const
  { return fPosition; };
  inline const G4ThreeVector GetDirection() const
  { return fDir; };
  G4double GetGlobalTime();
  G4double GetTotalEnergy();
  G4double GetKineticEnergy();
  G4String GetProcessName();
  G4String GetDetectorLocation();

private:

  G4ThreeVector fPosition,fDir;
  G4double _globaltime;
  G4double _TotalEnergy;
  G4double _KinEnergy;
  G4String _ProcessName;
  G4String _DetectorLocation;
};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<WCLiteTrajectoryPoint> WCLiteTrajectoryPointAllocator;
#else
  extern G4DLLIMPORT G4Allocator<WCLiteTrajectoryPoint> WCLiteTrajectoryPointAllocator;
#endif

inline void* WCLiteTrajectoryPoint::operator new(size_t)
{
   void *aTrajectoryPoint;
   aTrajectoryPoint = (void *) WCLiteTrajectoryPointAllocator.MallocSingle();
   return aTrajectoryPoint;
}

inline void WCLiteTrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   WCLiteTrajectoryPointAllocator.FreeSingle((WCLiteTrajectoryPoint *) aTrajectoryPoint);
}


#endif
