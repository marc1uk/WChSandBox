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
// $Id: faccPMThit.hh $
//
//
//
#ifndef faccPMThit_h
#define faccPMThit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"

//#include "tls.hh"

class G4VTouchable;

class faccPMThit : public G4VHit
{
  public:
 
    //faccPMThit();
    faccPMThit(const G4Step* aStep);
    virtual ~faccPMThit();
    faccPMThit(const faccPMThit &right);

    const faccPMThit& operator=(const faccPMThit &right);
    G4int operator==(const faccPMThit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();
    
    inline void SetHitID (G4int fHitIDin) {fHitID = fHitIDin;}
    inline G4int GetHitID() { return fHitID; }
  
    inline void SetHitTrackID (G4int fTrackIDin) {fTrackID = fTrackIDin;}
    inline G4int GetHitTrackID() { return fTrackID; }
    
    inline void SetHitParentID (G4int fParentIDin) {fParentID = fParentIDin;}
    inline G4int GetHitParentID() { return fParentID; }
    
    inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
    inline G4int GetPMTNumber() { return fPmtNumber; }
    
//    inline void IncPhotonCount(){fPhotons++;}
//    inline G4int GetPhotonCount(){return fPhotons;}
    
    inline void SetHitPos(const G4ThreeVector posin){fPos=posin;}
    inline G4ThreeVector GetHitPos(){return fPos;}
    
    inline void SetHitTime(G4double fTimein){fTime=fTimein;}
    inline G4double GetHitTime(){return fTime;}
    
    inline void SetPMTPhysVol(G4VPhysicalVolume* physVol){this->fPhysVol=physVol;}
    inline G4VPhysicalVolume* GetPMTPhysVol(){return fPhysVol;}
    
    inline void SetDrawit(G4bool b){fDrawit=b;}
    inline G4bool GetDrawit(){return fDrawit;}
    
    inline void SetCreationProcess(G4String fCreationProcessin){fCreationProcess=fCreationProcessin;}
    inline G4String GetCreationProcess(){return fCreationProcess;}
    
    inline void SetHitWavelength(G4double fhitWavelengthin){fhitWavelength=fhitWavelengthin;}
    inline G4double GetHitWavelength(){return fhitWavelength;}


  private:

    G4int fHitID;
    G4int fTrackID;
    G4int fParentID;
    G4int fPmtNumber;
//    G4int fPhotons;
    G4ThreeVector fPos;
    G4double fTime;
    G4VPhysicalVolume* fPhysVol;
    G4bool fDrawit;
    G4String fCreationProcess;
    G4double fhitWavelength;
    
};

typedef G4THitsCollection<faccPMThit> faccPMThitsCollection;

extern G4Allocator<faccPMThit> faccPMThitAllocator;

inline void* faccPMThit::operator new(size_t){
    return (void *) faccPMThitAllocator.MallocSingle();
}

inline void faccPMThit::operator delete(void *aHit){
  faccPMThitAllocator.FreeSingle((faccPMThit*) aHit);
}


#endif
