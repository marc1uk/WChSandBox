#ifndef WCLitePrimaryGeneratorAction_h
#define WCLitePrimaryGeneratorAction_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TTree.h"
#include "TRandom3.h"

#include <fstream>

class WCLiteDetectorConstruction;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class WCLitePrimaryGeneratorMessenger;

class WCLitePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  WCLitePrimaryGeneratorAction(WCLiteDetectorConstruction*);
  ~WCLitePrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);

  // Normal gun setting calls these functions to fill jhfNtuple and Root tree
  void SetVtx(G4ThreeVector i)     { vtx = i; };
  void SetBeamEnergy(G4double i)   { beamenergy = i; };
  void SetBeamDir(G4ThreeVector i) { beamdir = i; };
  void SetBeamPDG(G4int i)         { beampdg = i; };

  // These go with jhfNtuple
  G4int GetVecRecNumber(){return vecRecNumber;}
  G4int GetMode() {return mode;};
  G4int GetVtxVol() {return vtxvol;};
  G4ThreeVector GetVtx() {return vtx;}
  G4int GetNpar() {return npar;};
  G4int GetBeamPDG() {return beampdg;};
  G4double GetBeamEnergy() {return beamenergy;};
  G4ThreeVector GetBeamDir() {return beamdir;};
  G4int GetTargetPDG() {return targetpdg;};
  G4double GetTargetEnergy() {return targetenergy;};
  G4ThreeVector GetTargetDir() {return targetdir;};

  // older ...
  G4double GetNuEnergy() {return nuEnergy;};
  G4double GetEnergy() {return energy;};
  G4double GetXPos() {return xPos;};
  G4double GetYPos() {return yPos;};
  G4double GetZPos() {return zPos;};
  G4double GetXDir() {return xDir;};
  G4double GetYDir() {return yDir;};
  G4double GetZDir() {return zDir;};

private:
  WCLiteDetectorConstruction*      myDetector;
  G4ParticleGun*                  particleGun;
  G4GeneralParticleSource*        MyGPS;  //T. Akiri: GPS to run Laser
  WCLitePrimaryGeneratorMessenger* messenger;

  TTree *tcardfile;
  G4int mevt,ntrks,nneutrons,mmode,mbeam_id;
  G4double mvtxx,mvtxy,mvtxz,mE,mbeam_px,mbeam_py,mbeam_pz;
  G4int *mpid;
  G4double *mpx;
  G4double *mpy;
  G4double *mpz;
  G4double *mKE;

  /*
  G4int mpid[kmaxtrk];
  G4double mpx[kmaxtrk];
  G4double mpy[kmaxtrk];
  G4double mpz[kmaxtrk];
  G4double mKE[kmaxtrk];
  */

  G4int nMRDlayers;
  G4double mMRDtotEdep,mMRDtotEchdep; 
  
  G4int *mMRDhitlayer;
  G4int *mMRDhitorientation;
  G4double *mMRDhitEdep;
  G4double *mMRDhitEchdep;
  G4double *mMRDhitx;
  G4double *mMRDhity;
  G4double *mMRDhitz;
  
  G4int numlayers,mrdstrp_orient,msum;
  G4double edep,echdep;
  G4bool   bobfile;

  // Variables set by the messenger
  G4bool   useMulineEvt;
  G4bool   useNormalEvt;
  G4bool   useLaserEvt;  //T. Akiri: Laser flag
  std::fstream inputFile;
  G4String vectorFileName;
  G4bool   GenerateVertexInRock;

  // These go with jhfNtuple
  G4int mode;
  G4int vtxvol;
  G4ThreeVector vtx;
  G4int npar;
  G4int beampdg, targetpdg;
  G4ThreeVector beamdir, targetdir;
  G4double beamenergy, targetenergy;
  G4int vecRecNumber;

  G4double nneut,np;

  G4double nuEnergy;
  G4double energy;
  G4double xPos, yPos, zPos;
  G4double xDir, yDir, zDir;

  G4int    _counterRock; 
  G4int    _counterCublic; 

  TRandom3* mR;

public:

  inline void SetMulineEvtGenerator(G4bool choice) { useMulineEvt = choice; }
  inline G4bool IsUsingMulineEvtGenerator() { return useMulineEvt; }

  inline void SetNormalEvtGenerator(G4bool choice) { useNormalEvt = choice; }
  inline G4bool IsUsingNormalEvtGenerator()  { return useNormalEvt; }

  //T. Akiri: Addition of function for the laser flag
  inline void SetLaserEvtGenerator(G4bool choice) { useLaserEvt = choice; }
  inline G4bool IsUsingLaserEvtGenerator()  { return useLaserEvt; }

  inline void OpenVectorFile(G4String fileName) 
  {
    if ( inputFile.is_open() ) 
      inputFile.close();

    vectorFileName = fileName;
    inputFile.open(vectorFileName, std::fstream::in);
  }
  inline G4bool IsGeneratingVertexInRock() { return GenerateVertexInRock; }
  inline void SetGenerateVertexInRock(G4bool choice) { GenerateVertexInRock = choice; }

};

#endif


