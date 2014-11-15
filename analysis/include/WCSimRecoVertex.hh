#ifndef WCSIMRECOVERTEX_HH
#define WCSIMRECOVERTEX_HH

#include "TObject.h"

class WCSimRecoVertex : public TObject {

 public: 

  typedef enum EFitStatus {
   kOK  = 0x00,
   kFailSimpleVertex    = 0x01,
   kFailSimpleDirection = 0x02,
   kFailPointPosition   = 0x04,
   kFailPointDirection  = 0x08,
   kFailPointVertex     = 0x10,
   kFailExtendedVertex  = 0x20
  } FitStatus_t;

  WCSimRecoVertex();
  WCSimRecoVertex( Double_t x, Double_t y, Double_t z );
  WCSimRecoVertex( Double_t x, Double_t y, Double_t z,
                   Double_t px, Double_t py, Double_t pz );
  WCSimRecoVertex( Double_t x, Double_t y, Double_t z, Double_t t,
                   Double_t px, Double_t py, Double_t pz, 
                   Double_t fom, Int_t nsteps, Bool_t pass, Int_t status );
  WCSimRecoVertex( Double_t x, Double_t y, Double_t z, Double_t t,
                   Double_t px, Double_t py, Double_t pz, 
                   Double_t angle, Double_t length,
                   Double_t fom, Int_t nsteps, Bool_t pass, Int_t status );
  ~WCSimRecoVertex();

  void SetVertex( Double_t x, Double_t y, Double_t z, Double_t t);
  void SetVertex( Double_t x, Double_t y, Double_t z );
  void SetDirection( Double_t px, Double_t py, Double_t pz);
  void SetConeAngle( Double_t angle );
  void SetTrackLength( Double_t length );
  void SetFOM(Double_t fom, Int_t nsteps, Bool_t pass);
  void SetStatus(Int_t status);

  Double_t GetX() { return fX; }
  Double_t GetY() { return fY; }
  Double_t GetZ() { return fZ; }
  Double_t GetTime() { return fTime; }
  Bool_t FoundVertex() { return fFoundVertex; }

  Double_t GetDirX() { return fDirX; }
  Double_t GetDirY() { return fDirY; }
  Double_t GetDirZ() { return fDirZ; }
  Bool_t FoundDirection() { return fFoundDirection; }

  Double_t GetConeAngle() { return fConeAngle; }
  Double_t GetTrackLength() { return fTrackLength; }

  Double_t GetFOM() { return fFOM; }
  Int_t GetIterations(){ return fIterations; }
  Bool_t GetPass() { return fPass; }
  Int_t GetStatus(){ return fStatus; }

  void Reset();

 private:

  Double_t fX;
  Double_t fY;
  Double_t fZ;
  Double_t fTime;
  Bool_t fFoundVertex;

  Double_t fDirX;
  Double_t fDirY;
  Double_t fDirZ;
  Bool_t fFoundDirection;

  Double_t fConeAngle;
  Double_t fTrackLength;

  Double_t fFOM;
  Int_t fIterations;
  Bool_t fPass;

  Int_t fStatus;

  ClassDef(WCSimRecoVertex,0)

};

#endif







