#ifndef WCSIMRECODIGIT_HH
#define WCSIMRECODIGIT_HH

#include "TObject.h"

class WCSimRecoDigit : public TObject {

 public:

  WCSimRecoDigit( Int_t region,
                  Double_t x, Double_t y, Double_t z,
                  Double_t rawT, Double_t rawQ, Int_t rawPEtube,
                  Double_t calT, Double_t calQ );

  ~WCSimRecoDigit();

  Int_t GetRegion(){ return fRegion; }

  Double_t GetX() { return fX; }
  Double_t GetY() { return fY; }
  Double_t GetZ() { return fZ; }

  Double_t GetRawTime(){ return fRawTime; }
  Double_t GetRawQPEs(){ return fRawQPEs; }

  Double_t GetTime(){ return fCalTime; }
  Double_t GetQPEs(){ return fCalQPEs; }
  
  Int_t GetRawPEtube(){ return fRawPEtube; }

  Bool_t IsFiltered(){ return fIsFiltered; }

  void SetFilter( Bool_t pass = 1 ){ fIsFiltered = pass; }
  void ResetFilter() { SetFilter(0); }
  void PassFilter()  { SetFilter(1); }

 private:

  Int_t fRegion;

  Double_t fX;
  Double_t fY;
  Double_t fZ;

  Double_t fRawTime;
  Double_t fRawQPEs;

  Double_t fCalTime;
  Double_t fCalQPEs;
  
  Int_t fRawPEtube;

  Bool_t fIsFiltered;

  ClassDef(WCSimRecoDigit,0)

};

#endif
