#ifndef WCSIMTRUECAPTURE_HH
#define WCSIMTRUECAPTURE_HH

#include "WCSimRecoDigit.hh"
#include "TObject.h"

class WCSimTrueCapture : public TObject {

 public: 

  WCSimTrueCapture(Double_t x, Double_t y, Double_t z, Double_t t, Double_t E, Int_t capnum, Int_t capnuc, Int_t pid, Int_t ngamma, Int_t nphot); 

  ~WCSimTrueCapture();

  Double_t GetX() { return fX; }
  Double_t GetY() { return fY; }
  Double_t GetZ() { return fZ; }
  Double_t GetT() { return fT; }
  Double_t GetE() { return fE; }


  Int_t GetCaptNum() { return fCnum;}
  Int_t GetCaptNucleus() {return fCnucl;}
  Int_t GetNGamma() {return fNgamma;}
  Int_t GetNPhot() {return fNphot;}
  Int_t GetPID() {return fPID;}

 private: 

  Double_t fX;
  Double_t fY;
  Double_t fZ;
  Double_t fT;
  Double_t fE;

  Int_t fCnum;
  Int_t fCnucl;
  Int_t fNgamma;
  Int_t fNphot;
  Int_t fPID;

  ClassDef(WCSimTrueCapture,0)

};

#endif
