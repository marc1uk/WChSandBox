#ifndef WCSIMTRUELIGHT_HH
#define WCSIMTRUELIGHT_HH

#include "WCSimRecoDigit.hh"
#include "TObject.h"

class WCSimTrueLight : public TObject {

 public: 

  WCSimTrueLight(Double_t xs, Double_t ys, Double_t zs, Double_t ts, Double_t xe, Double_t ye, Double_t ze, Double_t te, Double_t wavelength, Int_t processs, Int_t isScat, Int_t parentID, Int_t trackID, Int_t ishit, Int_t capnum, Int_t PMTid); 

  ~WCSimTrueLight();

  Double_t GetXstart() { return fXs; }
  Double_t GetYstart() { return fYs; }
  Double_t GetZstart() { return fZs; }
  Double_t GetTstart() { return fTs; }
  Double_t GetXend() { return fXe; }
  Double_t GetYend() { return fYe; }
  Double_t GetZend() { return fZe; }
  Double_t GetTend() { return fTe; }
  Double_t GetWavelength() { return fWavelength; }

  Int_t GetIsScat() { return fisScattered;}
  Int_t GetIsHit() {return fisHit;}
  Int_t GetProcessStart() {return fprocessStart;}
  Int_t GetParentID() {return fparentID;}
  Int_t GetTrackID() {return ftrackID;}
  Int_t GetCaptureNum() {return fcaptNum;}
  Int_t GetPMTID() {return fPMTID;}

  Double_t GetPolarization1() { return fP1; }
  Double_t GetPolarization2() { return fP2; }

  WCSimRecoDigit* Convert2Reco();

 private: 

  Double_t fXs;
  Double_t fYs;
  Double_t fZs;
  Double_t fTs;
  Double_t fXe;
  Double_t fYe;
  Double_t fZe;
  Double_t fTe;
  Double_t fWavelength;

  Int_t fisScattered;
  Int_t fisHit;
  Int_t fprocessStart;
  Int_t fparentID;
  Int_t ftrackID;
  Int_t fPMTID;
  Int_t fcaptNum;


  Double_t fThet;
  Double_t fPhi;
  Double_t fP1;
  Double_t fP2;

  ClassDef(WCSimTrueLight,0)

};

#endif
