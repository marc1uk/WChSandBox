#ifndef WCSIMTRUEPART_HH
#define WCSIMTRUEPART_HH

#include "WCSimRecoDigit.hh"
#include "TObject.h"

class WCSimTruePart : public TObject {

 public: 

  WCSimTruePart(Double_t xs, Double_t ys, Double_t zs, Double_t ts, Double_t xe, Double_t ye, Double_t ze, Double_t te, Double_t pxs, Double_t pys, Double_t pzs, Double_t pxe, Double_t pye, Double_t pze, Double_t kes, Double_t kee, Int_t processs, Int_t processe, Int_t parentID, Int_t trackID, Int_t PID); 

  ~WCSimTruePart();

  Double_t GetXstart() { return fXs; }
  Double_t GetYstart() { return fYs; }
  Double_t GetZstart() { return fZs; }
  Double_t GetTstart() { return fTs; }
  Double_t GetXend() { return fXe; }
  Double_t GetYend() { return fYe; }
  Double_t GetZend() { return fZe; }
  Double_t GetTend() { return fTe; }
  Double_t GetPXstart() { return fPXs; }
  Double_t GetPYstart() { return fPYs; }
  Double_t GetPZstart() { return fPZs; }
  Double_t GetPXend() { return fPXe; }
  Double_t GetPYend() { return fPYe; }
  Double_t GetPZend() { return fPZe; }
  Double_t GetKEstart() { return fKEs; }
  Double_t GetKEend() { return fKEe; }

  Int_t GetProcessStart() {return fprocessStart;}
  Int_t GetProcessEnd() {return fprocessEnd;}
  Int_t GetParentID() {return fparentID;}
  Int_t GetTrackID() {return ftrackID;}
  Int_t GetPID() {return fPID;}

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
  Double_t fPXs;
  Double_t fPYs;
  Double_t fPZs;
  Double_t fPXe;
  Double_t fPYe;
  Double_t fPZe;
  Double_t fKEs;
  Double_t fKEe;

  Int_t fprocessStart;
  Int_t fprocessEnd;
  Int_t fparentID;
  Int_t ftrackID;
  Int_t fPID;

  ClassDef(WCSimTruePart,0)

};

#endif
