#include "WCSimTruePart.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"

ClassImp(WCSimTruePart)

  WCSimTruePart::  WCSimTruePart(Double_t xs, Double_t ys, Double_t zs, Double_t ts, Double_t xe, Double_t ye, Double_t ze, Double_t te, Double_t pxs, Double_t pys, Double_t pzs, Double_t pxe, Double_t pye, Double_t pze, Double_t kes, Double_t kee, Int_t processs, Int_t processe, Int_t parentID, Int_t trackID, Int_t PID) 
{
  fXs=xs;
  fYs=ys;
  fZs=zs;
  fTs=ts;
  fXe=xe;
  fYe=ye;
  fZe=ze;
  fTe=te;

  fPXs=pxs;
  fPYs=pys;
  fPZs=pzs;
  fPXe=pxe;
  fPYe=pye;
  fPZe=pze;

  fKEs=kes;
  fKEe=kee;
    
  fprocessStart=processs;
  fprocessEnd=processe;
  fparentID=parentID;
  ftrackID=trackID;
  fPID=PID;

  WCSimRecoObjectTable::Instance()->NewDigit();
}

WCSimTruePart::~WCSimTruePart()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}

WCSimRecoDigit* WCSimTruePart::Convert2Reco()
{
  //take the # per PMT to be 1 (not correct for true hits)
  
  WCSimRecoDigit* ndig = new WCSimRecoDigit(0,(this->GetXend()),(this->GetYend()),(this->GetZend()),(this->GetTend()),1,1,(this->GetTend()),1);

  return ndig;
}
