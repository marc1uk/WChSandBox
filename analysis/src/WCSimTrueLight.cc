#include "WCSimTrueLight.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"

ClassImp(WCSimTrueLight)

  WCSimTrueLight::WCSimTrueLight(Double_t xs, Double_t ys, Double_t zs, Double_t ts, Double_t xe, Double_t ye, Double_t ze, Double_t te, Double_t wavelength, Int_t processs, Int_t isScat, Int_t parentID, Int_t trackID, Int_t ishit, Int_t capnum, Int_t PMTid) 
{
  fXs=xs;
  fYs=ys;
  fZs=zs;
  fTs=ts;
  fXe=xe;
  fYe=ye;
  fZe=ze;
  fTe=te;
  fWavelength=wavelength;  

  fisScattered=isScat;
  fisHit=ishit;
  fprocessStart=processs;
  fparentID=parentID;
  ftrackID=trackID;
  fcaptNum=capnum;
  fPMTID=PMTid;

  fThet=-55555;
  fPhi=-55555;
  fP1=-55555;
  fP2=-55555;

  WCSimRecoObjectTable::Instance()->NewDigit();
}

WCSimTrueLight::~WCSimTrueLight()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}

WCSimRecoDigit* WCSimTrueLight::Convert2Reco()
{
  //take the # per PMT to be 1 (not correct for true hits)
  
  WCSimRecoDigit* ndig = new WCSimRecoDigit(0,(this->GetXend()),(this->GetYend()),(this->GetZend()),(this->GetTend()),1,1,(this->GetTend()),1);

  return ndig;
}
