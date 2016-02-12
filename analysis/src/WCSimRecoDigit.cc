#include "WCSimRecoDigit.hh"

#include "WCSimRecoObjectTable.hh"

ClassImp(WCSimRecoDigit)

WCSimRecoDigit::WCSimRecoDigit(Int_t region, Double_t x, Double_t y, Double_t z, Double_t rawT, Double_t rawQ, Int_t rawPEtube, Double_t calT, Double_t calQ)
{
  fRegion = region;

  fX = x;
  fY = y;
  fZ = z;
  
  fRawTime = rawT;
  fRawQPEs = rawQ;

  fCalTime = calT;  
  fCalQPEs = calQ;
  
  fRawPEtube = rawPEtube;

  fIsFiltered = 1; // okay by default

  WCSimRecoObjectTable::Instance()->NewDigit();
}

WCSimRecoDigit::~WCSimRecoDigit()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}
