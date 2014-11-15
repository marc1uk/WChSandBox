#include "WCSimTrueCapture.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"

ClassImp(WCSimTrueCapture)

  WCSimTrueCapture::WCSimTrueCapture(Double_t x, Double_t y, Double_t z, Double_t t, Double_t E, Int_t capnum, Int_t capnuc, Int_t pid, Int_t ngamma, Int_t nphot)
{
  fX=x;
  fY=y;
  fZ=z;
  fT=t;
  fE=E;

  fPID=pid;
  fNphot=nphot;
  fNgamma = ngamma;
  fCnucl = capnuc;
  fCnum = capnum;
  

  WCSimRecoObjectTable::Instance()->NewDigit();
}

WCSimTrueCapture::~WCSimTrueCapture()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}
