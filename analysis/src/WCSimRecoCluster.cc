#include "WCSimRecoCluster.hh"
#include "WCSimRecoDigit.hh"

#include "WCSimRecoObjectTable.hh"

#include <algorithm>

ClassImp(WCSimRecoCluster)

WCSimRecoCluster::WCSimRecoCluster()
{
  WCSimRecoObjectTable::Instance()->NewCluster();
}

WCSimRecoCluster::~WCSimRecoCluster()
{
  WCSimRecoObjectTable::Instance()->DeleteCluster();
}

void WCSimRecoCluster::Reset()
{  
  fDigitList.clear();
}

static bool CompareTimes(WCSimRecoDigit *rd1, WCSimRecoDigit *rd2)
{
  return ( rd1->GetTime() > rd2->GetTime() );
}

void WCSimRecoCluster::SortCluster()
{
  sort(fDigitList.begin(), fDigitList.end(), CompareTimes);
}

void WCSimRecoCluster::AddDigit(WCSimRecoDigit* digit)
{
  fDigitList.push_back(digit);
}

WCSimRecoDigit* WCSimRecoCluster::GetDigit(Int_t n)
{
  return (WCSimRecoDigit*)(fDigitList.at(n));
}
  
Int_t WCSimRecoCluster::GetNDigits()
{
  return fDigitList.size();
}
