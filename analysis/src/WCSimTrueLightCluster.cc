#include "WCSimTrueLight.hh"
#include "WCSimTrueLightCluster.hh"
#include "WCSimRecoCluster.hh"
#include "WCSimRecoObjectTable.hh"

#include <algorithm>

ClassImp(WCSimTrueLightCluster)

  WCSimTrueLightCluster::WCSimTrueLightCluster()
{
  WCSimRecoObjectTable::Instance()->NewCluster();
}

WCSimTrueLightCluster::~WCSimTrueLightCluster()
{
  WCSimRecoObjectTable::Instance()->DeleteCluster();
}

void WCSimTrueLightCluster::Reset()
{
  fTrueLightList.clear();
}

static bool CompareTimes(WCSimTrueLight *tl1, WCSimTrueLight *tl2)
{
  return ( tl1->GetTend() > tl2->GetTend() );
}

void WCSimTrueLightCluster::SortCluster()
{
  std::sort(fTrueLightList.begin(), fTrueLightList.end(), CompareTimes);
}

void WCSimTrueLightCluster::AddDigit(WCSimTrueLight* digit)
{
  fTrueLightList.push_back(digit);
}

WCSimTrueLight* WCSimTrueLightCluster::GetDigit(Int_t n)
{
  return (WCSimTrueLight*)(fTrueLightList.at(n));
}

Int_t WCSimTrueLightCluster::GetNDigits()
{
  return fTrueLightList.size();
}

WCSimRecoCluster* WCSimTrueLightCluster::Convert2RecoCluster()
{
  WCSimRecoCluster* _rcluster = new WCSimRecoCluster();

  for(int i=0; i<((int)fTrueLightList.size()); i++){

    WCSimRecoDigit* digiti =(fTrueLightList.at(i))->Convert2Reco();
    _rcluster->AddDigit(digiti);
  }

  return _rcluster;
}

// Local Variables: **
// c-basic-offset: **
// End: **
