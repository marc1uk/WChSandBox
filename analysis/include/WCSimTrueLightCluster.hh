#ifndef WCSIMTRUELIGHTCLUSTER_HH
#define WCSIMTRUELIGHTCLUSTER_HH

#include "WCSimTrueLight.hh"
#include "WCSimRecoCluster.hh"
#include "TObject.h"

class WCSimTrueLightCluster : public TObject {

 public: 

  WCSimTrueLightCluster();
  ~WCSimTrueLightCluster();

  void Reset();

  void SortCluster();
  void AddDigit(WCSimTrueLight* digit);
  WCSimTrueLight* GetDigit(Int_t n);
  WCSimRecoCluster* Convert2RecoCluster();
  Int_t GetNDigits();

 private: 

  std::vector<WCSimTrueLight*> fTrueLightList;

  ClassDef(WCSimTrueLightCluster,0)

};

#endif
