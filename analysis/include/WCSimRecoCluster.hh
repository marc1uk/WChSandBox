#ifndef WCSIMRECOCLUSTER_HH
#define WCSIMRECOCLUSTER_HH

#include "TObject.h"

#include <vector>

class WCSimRecoDigit;

class WCSimRecoCluster : public TObject {

 public:
  WCSimRecoCluster();                  
  ~WCSimRecoCluster();

  void Reset();
  void SortCluster();

  void AddDigit(WCSimRecoDigit* digit);
 
  WCSimRecoDigit* GetDigit(Int_t n);
  Int_t GetNDigits();

 private:

  std::vector<WCSimRecoDigit*> fDigitList;

  ClassDef(WCSimRecoCluster,0)

};

#endif







