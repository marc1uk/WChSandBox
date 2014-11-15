#ifndef WCSIMISOCHRON_HH
#define WCSIMISOCHRON_HH

#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraph.h"
#include "TObject.h"
#include "TVector3.h"
#include "WCSimRecoCluster.hh"
#include "WCSimWaterModel.hh"
//#include "RasterizeCircle.hh"
#include <iostream>
#include <vector>

using namespace std;

class WCSimIsochron : public TObject {

 public: 

  WCSimIsochron();
  WCSimIsochron(Double_t nBinsV, Double_t binscale);
  ~WCSimIsochron();

  void Reset();
  void ApplyTransform(vector< vector<double> > thehits, vector<double> VxHypothesis);
  void ProcessHit(vector<double> aHit, vector<double> VxHypothesis);
  void SetConstantIndexofRefraction(double newn);
  void SetConstantSpeedofParticle(double newc);
  void AddWaterModel(WCSimWaterModel* wm);

  Double_t CalcMaxAlpha(Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcMaxEmissionAngle(Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcS1(Double_t alpha, Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcS2(Double_t alpha, Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcAlpha(Double_t s1, Double_t dD, Double_t dT, Double_t theindex);

  TGraph* IndexVsS1(Double_t dD, Double_t dT);

 private: 

  WCSimWaterModel *_theWM;

  vector< vector<double> > isochronarray;

  double _n,_c;

  int nb2x,nb2y;

  int nbx;
  int nby;
  int nbz;
  double lx;
  double hx;
  double ly;
  double hy;
  double lz;
  double hz;

  bool iswatermodel;
  bool usewatermodel;

  ClassDef(WCSimIsochron,0)
};

#endif
