#ifndef WCSIMISOCHRONTRANSFORM_HH
#define WCSIMISOCHRONTRANSFORM_HH

#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraph.h"
#include "TObject.h"
#include "TVector3.h"
#include "WCSimRecoCluster.hh"
#include "WCSimWaterModel.hh"
#include "RasterizeCircle.hh"
#include <iostream>
#include <vector>

using namespace std;

class WCSimIsoChronTransform : public TObject {

 public: 

  WCSimIsoChronTransform();
  WCSimIsoChronTransform(Double_t nBinsV, Double_t binscale);
  ~WCSimIsoChronTransform();

  void Reset();
  void ApplyTransform(vector< vector<double> > thehits, vector<double> VxHypothesis);
  void ApplyTransformNoRing(vector< vector<double> > thehits, vector<double> VxHypothesis);
  void ApplyTransformPMTres(vector< vector<double> > thehits, vector<double> VxHypothesis, double theresolution);
  void ProcessHit(vector<double> aHit, vector<double> VxHypothesis);
  void ProcessHitPMTres(vector<double> hitcoordinates, vector<double> hypVtx, double theresolution);
  void ProcessHitNoRing(vector<double> hitcoordinates, vector<double> hypVtx);
  void SetConstantIndexofRefraction(double newn);
  void SetConstantSpeedofParticle(double newc);
  void AddWaterModel(WCSimWaterModel* wm);
  Int_t FillIsoChronRing_raster(TVector3 hvect, Double_t alpha, Double_t s1, Double_t theweight);
  Double_t CalcMaxAlpha(Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcMaxEmissionAngle(Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcS1(Double_t alpha, Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcS2(Double_t alpha, Double_t dD, Double_t dT, Double_t theindex);
  Double_t CalcAlpha(Double_t s1, Double_t dD, Double_t dT, Double_t theindex);

  TGraph* IndexVsS1(Double_t dD, Double_t dT);
  TH1D* getBCDist3D();
  TH1D* getBCDist2D();
  TH1D* getMaxAlphaDist();
  TH1D* getMaxThetaCDist();
  TH1D* getMinS1Dist();
  TH1D* getnDist();
  TH2D* getS1vsAlpha(double thresh);

  TH3D* Get3Dhisto(double bclimit);
  TH3D* GetRaw3dIsoChr();

  TH2D* XYProjection();
  TH2D* XZProjection(double bclimit);
  TH2D* YZProjection();

  TH2D* XYSlice();
  TH2D* XZSlice();
  TH2D* YZSlice();

 private: 

  WCSimWaterModel *_theWM;
  RasterizeCircle *_theIsochronRaster;

  vector< vector<double> > isochronarray;

  //  TH3D* _htrackreco;
  //  TH3D** _htrackreco_color;
  //  TH3D** _htrackreco_resolution;
  double _n,_c;
  TH1D* _thexb;
  TH1D* _theyb;
  TH1D* _thezb; 
  TH2D* theprojxz;
  TH2D* theprojxy;
  TH2D* theprojyz;

  TH1D* _maxAlphadist;
  TH1D* _maxThetaCdist;
  TH1D* _bcontentsdist3d;
  TH1D* _bcontentsdist2d;
  TH1D* _minS1dist;
  TH1D* _nDist;
  TH2D* _S1vsAlphadist;
  TH2D* _returnS1vsAlpha;

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

  ClassDef(WCSimIsoChronTransform,0)
};

#endif
