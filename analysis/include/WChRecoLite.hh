/* Andrey Elagin, September 5, 2014
 * LightReco is a light-weight standalone vertex and (in the near future) directionality
 * reconsruction code for 0vbb-decay events in liquid scintillator. The code is based on
 * quadruplet-based vertex-finding method by Michael Smy. Many lines are directly copied
 * from WCSimAnalysis package.
 *
 * See README.txt for instructions on how to use and notes on significant updates.
 */


#ifndef WCHRECOLITE_HH
#define WCHRECOLITE_HH

#include "WCSimRecoDigit.hh"
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TError.h"
#include "TMinuit.h"
#include "TRandom.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

class WChRecoLite : public TObject {

 public: 

  static WChRecoLite* Instance(); 

  Int_t FindVertex(Double_t x0, Double_t y0, Double_t z0, Double_t t0, Double_t x1, Double_t y1, Double_t z1, Double_t t1, Double_t x2, Double_t y2, Double_t z2, Double_t t2, Double_t x3, Double_t y3, Double_t z3, Double_t t3, Double_t& vxm, Double_t& vym, Double_t& vzm, Double_t& vtm, Double_t& vxp, Double_t& vyp, Double_t& vzp, Double_t& vtp);

  Int_t ChooseNextDigit(Double_t& xpos, Double_t& ypos, Double_t& zpos, Double_t& time);

  Int_t ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0, Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1, Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2, Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3);

  Int_t CalcVertexSeeds();

  Int_t TimePropertiesLnL(double & vtx_time, double & fom);

  void FitPointTimePropertiesLnL(Double_t& fit_time, Double_t& fom);

  Int_t SelectBestSeed(int evt_num);

  std::vector<double> GetSeedVtx(int wvert);

  void SetDigits(std::vector<double> iDigitX,std::vector<double> iDigitY,std::vector<double> iDigitZ,std::vector<double> iDigitT,std::vector<double> iDigitQ,std::vector<double> iDigitPE,std::vector<double> iDigitW,std::vector<double> iDigitV,std::vector<double> iDelta,std::vector<int> iSeedDigitList);

 private: 

  WChRecoLite();
  ~WChRecoLite();

  void Reset();

  double R_SPHERE; //sphere diameter [cm]
  double N_REF; //average index of refraction
  double C_VAC; //speed of light in vacuum [cm/ns]
  int NSeedsTarget; //number of quadruplets
  double TSIGMA; //total time spread (including detector TTS chromatic dispersions)

  double fBaseFOM; //Figure of merit. Borrowed from WCSim: 
                   //the higher it is the better
  double meanTime;
  double seedTime;

  int fNDigits;
  int fThisDigit;
  int fLastEntry;
  int fCounter;
  int fMinTime;
  
  // store photon hits after filtering cuts (e.g. position dependent cut to increase cherenkov fraction)
  std::vector<double> fDigitX;
  std::vector<double> fDigitY;
  std::vector<double> fDigitZ;
  std::vector<double> fDigitT;
  std::vector<double> fDigitQ;
  std::vector<double> fDigitPE;
  std::vector<double> fDigitW;
  std::vector<double> fDigitV;
  std::vector<double> fDelta; // time residual
  
  std::map<int, double> INDEX;

  TRandom RND;

  // store seed vertex calculated from quaruplets
  std::vector<double> vSeedVtxX;
  std::vector<double> vSeedVtxY;
  std::vector<double> vSeedVtxZ;
  std::vector<double> vSeedVtxTime;
  std::vector<int> vSeedDigitList;
  
  TFile* fFOM;
  TH1F* hT;
  TH1F* hDT0;
  TH1F* hDT;
  
  ClassDef(WChRecoLite,0)

};

#endif
