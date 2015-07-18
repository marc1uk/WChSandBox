#ifndef RASTERIZECIRCLE_HH
#define RASTERIZECIRCLE_HH

#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class RasterizeCircle : public TObject {

 public:
  RasterizeCircle(Int_t nbx, Double_t lx, Double_t hx, Int_t nby, Double_t ly, Double_t hy,Int_t nbz, Double_t lz, Double_t hz);
  ~RasterizeCircle();

  TH3D* returnhisto();
  Int_t Reset();
  Int_t Rasterize3D(std::vector<double> hvect, Double_t thealpha, Double_t theS1);
  Int_t Rasterize3D(std::vector<double> hvect, Double_t thealpha, Double_t theS1, Double_t wt);
  Int_t Rasterize3Dspeedorder(TVector3 hvect, Double_t thealpha, Double_t theS1, Int_t tfastb, Int_t tmedb, Int_t tslowb, Double_t wt);
  Int_t Rasterize3Dflat(Double_t R, Double_t tn0, Double_t tn1, Double_t tn2, Double_t wt);
  Int_t InnerLoops(Int_t init0bin, Int_t init1bin, Int_t init2bin, Int_t d1minBin, Int_t d2minBin, Int_t d1maxBin, Int_t d2maxBin, Double_t cV0,Double_t cV1, Double_t cV2, Double_t theR, Int_t incr1, Double_t wt);
  Int_t InnerMostLoop(Int_t init0bin, Int_t init1bin, Int_t init2bin, Int_t d1minBin, Int_t d2minBin, Int_t d1maxBin, Int_t d2maxBin, Double_t cV0,Double_t cV1, Double_t cV2, Double_t theR, Int_t incr2, Double_t wt);
  Int_t FillingAction(Int_t bin0, Int_t bin1, Int_t bin2, Double_t wt);
  Int_t SolveD1D2(Double_t theR, Double_t cV0, Double_t cV1, Double_t cV2, Double_t d0coor, Double_t &hid1, Double_t &hid2 , Double_t &lowd1, Double_t &lowd2);
  Int_t SolveForFixedCoordinate(Double_t theR, Double_t cVFixed, Double_t cV1, Double_t cV2, Double_t d0coor, Double_t d1coor, Double_t d2coor, Double_t &thed1, Double_t &thed2);
  TH2D* Rasterize2D(Int_t R, Int_t x0, Int_t y0);


 private:

  Int_t fRegion;

  Int_t fastb,medb,slowb;

  Int_t nbins;
  Double_t lr,hr;

  Double_t binsize;

  TH1D* binfinder;
  TH3D* circlehist;

  ClassDef(RasterizeCircle,0)

};

#endif
