#ifndef SANDBOXPMTCOVERAGE_HH
#define SANDBOXPMTCOVERAGE_HH

#include "TObject.h"

class SandBoxPMTcoverage : public TObject {

 public: 

  SandBoxPMTcoverage();

  ~SandBoxPMTcoverage();

  void SetBoxDimensions(double xl, double xh, double yl, double yh,
			double zl, double zh);

  void SetWallConfiguration(int mcase, int mmode, int mshape, double msize, double NbyN);

  bool isActiveHitOld(double mx, double my, double mz, int &whichPhotoSensor);

  bool isActiveHit(double mx, double my, double mz, int &whichPhotoSensor);

 private: 

  bool isActiveHitC(int mcase, double xx1, double xx2, int &whichPhotoSensor);

  Double_t NLAPPDs_per_rowcolumn;
  
  Double_t _xl,_yl,_zl,_xh,_yh,_zh;

  Double_t _NPDs_per_rowcolumn[6];
  Double_t _msize[6];
  Int_t _mmode[6];
  Int_t _mshape[6];

  ClassDef(SandBoxPMTcoverage,0)

};

#endif
