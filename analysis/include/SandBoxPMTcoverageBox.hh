#ifndef SANDBOXPMTCOVERAGEBOX_HH
#define SANDBOXPMTCOVERAGEBOX_HH

#include "TObject.h"

class SandBoxPMTcoverageBox : public TObject {

 public: 

  SandBoxPMTcoverageBox();

  ~SandBoxPMTcoverageBox();

  void SetBoxDimensions(double xl, double xh, double yl, double yh,
			double zl, double zh);

  void SetWallConfiguration(int mcase, int mmode, int mshape, double msize, double NbyN);

  bool isActiveHitOld(double mx, double my, double mz, int &whichPhotoSensor);

  bool isActiveHitLAPPD(double mx, double my, double mz, int &whichPhotoSensor);
  bool isActiveHitPMT(double mx, double my, double mz, int &whichPhotoSensor);

 private: 

  bool isActiveHitLAPPD_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &cxx1, double &cxx2);
  bool isActiveHitPMT_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &cxx1, double &cxx2);

  Double_t NLAPPDs_per_rowcolumn;
  
  Double_t _xl,_yl,_zl,_xh,_yh,_zh;

  Double_t _NLAPPDs_per_rowcolumn[6];
  Double_t _NPMTs_per_rowcolumn[6];
  Double_t _mPMTsize[6];
  Double_t _mLAPPDsize[6];
  Int_t _mmode[6];

  ClassDef(SandBoxPMTcoverageBox,0)

};

#endif
