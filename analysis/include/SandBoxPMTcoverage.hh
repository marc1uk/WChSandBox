#ifndef SANDBOXPMTCOVERAGE_HH
#define SANDBOXPMTCOVERAGE_HH

#include "TObject.h"

class SandBoxPMTcoverage : public TObject {

 public: 

  SandBoxPMTcoverage();

  ~SandBoxPMTcoverage();

//  void SetBoxDimensions(double xl, double xh, double yl, double yh,
//			double zl, double zh);

  void SetWallConfiguration(int mcase, int mmode, int mshapePMT, double msizePMT, double NrowPMT, double NcolPMT, double QE_PMT, int mshapeLAPPD, double msizeLAPPD, double NrowLAPPD, double NcolLAPPD, double QE_LAPPD);


  void SetBoxDimensions(int opt, double dim1, double dim2);

  int isActiveHit(double mx, double my, double mz, int &whichPhotoSensor, bool LAPPDs);

  int isActiveHitC(int mcase, double xx1, double xx2, int &whichPhotoSensor, bool LAPPDs);

  bool isActiveHitOld(double mx, double my, double mz, int &whichPhotoSensor);



 private: 

  bool isActiveHitLAPPD_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &cxx1, double &cxx2);
  bool isActiveHitPMT_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &cxx1, double &cxx2);

  Double_t NLAPPDs_per_rowcolumn;

  Int_t _mmode[7];
 
  Double_t _PMTs_per_row[7];
  Double_t _PMTs_per_column[7];
  Double_t _QE_PMT[7];
  Int_t _PMTs_mshape[7];
  Double_t _PMTs_msize[7];

  Double_t _LAPPDs_per_row[7];
  Double_t _LAPPDs_per_column[7];
  Double_t _QE_LAPPD[7];
  Int_t _LAPPDs_mshape[7];
  Double_t _LAPPDs_msize[7];

  Double_t _xl,_yl,_zl,_xh,_yh,_zh;  
  Int_t _opt;
  Double_t _Rh;

  ClassDef(SandBoxPMTcoverage,0)

};

#endif
