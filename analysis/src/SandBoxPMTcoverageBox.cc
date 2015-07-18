#include "SandBoxPMTcoverageBox.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"
#include <iostream>
#include <cmath>


ClassImp(SandBoxPMTcoverageBox)

SandBoxPMTcoverageBox::SandBoxPMTcoverageBox()
{
  _NLAPPDs_per_rowcolumn[0]=0;
  _NLAPPDs_per_rowcolumn[1]=0;
  _NLAPPDs_per_rowcolumn[2]=0;
  _NLAPPDs_per_rowcolumn[3]=0;
  _NLAPPDs_per_rowcolumn[4]=0;
  _NLAPPDs_per_rowcolumn[5]=0;

  _NPMTs_per_rowcolumn[0]=0;
  _NPMTs_per_rowcolumn[1]=0;
  _NPMTs_per_rowcolumn[2]=0;
  _NPMTs_per_rowcolumn[3]=0;
  _NPMTs_per_rowcolumn[4]=0;
  _NPMTs_per_rowcolumn[5]=0;

}

SandBoxPMTcoverageBox::~SandBoxPMTcoverageBox()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}

void SandBoxPMTcoverageBox::SetWallConfiguration(int mcase, int mmode, int mshape, double msize, double NbyN)
{
  // order of cases front,back,top,bottom,left,right;


  if(mcase==0){
    
    _mmode[0] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[0] = NbyN;
      _mLAPPDsize[0] = msize;
      
      std::cout<<"LAPPDs mode 0 "<<NbyN<<std::endl;
    } else{
      _NPMTs_per_rowcolumn[0] = NbyN;
      _mPMTsize[0] = msize;
    }
  }

  if(mcase==1){
    _mmode[1] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[1] = NbyN;
      _mLAPPDsize[1] = msize;

      std::cout<<"LAPPDs mode 1 "<<NbyN<<std::endl;
    } else{
      _NPMTs_per_rowcolumn[1] = NbyN;
      _mPMTsize[1] = msize;
    }
  }

  if(mcase==2){
    _mmode[2] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[2] = NbyN;
      _mLAPPDsize[2] = msize;
    } else{
      _NPMTs_per_rowcolumn[2] = NbyN;
      _mPMTsize[2] = msize;

      std::cout<<"LAPPDs mode 2 "<<NbyN<<std::endl;
    }
  }

  if(mcase==3){
    _mmode[3] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[3] = NbyN;
      _mLAPPDsize[3] = msize;
    } else{
      _NPMTs_per_rowcolumn[3] = NbyN;
      _mPMTsize[3] = msize;
    }
  }

  if(mcase==4){
     _mmode[4] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[4] = NbyN;
      _mLAPPDsize[4] = msize;
    } else{
      _NPMTs_per_rowcolumn[4] = NbyN;
      _mPMTsize[4] = msize;
    }
  }

  if(mcase==5){
     _mmode[5] = mmode;
    if(mshape==1){
      _NLAPPDs_per_rowcolumn[5] = NbyN;
      _mLAPPDsize[5] = msize;
    } else{
      _NPMTs_per_rowcolumn[5] = NbyN;
      _mPMTsize[5] = msize;
    }
  }
}


void SandBoxPMTcoverageBox::SetBoxDimensions(double xl, double xh, double yl, double yh,
				     double zl, double zh)
{
  _xl=xl;
  _xh=xh;
  _yl=yl;
  _yh=yh;
  _zl=zl;
  _zh=zh;
}


bool SandBoxPMTcoverageBox::isActiveHitLAPPD(double mx, double my, double mz, int &whichPhotoSensor){

  bool wasithit=false;
  int wP=-5555;
  int wcase=-4;
  double cxx1=-5555;
  double cxx2=-5555;
  if(mz==_zh){ wcase=0; wasithit=this->isActiveHitLAPPD_C(wcase, mx, my, wP, cxx1, cxx2); }
  if(mz==_zl){ wcase=1; wasithit=this->isActiveHitLAPPD_C(wcase, mx, my, wP, cxx1, cxx2); }
  if(my==_yh){ wcase=2; wasithit=this->isActiveHitLAPPD_C(wcase, mx, mz, wP, cxx1, cxx2); }
  if(my==_yl){ wcase=3; wasithit=this->isActiveHitLAPPD_C(wcase, mx, mz, wP, cxx1, cxx2); }
  if(mx==_xh){ wcase=4; wasithit=this->isActiveHitLAPPD_C(wcase, mz, my, wP, cxx1, cxx2); }
  if(mx==_xl){ wcase=5; wasithit=this->isActiveHitLAPPD_C(wcase, mz, my, wP, cxx1, cxx2); }

  whichPhotoSensor=wP;
  return wasithit;
}


bool SandBoxPMTcoverageBox::isActiveHitPMT(double mx, double my, double mz, int &whichPhotoSensor){

  bool wasithit=false;
  int wP=-5555;
  double cxx1=-5555;
  double cxx2=-5555;
  int wcase=-4;
  if(mz==_zh){ wcase=0; wasithit=this->isActiveHitPMT_C(wcase, mx, my, wP, cxx1, cxx2); }
  if(mz==_zl){ wcase=1; wasithit=this->isActiveHitPMT_C(wcase, mx, my, wP, cxx1, cxx2); }
  if(my==_yh){ wcase=2; wasithit=this->isActiveHitPMT_C(wcase, mx, mz, wP, cxx1, cxx2); }
  if(my==_yl){ wcase=3; wasithit=this->isActiveHitPMT_C(wcase, mx, mz, wP, cxx1, cxx2); }
  if(mx==_xh){ wcase=4; wasithit=this->isActiveHitPMT_C(wcase, mz, my, wP, cxx1, cxx2); }
  if(mx==_xl){ wcase=5; wasithit=this->isActiveHitPMT_C(wcase, mz, my, wP, cxx1, cxx2); }

  whichPhotoSensor=wP;
  return wasithit;
}


bool SandBoxPMTcoverageBox::isActiveHitLAPPD_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &centerxx1, double &centerxx2){

  bool washit=false;
  int PSnum=1;

  if(_mmode[mcase]==2) return true;
  if(_mmode[mcase]==0) return false;

  if(_mmode[mcase]==1 && _NLAPPDs_per_rowcolumn[mcase]>0){

    int npparity;
    if((int)_NLAPPDs_per_rowcolumn[mcase]%2==0) npparity=0;
    else npparity=1;

    double NPDs = _NLAPPDs_per_rowcolumn[mcase]*_NLAPPDs_per_rowcolumn[mcase];
    double lrowcolumn = 3000.0/_NLAPPDs_per_rowcolumn[mcase];

    double wxx1,wxx2;

    if(npparity==0){
      wxx1 = (double) floor(xx1/lrowcolumn);
      wxx2 = (double) floor(xx2/lrowcolumn);
      centerxx1 = (lrowcolumn*wxx1) + (lrowcolumn/2.0);
      centerxx2 = (lrowcolumn*wxx2) + (lrowcolumn/2.0);
    }

    if(npparity==1){
      wxx1 = (double) round(xx1/lrowcolumn);
      wxx2 = (double) round(xx2/lrowcolumn);
      centerxx1 = (lrowcolumn*wxx1);
      centerxx2 = (lrowcolumn*wxx2);
    }
    
    bool isX1range=false;
    bool isX2range=false;

    double d1 = fabs(centerxx1-xx1);
    double d2 = fabs(centerxx2-xx2);
    if( (d1<(_mLAPPDsize[mcase]/2)) && (d2<(_mLAPPDsize[mcase]/2) )){
      washit=true;
      PSnum=(int)(3*wxx1 + wxx2) + mcase*100;
    } 
    
    whichPhotoSensor=PSnum;  
  }  
  
  return washit;
}


bool SandBoxPMTcoverageBox::isActiveHitPMT_C(int mcase, double xx1, double xx2, int &whichPhotoSensor, double &centerxx1, double &centerxx2){

  bool washit=false;
  int PSnum=1;

  if(_mmode[mcase]==2) return true;
  if(_mmode[mcase]==0) return false;

  if(_mmode[mcase]==1 && _NPMTs_per_rowcolumn[mcase]>0){

    int npparity;
    if((int)_NPMTs_per_rowcolumn[mcase]%2==0) npparity=0;
    else npparity=1;

    double NPDs = _NPMTs_per_rowcolumn[mcase]*_NPMTs_per_rowcolumn[mcase];
    double lrowcolumn = 3000.0/_NPMTs_per_rowcolumn[mcase];

    double wxx1,wxx2;

    if(npparity==0){
      wxx1 = (double) floor(xx1/lrowcolumn);
      wxx2 = (double) floor(xx2/lrowcolumn);
      centerxx1 = (lrowcolumn*wxx1) + (lrowcolumn/2.0);
      centerxx2 = (lrowcolumn*wxx2) + (lrowcolumn/2.0);
    }

    if(npparity==1){
      wxx1 = (double) round(xx1/lrowcolumn);
      wxx2 = (double) round(xx2/lrowcolumn);
      centerxx1 = (lrowcolumn*wxx1);
      centerxx2 = (lrowcolumn*wxx2);
    }
    
    bool isX1range=false;
    bool isX2range=false;
    
    double r1sqr = (centerxx1-xx1)*(centerxx1-xx1);
    double r2sqr = (centerxx2-xx2)*(centerxx2-xx2);
    double distance = sqrt( r1sqr + r2sqr );
    if(distance<(_mPMTsize[mcase]/2)){
      washit=true;
      PSnum=(int)(_NPMTs_per_rowcolumn[mcase]*wxx1 + wxx2) + mcase*100;
    }
    

    whichPhotoSensor=PSnum;  
  }
   
  return washit;
}



bool SandBoxPMTcoverageBox::isActiveHitOld(double mx, double my, double mz, int &whichPhotoSensor){


  double NLAPPDs = NLAPPDs_per_rowcolumn*NLAPPDs_per_rowcolumn;

  double lrowcolumn = 3000.0/NLAPPDs_per_rowcolumn;

  double wx = (double) floor(mx/lrowcolumn);
  double wy = (double) floor(my/lrowcolumn);
  double wz = 0.;
  
  double centerx = (lrowcolumn*wx) + (lrowcolumn/2.0);
  double centery = (lrowcolumn*wy) + (lrowcolumn/2.0);
  double centerz = 1500.;

  bool applyQE=false;
  double QEval=0;  
  bool isXrange=false;
  bool isYrange=false;
  bool isZrange=false;

  if(mz==1500) isZrange=true;
  if( (mx<(centerx+101.6)) && (mx>(centerx-101.6)) ) isXrange=true;
  if( (my<(centery+101.6)) && (my>(centery-101.6)) ) isYrange=true;

  bool isHit = false;
  int PSnum=-1;
  if( isXrange && isYrange && isZrange ) { isHit=true;  PSnum=(int)(3*wy + wx); }

  whichPhotoSensor = PSnum;
  return isHit;
}



