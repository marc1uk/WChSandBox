#include "SandBoxPMTcoverage.hh"
#include "TObject.h"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoObjectTable.hh"
#include <iostream>
#include <cmath>


ClassImp(SandBoxPMTcoverage)

SandBoxPMTcoverage::SandBoxPMTcoverage()
{

}

SandBoxPMTcoverage::~SandBoxPMTcoverage()
{
  WCSimRecoObjectTable::Instance()->DeleteDigit();
}

void SandBoxPMTcoverage::SetWallConfiguration(int mcase, int mmode, int mshape, double msize, double NbyN)
{
  // order of cases front,back,top,bottom,left,right;

  if(mcase==0){
    _NPDs_per_rowcolumn[0] = NbyN;
    _mmode[0] = mmode;
    _msize[0] = msize;
    _mshape[0] = mshape;
  }

  if(mcase==1){
    _NPDs_per_rowcolumn[1] = NbyN;
    _mmode[1] = mmode;
    _msize[1] = msize;
    _mshape[1] = mshape;
  }

  if(mcase==2){
    _NPDs_per_rowcolumn[2] = NbyN;
    _mmode[2] = mmode;
    _msize[2] = msize;
    _mshape[2] = mshape;
  }

  if(mcase==3){
    _NPDs_per_rowcolumn[3] = NbyN;
    _mmode[3] = mmode;
    _msize[3] = msize;
    _mshape[3] = mshape;
  }

  if(mcase==4){
    _NPDs_per_rowcolumn[4] = NbyN;
    _mmode[4] = mmode;
    _msize[4] = msize;
    _mshape[4] = mshape;
  }

  if(mcase==5){
    _NPDs_per_rowcolumn[5] = NbyN;
    _mmode[5] = mmode;
    _msize[5] = msize;
    _mshape[5] = mshape;
  }
}


void SandBoxPMTcoverage::SetBoxDimensions(double xl, double xh, double yl, double yh,
				     double zl, double zh)
{
  _xl=xl;
  _xh=xh;
  _yl=yl;
  _yh=yh;
  _zl=zl;
  _zh=zh;
}


bool SandBoxPMTcoverage::isActiveHit(double mx, double my, double mz, int &whichPhotoSensor){

  bool wasithit=false;
  int wP=-5555;
  int wcase=-4;
  if(mz==_zh){ wcase=0; wasithit=this->isActiveHitC(wcase, mx, my, wP); }
  if(mz==_zl){ wcase=1; wasithit=this->isActiveHitC(wcase, mx, my, wP); }
  if(my==_yh){ wcase=2; wasithit=this->isActiveHitC(wcase, mx, mz, wP); }
  if(my==_yl){ wcase=3; wasithit=this->isActiveHitC(wcase, mx, mz, wP); }
  if(mx==_xh){ wcase=4; wasithit=this->isActiveHitC(wcase, mz, my, wP); }
  if(mx==_xl){ wcase=5; wasithit=this->isActiveHitC(wcase, mz, my, wP); }

  whichPhotoSensor=wP;
  return wasithit;
}


bool SandBoxPMTcoverage::isActiveHitC(int mcase, double xx1, double xx2, int &whichPhotoSensor){

  bool washit=false;
  int PSnum=1;

  if(_mmode[mcase]==2) return true;
  if(_mmode[mcase]==0) return false;

  if(_mmode[mcase]==1){

    int npparity;
    if((int)_NPDs_per_rowcolumn[mcase]%2==0) npparity=0;
    else npparity=1;

    double NPDs = _NPDs_per_rowcolumn[mcase]*_NPDs_per_rowcolumn[mcase];
    double lrowcolumn = 3000.0/_NPDs_per_rowcolumn[mcase];

    double wxx1,wxx2,centerxx1,centerxx2;

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
    
    if(_mshape[mcase]==1)
      {
	double d1 = fabs(centerxx1-xx1);
	double d2 = fabs(centerxx2-xx2);
	if( (d1<(_msize[mcase]/2)) && (d2<(_msize[mcase]/2) )){
	  washit=true;
	  PSnum=(int)(3*wxx1 + wxx2) + mcase*100;
	} 
      }

    if(_mshape[mcase]==0) 
      {
      
	double r1sqr = (centerxx1-xx1)*(centerxx1-xx1);
	double r2sqr = (centerxx2-xx2)*(centerxx2-xx2);
	double distance = sqrt( r1sqr + r2sqr );
	if(distance<(_msize[mcase]/2)){
	  washit=true;
	  PSnum=(int)(3*wxx1 + wxx2) + mcase*100;
	}
      }
    whichPhotoSensor=PSnum;  
  }
  
  
    return washit;
}


bool SandBoxPMTcoverage::isActiveHitOld(double mx, double my, double mz, int &whichPhotoSensor){


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



