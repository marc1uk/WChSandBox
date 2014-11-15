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


void SandBoxPMTcoverage::SetWallConfiguration(int mcase, int mmode, int mshapePMT, double msizePMT, double NrowPMT, double NcolPMT, double QE_PMT, int mshapeLAPPD, double msizeLAPPD, double NrowLAPPD, double NcolLAPPD, double QE_LAPPD)
{
  // order of cases : cylindar, xh, xl, yh, yl, zh, zl; 
  for(int ii=0;ii<7;ii++){
    if(mcase==ii){
      _mmode[ii] = mmode;
      _QE_PMT[ii] = QE_PMT;
      _QE_LAPPD[ii] = QE_LAPPD;
      _PMTs_per_row[ii] = NrowPMT;
      _PMTs_per_column[ii] = NcolPMT;
      _PMTs_msize[ii] = msizePMT;
      _PMTs_mshape[ii] = mshapePMT;
      _LAPPDs_per_row[ii] = NrowLAPPD;
      _LAPPDs_per_column[ii] = NcolLAPPD;
      _LAPPDs_msize[ii] = msizeLAPPD;
      _LAPPDs_mshape[ii] = mshapeLAPPD;
    }
  }
}

void SandBoxPMTcoverage::SetBoxDimensions(int opt, double dim1, double dim2)
{
  _opt=opt;

  if(_opt==0){
  _Rh=dim2;
  _zh=dim1;
  _zl=-dim1;
  }

  if(_opt==1){
  _xh=dim1;
  _xl=-dim1;
  _yh=dim1;
  _yl=-dim1;
  _zh=dim1;
  _zl=-dim1;
  }
}


int SandBoxPMTcoverage::isActiveHit(double mx, double my, double mz, int &whichPhotoSensor, bool LAPPDs){

  int wasithit=0;
  int wP=-5555;
  int wcase=-4;
  double mtheta = atan2(my,mx)*_Rh;
  if((mx*mx + my*my)==(_Rh*_Rh) && _opt==0){ wcase=0; wasithit=this->isActiveHitC(wcase, mz, mtheta, wP, LAPPDs); }
  if(mx==_xh && _opt==1){ wcase=1; wasithit=this->isActiveHitC(wcase, my, mz, wP, LAPPDs); }
  if(mx==_xl && _opt==1){ wcase=2; wasithit=this->isActiveHitC(wcase, my, mz, wP, LAPPDs); }
  if(my==_yh && _opt==1){ wcase=3; wasithit=this->isActiveHitC(wcase, mx, mz, wP, LAPPDs); }
  if(my==_yl && _opt==1){ wcase=4; wasithit=this->isActiveHitC(wcase, mx, mz, wP, LAPPDs); }
  if(mz==_zh){ wcase=5; wasithit=this->isActiveHitC(wcase, mx, my, wP, LAPPDs); }
  if(mz==_zl){ wcase=6; wasithit=this->isActiveHitC(wcase, mx, my, wP, LAPPDs); }

  whichPhotoSensor=wP;
  return wasithit;
}


int SandBoxPMTcoverage::isActiveHitC(int mcase, double xx1, double xx2, int &whichPhotoSensor, bool LAPPDs){

  int washit=0;
  int PSnum=1;

  if(_mmode[mcase]==2) return 2;
  if(_mmode[mcase]==0) return 0;

  if(_mmode[mcase]==1){

    int npparityrow;
    if((int)_PMTs_per_row[mcase]%2==0) npparityrow=0;
    else npparityrow=1;
    
    int npparitycolumn;
    if((int)_PMTs_per_column[mcase]%2==0) npparitycolumn=0;
    else npparitycolumn=1;

    double PMTs = _PMTs_per_row[mcase]*_PMTs_per_column[mcase];
    double lrow;
    double lcolumn;

    if( _opt==0 && mcase==0 ) {
      lrow = (2*_zh)/_PMTs_per_row[mcase];
      lcolumn = (2*3.1416*_Rh)/_PMTs_per_column[mcase];
    }
    if( _opt==0 && (mcase==5 || mcase==6) ){ 
      lrow = (2*_Rh)/_PMTs_per_row[mcase];
      lcolumn = (2*_Rh)/_PMTs_per_column[mcase]; 
    }
    if( _opt==1 ){
      lrow = (2*_xh)/_PMTs_per_row[mcase];
      lcolumn = (2*_xh)/_PMTs_per_column[mcase];
    }

    double wxx1,wxx2,centerxx1,centerxx2;

    if(npparityrow==0){
      wxx1 = (double) floor(xx1/lrow);
      centerxx1 = (lrow*wxx1) + (lrow/2.0);
    }

    if(npparityrow==1){
      wxx1 = (double) round(xx1/lrow);
      centerxx1 = (lrow*wxx1);
    }
    
    if(npparitycolumn==0){
      wxx2 = (double) floor(xx2/lcolumn);
      centerxx2 = (lcolumn*wxx2) + (lcolumn/2.0);
    }

    if(npparitycolumn==1){
      wxx2 = (double) round(xx2/lcolumn);
      centerxx2 = (lcolumn*wxx2);
    }

    bool isX1range=false;
    bool isX2range=false;
    double rd = (rand()/(double)RAND_MAX);

    if(_PMTs_mshape[mcase]==1)
    {
      double d1 = fabs(centerxx1-xx1);
      double d2 = fabs(centerxx2-xx2);
      double d3 = (fabs(centerxx1)+_PMTs_msize[mcase]/2)*(fabs(centerxx1)
			+_PMTs_msize[mcase]/2) + (fabs(centerxx2)+_PMTs_msize[mcase]/2)
			*(fabs(centerxx2)+_PMTs_msize[mcase]/2);
      if( d1<(_PMTs_msize[mcase]/2) && d2<(_PMTs_msize[mcase]/2) 
			&& (d3<(_Rh*_Rh) || mcase==0 || _opt==1) && rd<(_QE_PMT[mcase]/100.) ){
	washit=1;
	PSnum=(int)(wxx1 + wxx2*100) + mcase*10000;
      }	
    }

    if(_PMTs_mshape[mcase]==0) 
    {
      double r1sqr = (centerxx1-xx1)*(centerxx1-xx1);
      double r2sqr = (centerxx2-xx2)*(centerxx2-xx2);
      double distance = sqrt( r1sqr + r2sqr );
      double d3 = sqrt((centerxx1*centerxx1)+(centerxx2*centerxx2))+_PMTs_msize[mcase]/2;

      if( distance<(_PMTs_msize[mcase]/2) && (d3<_Rh || mcase==0 || _opt==1) && rd<(_QE_PMT[mcase]/100.) ){
	washit=1;
	PSnum=(int)(wxx1 + wxx2*100) + mcase*10000;
      }
    }

    if(washit==0 && LAPPDs==1){

      if((int)_LAPPDs_per_row[mcase]%2==0) npparityrow=0;
      else npparityrow=1;

      if((int)_LAPPDs_per_column[mcase]%2==0) npparitycolumn=0;
      else npparitycolumn=1;

      double LAPPDs = _LAPPDs_per_row[mcase]*_LAPPDs_per_column[mcase];

      if( _opt==0 && mcase==0 ){
        lrow = (2*_zh)/_LAPPDs_per_row[mcase];
        lcolumn = (2*3.1416*_Rh)/_LAPPDs_per_column[mcase];
      }
      if( _opt==0 && (mcase==5 || mcase==6) ) {
        lrow = (2*_Rh)/_LAPPDs_per_row[mcase];
        lcolumn = (2*_Rh)/_LAPPDs_per_column[mcase];
      }
      if(_opt==1){
        lrow = (2*_xh)/_LAPPDs_per_row[mcase];
        lcolumn = (2*_xh)/_LAPPDs_per_column[mcase];
      }

      if(npparityrow==0){
        wxx1 = (double) round(xx1/lrow);
        centerxx1 = (lrow*wxx1);
      }

      if(npparityrow==1){
        wxx1 = (double) floor(xx1/lrow);
        centerxx1 = (lrow*wxx1) + lrow/2.0;
      }

      if(npparitycolumn==0){
        wxx2 = (double) round(xx2/lcolumn);
        centerxx2 = (lcolumn*wxx2);
      }

      if(npparitycolumn==1){
        wxx2 = (double) floor(xx2/lcolumn);
        centerxx2 = (lcolumn*wxx2) + lcolumn/2.0;
      }
 
      isX1range=false;
      isX2range=false;
      double rd = (rand()/(double)RAND_MAX);

      if(_LAPPDs_mshape[mcase]==1)
      {
        double d1 = fabs(centerxx1-xx1);
        double d2 = fabs(centerxx2-xx2);
        double d3 = (fabs(centerxx1)+_LAPPDs_msize[mcase]/2)*(fabs(centerxx1)
			+_LAPPDs_msize[mcase]/2) + (fabs(centerxx2)+_LAPPDs_msize[mcase]/2)
			*(fabs(centerxx2)+_LAPPDs_msize[mcase]/2);
        if( d1<(_LAPPDs_msize[mcase]/2) && d2<(_LAPPDs_msize[mcase]/2) && fabs(fmod(wxx1+wxx2,2))==1
		&& ( (d3<(_Rh*_Rh) && (mcase==5 || mcase==6) && _opt==0) || (fabs(centerxx1)+lrow/2<_zh 
		&& mcase==0) || (fabs(centerxx1)+lrow/2<_xh && fabs(centerxx2)+lcolumn/2<_xh && _opt==1) ) 
		&& rd<(_QE_LAPPD[mcase]/100.) ){
          washit=2;
          PSnum=(int)(wxx1 + wxx2*100) + mcase*10000;
        }
      }

      if(_LAPPDs_mshape[mcase]==0)
      {
        double r1sqr = (centerxx1-xx1)*(centerxx1-xx1);
        double r2sqr = (centerxx2-xx2)*(centerxx2-xx2);
        double distance = sqrt( r1sqr + r2sqr );
        double d3 = sqrt((centerxx1*centerxx1)+(centerxx2*centerxx2))+_LAPPDs_msize[mcase]/2;

        if(distance<(_LAPPDs_msize[mcase]/2) && fabs(fmod(wxx1+wxx2,2))==1 
		&& ( (d3<(_Rh*_Rh) && (mcase==5 || mcase==6) && _opt==0) || (fabs(centerxx1)+lrow/2<_zh 
		&& mcase==0) || (fabs(centerxx1)+lrow/2<_xh && fabs(centerxx2)+lcolumn/2<_xh && _opt==1) ) 
		&& rd<(_QE_LAPPD[mcase]/100.) ){
          washit=2;
          PSnum=(int)(wxx1 + wxx2*100) + mcase*10000;
        }
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



