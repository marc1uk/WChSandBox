#include "WCSimIsoChronTransform.hh"
#include "TObject.h"
#include "TMath.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include "TGraph.h"
#include "WCSimRecoCluster.hh"
#include "WCSimWaterModel.hh"
#include "TVector3.h"
#include <vector>
#include "RasterizeCircle.hh"

using namespace std;

ClassImp(WCSimIsoChronTransform);

WCSimIsoChronTransform::WCSimIsoChronTransform(Double_t nBinsV, Double_t binscale){

  _n=1.34689;
  _c=299.792458;

  if((int)nBinsV%2==0){

    nbx=(int) nBinsV;
    nby=(int) nBinsV;
    nbz=(int) nBinsV;

    double halfbin = binscale/2.;
    double nbinshalf = nBinsV/2.;
    
    lx = -(nbinshalf*binscale)-halfbin;
    hx = (nbinshalf*binscale)-halfbin;
    ly = -(nbinshalf*binscale)-halfbin;
    hy = (nbinshalf*binscale)-halfbin;
    lz = -(nbinshalf*binscale)-halfbin;
    hz = (nbinshalf*binscale)-halfbin;

    cout<<lx<<" "<<hx<<endl;
  }
  else{

    nbx=(int) nBinsV;
    nby=(int) nBinsV;
    nbz=(int) nBinsV;

    double halfbin = binscale/2.;
    double nbinshalf = (nBinsV-1)/2.;
    
    lx = -(nbinshalf*binscale)-halfbin;
    hx = (nbinshalf*binscale)+halfbin;
    ly = -(nbinshalf*binscale)-halfbin;
    hy = (nbinshalf*binscale)+halfbin;
    lz = -(nbinshalf*binscale)-halfbin;
    hz = (nbinshalf*binscale)+halfbin;

    cout<<lx<<" "<<hx<<endl;
  }
  
   _theIsochronRaster = new RasterizeCircle(nbx,lx,hx,nby,ly,hy,nbz,lz,hz);

  _thexb = new TH1D("xb","xb",nbx,lx,hx);
  _theyb = new TH1D("yb","yb",nby,ly,hy);
  _thezb = new TH1D("zb","zb",nbz,lz,hz);

  _maxAlphadist = new TH1D("maxalpha","maxalpha",200,-3.14159,3.14159);
  _maxThetaCdist = new TH1D("maxthetac","maxthetac",12000,-3.14159,3.14159);
  _bcontentsdist3d = new TH1D("bch3","bch3",50000,-0.5,50000.5);
  _bcontentsdist2d = new TH1D("bch2","bch2",50000,-0.5,50000.5);
  _minS1dist= new TH1D("ms1","ms1",50000,0.,5000.);
  _nDist = new TH1D("ndist","ndist",400,1.2,1.4);

  nb2x=1000;
  nb2y=1000;
  _S1vsAlphadist = new TH2D("s1va","s1va",nb2x,0.,50000.,nb2y,0,3.);
  //  _S1vsAlphadist = new TH2D("s1va","s1va",nb2x,0.,50000.,nb2y,0,1.57);
  _returnS1vsAlpha = new TH2D("Rs1va","Rs1va",nb2x,0.,50000.,nb2y,0,3.);
  //  _returnS1vsAlpha = new TH2D("Rs1va","Rs1va",nb2x,0.,50000.,nb2y,0,1.57);

  theprojxz = new TH2D("sohpxz","sohpxz",nbx,lx,hx,nbz,lz,hz);;
  theprojxy = new TH2D("sohpxy","sohpxy",nbx,lx,hx,nby,ly,hy);;
  theprojyz = new TH2D("sohpyz","sohpyz",nby,ly,hy,nbz,lz,hz);;

}


WCSimIsoChronTransform::WCSimIsoChronTransform(){

  _n=1.34689;
  _c=299.792458;
  
  nbx=110;
  nby=200;
  nbz=200;
  lx=-95.;
  hx=1005.;
  ly=-995.;
  hy=1005.;
  lz=-995.;
  hz=1005.;
  
  _theIsochronRaster = new RasterizeCircle(nbx,lx,hx,nby,ly,hy,nbz,lz,hz);

  _thexb = new TH1D("xb","xb",nbx,lx,hx);
  _theyb = new TH1D("yb","yb",nby,ly,hy);
  _thezb = new TH1D("zb","zb",nbz,lz,hz);

  _maxAlphadist = new TH1D("maxalpha","maxalpha",200,-3.14159,3.14159);
  _maxThetaCdist = new TH1D("maxthetac","maxthetac",12000,-3.14159,3.14159);
  _bcontentsdist3d = new TH1D("bch3","bch3",50000,-0.5,50000.5);
  _bcontentsdist2d = new TH1D("bch2","bch2",50000,-0.5,50000.5);
  _minS1dist= new TH1D("ms1","ms1",50000,0.,5000.);
  _nDist = new TH1D("ndist","ndist",400,1.2,1.4);

  nb2x=1000;
  nb2y=1000;
  _S1vsAlphadist = new TH2D("s1va","s1va",nb2x,0.,50000.,nb2y,0,3.);
  //  _S1vsAlphadist = new TH2D("s1va","s1va",nb2x,0.,5000.,nb2y,0,1.57);
  //  _returnS1vsAlpha = new TH2D("Rs1va","Rs1va",nb2x,0.,50000.,nb2y,0,1.57);
  _returnS1vsAlpha = new TH2D("Rs1va","Rs1va",nb2x,0.,50000.,nb2y,0,3.);

  theprojxz = new TH2D("sohpxz","sohpxz",nbx,lx,hx,nbz,lz,hz);
  theprojxy = new TH2D("sohpxy","sohpxy",nbx,lx,hx,nby,ly,hy);
  theprojyz = new TH2D("sohpyz","sohpyz",nby,ly,hy,nbz,lz,hz);

}

  WCSimIsoChronTransform::~WCSimIsoChronTransform()
{

  //  delete _htrackreco;
  delete _maxAlphadist;
  delete _maxThetaCdist;
  delete _bcontentsdist2d;
  delete _bcontentsdist3d;
  delete _minS1dist;
  delete _S1vsAlphadist;
  delete _returnS1vsAlpha;
  delete _nDist;
  delete _thexb;
  delete _theyb;
  delete _thezb;
  delete theprojxz;
  delete theprojxy;
  delete theprojyz;

}

void WCSimIsoChronTransform::Reset()
{
  _theIsochronRaster->Reset();

  _S1vsAlphadist->Reset();
  _returnS1vsAlpha->Reset();
  theprojxz->Reset();
  theprojxy->Reset();
  theprojyz->Reset();

  _maxAlphadist->Reset(); 
  _maxThetaCdist->Reset();
  _bcontentsdist3d->Reset();
  _bcontentsdist2d->Reset();
  _minS1dist->Reset();
  _nDist->Reset(); 
}



void WCSimIsoChronTransform::ApplyTransform(vector< vector<double> >  thehits, vector<double> hypVtx)
{
  _S1vsAlphadist->Reset();
  _maxAlphadist->Reset();
  _maxThetaCdist->Reset();
  _minS1dist->Reset();

  for(int i=0; i<(int)(thehits.size()); i++){
  //  for(int i=0; i<1500; i++){

    if(i%1000==0) cout<<"reconstructing from hits...hitcount: "<<i<<endl;

    vector<double> hcoo = thehits.at(i);
    this->ProcessHit(hcoo,hypVtx);
  }
}



void WCSimIsoChronTransform::ApplyTransformPMTres(vector< vector<double> >  thehits, vector<double> hypVtx, double theresolution)
{
  _S1vsAlphadist->Reset();
  _maxAlphadist->Reset();
  _maxThetaCdist->Reset();
  _minS1dist->Reset();

  for(int i=0; i<(int)(thehits.size()); i++){
  //  for(int i=0; i<1500; i++){

    if(i%1000==0) cout<<"reconstructing from hits...hitcount: "<<i<<endl;

    vector<double> hcoo = thehits.at(i);
    this->ProcessHitPMTres(hcoo,hypVtx,theresolution);
  }
}




void WCSimIsoChronTransform::ApplyTransformNoRing(vector< vector<double> >  thehits, vector<double> hypVtx)
{
  
  _S1vsAlphadist->Reset();
  /*
  _maxAlphadist->Reset();
  _maxThetaCdist->Reset();
  _minS1dist->Reset();
  */
  for(int i=0; i<(int)(thehits.size()); i++){
  //  for(int i=0; i<1500; i++){

    if(i%100000==0) cout<<"reconstructing from hits...hitcount: "<<i<<endl;

    vector<double> hcoo = thehits.at(i);
    this->ProcessHitNoRing(hcoo,hypVtx);
  }
}



  void WCSimIsoChronTransform::ProcessHit(vector<double> hitcoordinates, vector<double> hypVtx)
{

  TVector3 hitcoor(hitcoordinates.at(0),hitcoordinates.at(1),hitcoordinates.at(2));
  TVector3 hvect((hitcoordinates.at(0)-hypVtx.at(0)),(hitcoordinates.at(1)-hypVtx.at(1)),(hitcoordinates.at(2)-hypVtx.at(2)));

  double dD=hvect.Mag();
  double dT=hitcoordinates.at(3)-hypVtx.at(3);
  double theR; 
  double minalpha=0;
  double maxalpha = this->CalcMaxAlpha(dD,dT,_n);
  double maxemmang = this->CalcMaxEmissionAngle(dD,dT,_n);  
  double minS1val = this->CalcS1(maxalpha,dD,dT,_n); 

  /*
  cout<<"VERTEX: "<<hypVtx.at(0)<<" "<<hypVtx.at(1)<<" "<<hypVtx.at(2)<<" "<<hypVtx.at(3)<<endl;
  cout<<"Hit: "<<hitcoordinates.at(0)<<" "<<hitcoordinates.at(1)<<" "<<hitcoordinates.at(2)<<" "<<hitcoordinates.at(3)<<endl;
  cout<<"Calculated Stuff: "<<dD<<" "<<dT<<" "<<minS1val<<" "<<maxalpha<<endl;
  */

  _maxAlphadist->Fill(maxalpha);
  _maxThetaCdist->Fill(maxemmang);
  _minS1dist->Fill(minS1val);

  _S1vsAlphadist->Fill(minS1val,maxalpha);

  /*
  if(minS1val!=-55555){
    cout<<minS1val<<" "<<maxalpha<<endl;
  }
  */

  double alphastep=(maxalpha-minalpha)/1000.;
  theR=0;

  if(maxalpha>0){
    for(int i=999; i<1000; i++){
 
      double thealpha=minalpha+(i*alphastep);
      double theS1 = this->CalcS1(thealpha,dD,dT,_n);
      
      //      cout<<"alpha and s1: "<<thealpha<<" "<<theS1<<endl;

      this->FillIsoChronRing_raster(hvect, thealpha, theS1, 1.);
    }
  }
}





void WCSimIsoChronTransform::ProcessHitPMTres(vector<double> hitcoordinates, vector<double> hypVtx, double theresolution)
{

  TVector3 hitcoor(hitcoordinates.at(0),hitcoordinates.at(1),hitcoordinates.at(2));
  TVector3 hvect((hitcoordinates.at(0)-hypVtx.at(0)),(hitcoordinates.at(1)-hypVtx.at(1)),(hitcoordinates.at(2)-hypVtx.at(2)));

  double stepsize = theresolution/10.;

  for(int j=0; j<21; j++){

    double smR = -theresolution + j*stepsize;
    double mwt = ( 1/(sqrt(2*3.14159265*theresolution*theresolution)) )* exp(-(smR*smR/(2*theresolution*theresolution)) );

    //    cout<< smR <<" "<< mwt << endl;

    double dD=hvect.Mag();
    double dT=hitcoordinates.at(3)-hypVtx.at(3) + smR;
    double theR; 
    double minalpha=0;
    double maxalpha = this->CalcMaxAlpha(dD,dT,_n);
    double maxemmang = this->CalcMaxEmissionAngle(dD,dT,_n);  
    double minS1val = this->CalcS1(maxalpha,dD,dT,_n); 
    
    _maxAlphadist->Fill(maxalpha);
    _maxThetaCdist->Fill(maxemmang);
    _minS1dist->Fill(minS1val);
    _S1vsAlphadist->Fill(minS1val,maxalpha);
    
    /* if(minS1val!=-55555){
      cout<<minS1val<<" "<<maxalpha<<endl;
      }*/

    double alphastep=(maxalpha-minalpha)/1000.;
    theR=0;
    
    if(maxalpha>0){
      for(int i=999; i<1000; i++){
	
	double thealpha=minalpha+(i*alphastep);
	double theS1 = this->CalcS1(thealpha,dD,dT,_n);
	
	//      cout<<"alpha and s1: "<<thealpha<<" "<<theS1<<endl;
	
	double mwt=1.;
	this->FillIsoChronRing_raster(hvect, thealpha, theS1, mwt);	
      }
    }
  }
}





  void WCSimIsoChronTransform::ProcessHitNoRing(vector<double> hitcoordinates, vector<double> hypVtx)
{
  
  TVector3 hitcoor(hitcoordinates.at(0),hitcoordinates.at(1),hitcoordinates.at(2));
  TVector3 hvect((hitcoordinates.at(0)-hypVtx.at(0)),(hitcoordinates.at(1)-hypVtx.at(1)),(hitcoordinates.at(2)-hypVtx.at(2)));

  double dD=hvect.Mag();
  double dT=hitcoordinates.at(3)-hypVtx.at(3);
  double theR; 
  double minalpha=0;
  double maxalpha = this->CalcMaxAlpha(dD,dT,_n);
  double maxemmang = this->CalcMaxEmissionAngle(dD,dT,_n);  
  double minS1val = this->CalcS1(maxalpha,dD,dT,_n); 

  _maxAlphadist->Fill(maxalpha);
  _maxThetaCdist->Fill(maxemmang);
  _minS1dist->Fill(minS1val);

  //  if(dD<minS1val) cout<<"dD: "<<dD<<" minS1val: "<<minS1val<<endl;

  if( (dD>=minS1val) && (minS1val<2500.) && (minS1val>500.) ){
    _S1vsAlphadist->Fill(minS1val,(cos(maxalpha)-(1/tan(0.72))*sin(maxalpha)));
  }
  /*
  if(minS1val!=-55555){
    cout<<minS1val<<" "<<maxalpha<<endl;
  }
  */
  double alphastep=(maxalpha-minalpha)/1000.;
  theR=0;
}





Int_t WCSimIsoChronTransform::FillIsoChronRing_raster(TVector3 hvect, Double_t thealpha, Double_t theS1, Double_t theweight){

  vector<double> nhvect;
  nhvect.push_back(hvect.X());
  nhvect.push_back(hvect.Y());
  nhvect.push_back(hvect.Z());

  //  cout<<"in the raster step: "<<nhvect.at(0)<<" "<<nhvect.at(1)<<" "<<nhvect.at(2)<<" "<<thealpha<<" "<<theS1<<endl;

  _theIsochronRaster->Rasterize3D(nhvect,thealpha,theS1);

  return 1;
}



Double_t WCSimIsoChronTransform::CalcMaxAlpha(Double_t dD, Double_t dT, Double_t ni){

   double qA=4*(ni*ni*ni*ni)*dD*dD;
   double qB=-8*(ni*ni)*(_c*dT)*dD;
   double qC= (4*dT*dT*_c*_c) - (4*((ni*ni) - 1)*( (ni*ni*dD*dD) - (dT*dT*_c*_c) ));


   double maxa=-555555;
   if( (qB*qB) > (4*qA*qC) ){
     double cosa = (-qB + sqrt(qB*qB - 4*qA*qC))/(2*qA);
     maxa = acos(cosa);

     double cosalphap = cos(maxa + 0.01);
     double qA1=((ni*ni)-1);
     double qB1=2*((_c*dT)-(ni*ni*dD*cosalphap));
     double qC1=(ni*ni*dD*dD)-(_c*_c*dT*dT);

     double cosalpham = cos(maxa - 0.01);
     qA1=((ni*ni)-1);
     qB1=2*((_c*dT)-(ni*ni*dD*cosalpham));
     qC1=(ni*ni*dD*dD)-(_c*_c*dT*dT);

   } else{
     //   cout<<"CalcMaxAlpha::Causally impossible "<<qB*qB<<" "<<(4*qA*qC)<<endl;
   }

   return maxa;
}


Double_t WCSimIsoChronTransform::CalcMaxEmissionAngle(Double_t dD, Double_t dT, Double_t ni){

  double maxalpha=this->CalcMaxAlpha(dD,dT,ni);
  double mins1 = this->CalcS1(maxalpha,dD,dT,ni);
  double ds2x = dD - mins1*cos(maxalpha);
  double ds2y = mins1*sin(maxalpha);
  double ds2 = sqrt( ds2x*ds2x + ds2y*ds2y );
  double sinwangle;
  if(ds2>0.) sinwangle = dD*sin(maxalpha)/ds2;
  double wangle = asin(sinwangle); 
   //  double maxemissangle = (3.14159265358979 - wangle);
  double maxemissangle=-55555;
  if(mins1>0) maxemissangle= wangle;

  //  cout<<sinwangle<<" "<<wangle<<" "<<maxemissangle<<" "<<mins1<<" "<<ds2<<endl;
 
  return maxemissangle;
}

Double_t WCSimIsoChronTransform::CalcAlpha(Double_t s1, Double_t dD, Double_t dT, Double_t ni){

  double s2 = (_c*dT -s1)/ni;
  double cosalpha;
  if(s1!=0) cosalpha = (s1*s1 + dD*dD - s2*s2)/(2*s1*dD);
  double alpha = acos(cosalpha);

  return alpha;
}


Double_t WCSimIsoChronTransform::CalcS1(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
{
   double theS1=-55555.;

   if(alpha>=0){
     double cosalpha = cos(alpha);
     
     double qA=((ni*ni)-1);
     double qB=2*((_c*dT)-(ni*ni*dD*cosalpha));
     double qC=(ni*ni*dD*dD)-(_c*_c*dT*dT);
     double sqrtterm=qB*qB - 4*qA*qC;
     if( fabs(sqrtterm)<0.0001 ) sqrtterm=0;     
     if( sqrtterm>=0 ){    

       theS1 = (-qB + sqrt(sqrtterm))/(2*qA);
     }
     else{ 
       //       cout<<"CalcS1::Causally Impossible "<<(qB*qB)<<" "<<(4*qA*qC)<<" "<<(qB*qB - 4*qA*qC)<<endl; 
     }
   }
   return theS1;
}

Double_t WCSimIsoChronTransform::CalcS2(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
{
  double mins1 = this->CalcS1(alpha,dD,dT,ni);
  double ds2x = dD - mins1*cos(alpha);
  double ds2y = mins1*sin(alpha);
  double theS2 = sqrt( ds2x*ds2x + ds2y*ds2y );

  return theS2;
}


//Double_t WCSimIsoChronTransform::FindIndex(Double_t dD, Double_t dT)
TGraph* WCSimIsoChronTransform::IndexVsS1(Double_t dD, Double_t dT)
{

  double minimalS1=1000000.;
  double minimalIndex=0.;
  double theindex_i=1.30;

  TGraph* s1vsindex = new TGraph(100);

  for(int i=0; i<100; i++)
    {
      double maxalpha_i=this->CalcMaxAlpha(dD,dT, theindex_i);
      double mins1_i=this->CalcS1(maxalpha_i,dD,dT, theindex_i);
      // double eang_i=this->CalcMaxEmissionAngle(dD,dT);
      bool past_thresh=false;
      s1vsindex->SetPoint(i,theindex_i,mins1_i);
      //      cout<<"fitting..."<<theindex_i<<" "<<mins1_i<<" "<<maxalpha_i<<" "<<eang_i<<" "<<acos(1/theindex_i)<<endl;

      theindex_i+=0.001;
    }
 
  return s1vsindex;
}

TH2D* WCSimIsoChronTransform::XYProjection()
{
  return theprojxy;
}

TH2D* WCSimIsoChronTransform::XZProjection(double bclimit)
{
  cout<<"starting XZProjection..."<<endl;

   for(int i=0; i<nbx; i++){
     //     cout<<i<<endl;
     for(int j=0; j<nbz; j++){

       double zsum=0;
       double zmax=0;

       for(int k=0; k<nby; k++){

	 double nbc = (_theIsochronRaster->returnhisto())->GetBinContent((i+1),(k+1),(j+1));
	 //	 cout<<"here now..."<<i<<" "<<j<<" "<<k<<" "<<nbc<<endl;
	 if(nbc>=0) _bcontentsdist3d->Fill(nbc);
	 //	 if(nbc>0) cout<<i<<" "<<j<<" "<<k<<" "<<nbc<<endl;
	 if(nbc>bclimit) zsum+=nbc;
	 if(nbc>zmax) zmax=nbc;
       }

       //       if(zsum==0) zsum=1;
       theprojxz->SetBinContent((i+1),(j+1),zsum);
       //      if(zmax>400) sohp->SetBinContent((i+1),(j+1),zmax);
       //       else sohp->SetBinContent((i+1),(j+1),2);
       //       if(zsum==0) sohp->SetBinContent((i+1),(j+1),2);
     }
   }

  return theprojxz;
}


TH3D* WCSimIsoChronTransform::GetRaw3dIsoChr()
{

  return (_theIsochronRaster->returnhisto());
}


TH3D* WCSimIsoChronTransform::Get3Dhisto(double bclimit)
{
  //  cout<<"HERE I AMMMM" <<_htrackreco->GetEntries()<<endl;
  TH3D* newhist = (TH3D*) (_theIsochronRaster->returnhisto())->Clone("clone3d");
  newhist->Reset();
  
  for(int i=1; i<134; i++){
    for(int j=1; j<134; j++){
      for(int k=1; k<134; k++){

	double thebc = (_theIsochronRaster->returnhisto())->GetBinContent(i,j,k);
	if(thebc>bclimit){
	  newhist->SetBinContent(i,j,k,thebc);
	}


      }
    }
  }


  return newhist;

}

TH2D* WCSimIsoChronTransform::YZProjection()
{
  return theprojyz;
}

TH2D* WCSimIsoChronTransform::XYSlice()
{
  TH2D* theslice;
  return theslice;
}

TH2D* WCSimIsoChronTransform::XZSlice()
{
  TH2D* theslice;
  return theslice;
}

TH2D* WCSimIsoChronTransform::YZSlice()
{
  TH2D* theslice;
  return theslice;
}

TH1D* WCSimIsoChronTransform::getBCDist3D()
{
  return _bcontentsdist3d;
}

TH1D* WCSimIsoChronTransform::getBCDist2D()
{
  return _bcontentsdist2d;
}

TH1D* WCSimIsoChronTransform::getMaxAlphaDist()
{
  return _maxAlphadist;
}

TH1D* WCSimIsoChronTransform::getMaxThetaCDist()
{
  return _maxThetaCdist;
}


TH1D* WCSimIsoChronTransform::getMinS1Dist()
{
  return _minS1dist;
}

TH1D* WCSimIsoChronTransform::getnDist()
{
  return _nDist;
}

TH2D* WCSimIsoChronTransform::getS1vsAlpha(double thresh)
{
  _returnS1vsAlpha->Reset();
  _bcontentsdist2d->Reset();

  cout<<"integral "<<_S1vsAlphadist->Integral()<<endl;

  for(int i=0; i<nb2x; i++){
    for(int j=0; j<nb2y; j++){

      double mbc = _S1vsAlphadist->GetBinContent(i,j);
      if(mbc>thresh){

	_returnS1vsAlpha->SetBinContent(i,j,mbc);
	_bcontentsdist2d->Fill(mbc);
      }
    }
  }


  return _returnS1vsAlpha;
  }

void WCSimIsoChronTransform::SetConstantIndexofRefraction(double newn)
{
  _n=newn;
}

void WCSimIsoChronTransform::SetConstantSpeedofParticle(double newc)
{
  _c=newc;
}

void WCSimIsoChronTransform::AddWaterModel(WCSimWaterModel* wm)
{
  iswatermodel=true;
  _theWM = wm;
}
