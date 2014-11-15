#include "WCSimIsochron.hh"
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
//#include "RasterizeCircle.hh"

using namespace std;

ClassImp(WCSimIsochron);

WCSimIsochron::WCSimIsochron(){

  _n=1.34689;
  _c=299.792458;

}



  WCSimIsochron::~WCSimIsochron()
{



}

void WCSimIsochron::Reset()
{

}



void WCSimIsochron::ApplyTransform(vector< vector<double> >  thehits, vector<double> hypVtx)
{
  for(int i=0; i<(int)(thehits.size()); i++){

    if(i%1000==0) cout<<"reconstructing from hits...hitcount: "<<i<<endl;

    vector<double> hcoo = thehits.at(i);
    this->ProcessHit(hcoo,hypVtx);
  }
}


  void WCSimIsochron::ProcessHit(vector<double> hitcoordinates, vector<double> hypVtx)
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

  double alphastep=(maxalpha-minalpha)/1000.;
  theR=0;

  if(maxalpha>0){
    for(int i=999; i<1000; i++){
 
      double thealpha=minalpha+(i*alphastep);
      double theS1 = this->CalcS1(thealpha,dD,dT,_n);
      
      //      this->FillIsoChronRing_raster(hvect, thealpha, theS1, 1.);
    }
  }
}




Double_t WCSimIsochron::CalcMaxAlpha(Double_t dD, Double_t dT, Double_t ni){

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


Double_t WCSimIsochron::CalcMaxEmissionAngle(Double_t dD, Double_t dT, Double_t ni){

  double maxalpha=this->CalcMaxAlpha(dD,dT,ni);
  double mins1 = this->CalcS1(maxalpha,dD,dT,ni);
  double ds2x = dD - mins1*cos(maxalpha);
  double ds2y = mins1*sin(maxalpha);
  double ds2 = sqrt( ds2x*ds2x + ds2y*ds2y );
  double sinwangle;
  if(ds2>0.) sinwangle = dD*sin(maxalpha)/ds2;
  double wangle = asin(sinwangle); 

  double maxemissangle=-55555;
  if(mins1>0) maxemissangle= wangle;

  return maxemissangle;
}

Double_t WCSimIsochron::CalcAlpha(Double_t s1, Double_t dD, Double_t dT, Double_t ni){

  double s2 = (_c*dT -s1)/ni;
  double cosalpha;
  if(s1!=0) cosalpha = (s1*s1 + dD*dD - s2*s2)/(2*s1*dD);
  double alpha = acos(cosalpha);

  return alpha;
}


Double_t WCSimIsochron::CalcS1(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
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
     }
   }
   return theS1;
}

Double_t WCSimIsochron::CalcS2(Double_t alpha, Double_t dD, Double_t dT, Double_t ni)
{
  double mins1 = this->CalcS1(alpha,dD,dT,ni);
  double ds2x = dD - mins1*cos(alpha);
  double ds2y = mins1*sin(alpha);
  double theS2 = sqrt( ds2x*ds2x + ds2y*ds2y );

  return theS2;
}


TGraph* WCSimIsochron::IndexVsS1(Double_t dD, Double_t dT)
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


void WCSimIsochron::SetConstantIndexofRefraction(double newn)
{
  _n=newn;
}

void WCSimIsochron::SetConstantSpeedofParticle(double newc)
{
  _c=newc;
}

void WCSimIsochron::AddWaterModel(WCSimWaterModel* wm)
{
  iswatermodel=true;
  _theWM = wm;
}
