#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMath.h"
#include "RasterizeCircle.hh"
#include <iostream>
#include <fstream>
using namespace std;


ClassImp(RasterizeCircle)

  RasterizeCircle::RasterizeCircle(Int_t nbx, Double_t lx, Double_t hx, Int_t nby, Double_t ly, Double_t hy,Int_t nbz, Double_t lz, Double_t hz)
{
  
  lr=lx;
  hr=hx;
  nbins=nbx;
  binsize = (hr-lr)/((Double_t)nbins);

  binfinder = new TH1D("fb","fb",nbins,lr,hr);
  circlehist = new TH3D("tt","tt",nbx,lx,hx,nby,ly,hy,nbz,lz,hz);
}


Int_t RasterizeCircle::Reset(){

  circlehist->Reset();

  return 1;
}


TH3D* RasterizeCircle::returnhisto(){

  return circlehist;

}


RasterizeCircle::~RasterizeCircle()
{
  delete binfinder;
  delete circlehist;
}



Int_t RasterizeCircle::Rasterize3D(vector<Double_t> hvect, Double_t thealpha, Double_t theS1){

  return (this->Rasterize3D(hvect, thealpha, theS1, 1));
}



Int_t RasterizeCircle::Rasterize3D(vector<Double_t> hvect, Double_t thealpha, Double_t theS1, Double_t wt){

  Double_t theR = theS1*sin(thealpha);

  //  cout<<"Rasterizing: "<<theR<<endl;

  Double_t rhv0=0; Double_t rhv1=0; Double_t rhv2=0;
  Int_t mfastb=-1; Int_t mmedb=-1; Int_t mslowb=-1;

  if( (hvect.at(0)==0) && (hvect.at(1)==0) && (hvect.at(2)==0) ){
    cout<<"ERROR! At least one coordinate of hvect must be non-zero"<<endl;
    return 1;
  }

  if( ((hvect.at(0)!=0) && (hvect.at(1)==0) && (hvect.at(2)==0)) || ((hvect.at(0)==0) && (hvect.at(1)!=0) && (hvect.at(2)==0))
      || ((hvect.at(0)==0) && (hvect.at(1)==0) && (hvect.at(2)!=0)) ){

    this->Rasterize3Dflat(theR, hvect.at(0), hvect.at(1), hvect.at(2), wt);
    return 1;
  } else{

    for(Int_t i=0; i<3; i++){

      if( fabs(hvect.at(i))>fabs(rhv2) ){
	
	rhv0=rhv1;
	rhv1=rhv2;
	rhv2=hvect.at(i);
	mslowb=mmedb;
	mmedb=mfastb;
	mfastb=i;
      } else if ( fabs(hvect.at(i))>fabs(rhv1) ){

	rhv0=rhv1;
	rhv1=hvect.at(i);
	mslowb=mmedb;
	mmedb=i;
      } else if( fabs(hvect.at(i))>fabs(rhv0) ){
	rhv0=hvect.at(i);
	mslowb=i;
      }
    }

    TVector3 rankedhvect(rhv0,rhv1,rhv2);
    this->Rasterize3Dspeedorder(rankedhvect, thealpha, theS1, mfastb, mmedb, mslowb, wt);

    return 1;
  }

  return 1;
}




Int_t RasterizeCircle::Rasterize3Dspeedorder(TVector3 hvect, Double_t thealpha, Double_t theS1, Int_t tfastb, Int_t tmedb, Int_t tslowb, Double_t wt)
{
  fastb = tfastb;
  medb = tmedb;
  slowb = tslowb;

  Double_t theR = theS1*sin(thealpha);
  Double_t theC = theS1*cos(thealpha);
  Double_t pmag = hvect.Mag();


  // Calculate ranges for the looping
  TVector3 cVect = (theC/pmag)*hvect;
  Double_t cV0 = cVect.X();   
  Double_t cV1 = cVect.Y();  
  Double_t cV2 = cVect.Z();
  TVector3 v0axis(1.,0.,0.);  
  TVector3 v1axis(0.,1.,0.);  
  TVector3 v2axis(0.,0.,1.);
  TVector3 cVectXYplane(cVect.X(),cVect.Y(),0);
  TVector3 cVectXZplane(cVect.X(),0,cVect.Z());
  TVector3 cVectYZplane(0,cVect.Y(),cVect.Z());

  Double_t dBeta0 = cVect.Angle(v0axis);  
  Double_t dBeta1 = cVect.Angle(v1axis);  
  Double_t dBeta2 = cVect.Angle(v2axis);  

  //  cout<<"PPPOPPOPOPO "<<" theC="<<theC<<" pmag="<<pmag<<" cV0="<<cV0<<" R="<<theR<<" XYmag="<<(cVectXYplane.Mag())<<endl;

  Double_t d0max = cV0 + theR*fabs(sin( dBeta0 ));
  Double_t d0min = cV0 - theR*fabs(sin( dBeta0 ));
  Double_t d1max = cV1 + theR*fabs(sin( dBeta1 ));
  Double_t d1min = cV1 - theR*fabs(sin( dBeta1 ));
  Double_t d2max = cV2 + theR*fabs(sin( dBeta2 ));
  Double_t d2min = cV2 - theR*fabs(sin( dBeta2 ));
  
  Int_t d0minBin = binfinder->FindBin(d0min);
  Int_t d1minBin = binfinder->FindBin(d1min);
  Int_t d2minBin = binfinder->FindBin(d2min);
  Int_t d0maxBin = binfinder->FindBin(d0max);
  Int_t d1maxBin = binfinder->FindBin(d1max);
  Int_t d2maxBin = binfinder->FindBin(d2max);

  Int_t n0steps = d0maxBin - d0minBin;  

  //  cout<<"mins and maxes "<<d0min<<" "<<d0max<<" "<<d1min<<" "<<d1max<<" "<<d2min<<" "<<d2max<<endl;
  //  cout<<"MAX and MIN BINS: "<<d0minBin<<" "<<d0maxBin<<" "<<d1minBin<<" "<<d1maxBin<<" "<<d2minBin<<" "<<d2maxBin<<" "<<dBeta1<<" "<<dBeta2<<endl;

  //calculate corresponding y and z for this particular x value
  Double_t init0 = d0min;
  Double_t init1 = ((2*cV1*(cV0*cV0+cV1*cV1+cV2*cV2-cV0*init0))/(cV2*cV2))/(2*((1 + (cV1*cV1)/(cV2*cV2))));
  Double_t init2 = (cV0*cV0+cV1*cV1+cV2*cV2-cV0*init0-cV1*init1)/cV2;

  Double_t final1= ((2*cV1*(cV0*cV0+cV1*cV1+cV2*cV2-cV0*d0max))/(cV2*cV2))/(2*((1 + (cV1*cV1)/(cV2*cV2))));
  Double_t final2= (cV0*cV0+cV1*cV1+cV2*cV2-cV0*d0max-cV1*final1)/cV2;

  Int_t init0bin = binfinder->FindBin(init0);
  Int_t init1bin = binfinder->FindBin(init1);
  Int_t init2bin = binfinder->FindBin(init2);
  Int_t final1bin = binfinder->FindBin(final1);
  Int_t final2bin = binfinder->FindBin(final2);

  //  cout<<"d0min  "<<d0min<<" init0 "<<init0<<" init1 "<<init1<<" init2 "<<init2<<" sqrt "<< sqrt((init0-cV0)*(init0-cV0) + (init1-cV1)*(init1-cV1) + (init2-cV2)*(init2-cV2)) <<" "<<theR<<endl;
  // Initialize loop variables

  Double_t init0bc,init1bc,init2bc;
  Double_t new1pos,new2pos,new0bc;
  Double_t dummy1,dummy2;

  //Looping /////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  for(Int_t i=0; i<n0steps+1; i++){

    init0bc = binfinder->GetBinCenter(init0bin);

    //    cout<<i<<" "<<init0bc<<" "<<init0bin<<endl;
      
    if( (i>0) && (i<n0steps) ){
    // calculate corresponding y and z, ybin, zbin for this particular x value...
      this->SolveD1D2(theR, cV0, cV1, cV2, init0bc, dummy1, dummy2, new1pos, new2pos);
      init1bin = binfinder->FindBin(new1pos);
      init2bin = binfinder->FindBin(new2pos);
    }
    if(i==(n0steps)){
      init1bin = final1bin;
      init2bin = final2bin;
      //      cout<<"HERE "<<init0bin<<" "<<init1bin<<" "<<init2bin<<endl;
    }

    init1bc = binfinder->GetBinCenter(init1bin);
    init2bc = binfinder->GetBinCenter(init2bin);
    
    this->InnerLoops(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, -1, wt);
    this->InnerLoops(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, 1, wt);

    if( (i>0) && (i<n0steps) ){
    // calculate corresponding y and z, ybin, zbin for this particular x value...
      this->SolveD1D2(theR, cV0, cV1, cV2, init0bc, new1pos, new2pos, dummy1, dummy2);
      init1bin = binfinder->FindBin(new1pos);
      init2bin = binfinder->FindBin(new2pos);
      this->InnerLoops(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, 1, wt);
      this->InnerLoops(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, -1, wt);
    }
    init0bin++;
  }

  // Done Looping /////////////////////////////////////////////////////////////////////////////////  

  //  cout<<"DONE ALTOGETHER"<<endl;
  return 1;

}



Int_t RasterizeCircle::InnerLoops(Int_t init0bin, Int_t init1bin, Int_t init2bin, Int_t d1minBin, Int_t d2minBin, Int_t d1maxBin, Int_t d2maxBin, Double_t cV0, Double_t cV1, Double_t cV2, Double_t theR, Int_t incr1, Double_t wt)
{

  bool binchange0=false;    bool binchange1=false;    bool binchange2=false;


  Double_t new0pos,new1pos,new2pos;
  Double_t new1bc,new2bc;

  Double_t init0bc = binfinder->GetBinCenter(init0bin);
  Double_t init1bc = binfinder->GetBinCenter(init1bin);
  Double_t init2bc = binfinder->GetBinCenter(init2bin);
  
  while( !binchange0 && !binchange1 && (init1bin<=d1maxBin) && (init1bin>=d1minBin) && (init2bin<=d2maxBin) && (init2bin>=d2minBin) ){
    
    binchange1=false;
    binchange2=false;

    this->InnerMostLoop(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, 1, wt);
    this->InnerMostLoop(init0bin, init1bin, init2bin, d1minBin, d2minBin, d1maxBin, d2maxBin, cV0, cV1, cV2, theR, -1, wt);

    binchange0=false;    binchange1=false;    binchange2=false;
    
    //increment d1-bin in gradient direction...check if circle still passes through the new bin...
    // if so, repeat the d1 loop, loop and fill d2-bins for new d1 bin

    init1bin = init1bin + incr1;
    new1bc = binfinder->GetBinCenter(init1bin);

    this->SolveForFixedCoordinate(theR, cV1, cV0, cV2, new1bc, init0bc, init2bc, new0pos, new2pos);
    
    if( fabs(new0pos-init0bc) > (binsize/2.) ) binchange0 = true;
    //    if( fabs(new2pos-new2bc) > (binsize/2.) ) binchange2 = true;
    
    //if circle does not pass through the new bin, check if an additional step of the the d2 bin works...
    // if not we're done looping in d1 an d2
    if( binchange0 || binchange1 ){
      
      //      cout<<"blrrg "<<binchange0<<" "<<binchange1<<endl;
      
      init2bin = init2bin + 1; //this is an issue
      new2bc = binfinder->GetBinCenter(init2bin);
      this->SolveForFixedCoordinate(theR, cV2, cV0, cV1, new2bc, init0bc, init1bc, new0pos, new1pos);       

      if( fabs(new0pos-init0bc) < (binsize/2.) ) binchange0 = false;
      if( fabs(new1pos-init1bc) < (binsize/2.) ) binchange1 = false;
      //      cout<<init2bin<<" W00t "<<binchange0<<" "<<binchange1<<endl;
    }
  }

  return 1;
}




Int_t RasterizeCircle::InnerMostLoop(Int_t init0bin, Int_t init1bin, Int_t init2bin, Int_t d1minBin, Int_t d2minBin, Int_t d1maxBin, Int_t d2maxBin, Double_t cV0, Double_t cV1, Double_t cV2, Double_t theR, Int_t incr2, Double_t wt){

  Double_t init0bc = binfinder->GetBinCenter(init0bin);
  Double_t init1bc = binfinder->GetBinCenter(init1bin);

  bool binchange0= false;
  bool binchange1=false;
  Int_t new2bin=init2bin;
  Double_t new2bc;
  Double_t new0pos,new1pos;
    
  while( !binchange1 && !binchange0 && (new2bin<=d2maxBin) && (new2bin>=d2minBin)  ){
      
      // Filling Histo Here
      
    //      cout<<"bins "<<init0bin<<" "<<init1bin<<" "<<new2bin<<endl;
    this->FillingAction(init0bin, init1bin, new2bin, wt);
      
      new2bin = new2bin + incr2;
      new2bc = binfinder->GetBinCenter(new2bin);
      this->SolveForFixedCoordinate(theR, cV2, cV0, cV1, new2bc, init0bc, init1bc, new0pos, new1pos);

	//      this->FindNearestCirclePoint(init0bc,init1bc,new2bc,cV0,cV1,cV2,theR,new0pos,new1pos,new2pos);
      
      if( fabs(new0pos-init0bc) > (binsize/2.) ) binchange0 = true;
      if( fabs(new1pos-init1bc) > (binsize/2.) ) binchange1 = true;
    }

  return 1;
}




Int_t  RasterizeCircle::FillingAction(Int_t bin0, Int_t bin1, Int_t bin2, Double_t wt)
{

  Int_t xbin,ybin,zbin;
  
  if(fastb==0){
    if(medb==1){
      xbin=bin2; ybin=bin1; zbin=bin0;	  
    } else{
      xbin=bin2; ybin=bin0; zbin=bin1;
    }
  }
  
  if(fastb==1){
    if(medb==0){
      ybin=bin2; xbin=bin1; zbin=bin0;
    } else{
      ybin=bin2; xbin=bin0; zbin=bin1;
    }
  }
  
  if(fastb==2){
    if(medb==0){
      zbin=bin2; ybin=bin1; xbin=bin0;
    } else{
      zbin=bin2; ybin=bin0; xbin=bin1;
    }
  }
  
  Double_t bc = circlehist->GetBinContent(xbin,ybin,zbin);
  circlehist->SetBinContent(xbin,ybin,zbin,(bc+wt));

  return 1;
}



Int_t RasterizeCircle::SolveD1D2(Double_t theR, Double_t cV0, Double_t cV1, Double_t cV2, Double_t d0coor, Double_t &hid1, Double_t &hid2, Double_t &lowd1, Double_t &lowd2)
{
  
  Double_t cterm = (cV0*cV0 + cV1*cV1 + cV2*cV2 - cV0*d0coor);
  Double_t theC = -(theR*theR - cV0*cV0 - cV1*cV1 - cV2*cV2 -(d0coor*d0coor - 2*d0coor*cV0) + 2*cterm - ((cterm*cterm)/(cV2*cV2)));
  Double_t theA = (1 + (cV1*cV1)/(cV2*cV2));
  Double_t theB = -(2*cV1*cterm)/(cV2*cV2);
  
  hid1=-55555;  lowd1=-5555;  lowd2=-55555;  hid2=-5555;  

  if( (theB*theB - 4*theA*theC) >0 ){
    hid1 = (-theB + sqrt(theB*theB - 4*theA*theC))/(2*theA);
    lowd1 = (-theB - sqrt(theB*theB - 4*theA*theC))/(2*theA);
  }
  else{
    //    std::cout<<"AAAARRGGHH!!!!  B2-4AC="<<(theB*theB - 4*theA*theC)<<" B="<<theB<<" A="<<theA<<" C="<<theC<<" d0coor="<<d0coor<<std::endl;
  }
  
  if( (theB*theB - 4*theA*theC) >0 ){
    hid2= (cterm-cV1*hid1)/cV2;
    lowd2= (cterm-cV1*lowd1)/cV2;
  }

  return 1;
}


Int_t RasterizeCircle::SolveForFixedCoordinate(Double_t theR, Double_t cVfixed, Double_t cVA, Double_t cVB, Double_t fixedcoor, Double_t coorA, Double_t coorB, Double_t &theAcoor, Double_t &theBcoor)
{
  
  Double_t cterm = (cVfixed*cVfixed + cVA*cVA + cVB*cVB - cVfixed*fixedcoor);
  Double_t theC = -(theR*theR - cVfixed*cVfixed - cVA*cVA - cVB*cVB -(fixedcoor*fixedcoor - 2*fixedcoor*cVfixed) + 2*cterm - ((cterm*cterm)/(cVB*cVB)));
  Double_t theA = (1 + (cVA*cVA)/(cVB*cVB));
  Double_t theB = -(2*cVA*cterm)/(cVB*cVB);
  
  Double_t hiA=-55555;  Double_t lowA=-5555;   Double_t lowB=-55555;   Double_t hiB=-5555;
 
  
  if( (theB*theB - 4*theA*theC) >0 ){
    hiA = (-theB + sqrt(theB*theB - 4*theA*theC))/(2*theA);
    lowA = (-theB - sqrt(theB*theB - 4*theA*theC))/(2*theA);
  }
  else{
    //    std::cout<<"AAAARRGGHH!!!!  B2-4AC="<<(theB*theB - 4*theA*theC)<<" B="<<theB<<" A="<<theA<<" C="<<theC<<" fixedcoor="<<fixedcoor<<std::endl;
  }
  
  if( (theB*theB - 4*theA*theC) >0 ){
    hiB= (cterm-cVA*hiA)/cVB;
    lowB= (cterm-cVA*lowA)/cVB;
  }

  if( sqrt( (hiA-coorA)*(hiA-coorA) +  (hiB-coorB)*(hiB-coorB) ) < sqrt( (lowA-coorA)*(lowA-coorA) +  (lowB-coorB)*(lowB-coorB) ) ){
    theAcoor=hiA;
    theBcoor=hiB;
  } else{
    theAcoor=lowA;
    theBcoor=lowB;
  }


  return 1;
}




Int_t RasterizeCircle::Rasterize3Dflat(Double_t tradius, Double_t tn0, Double_t tn1, Double_t tn2, Double_t wt)
{

  // Haven't finished writing this class...It currently does nothing.

  Int_t radius = (int)(tradius/binsize);
  Int_t n0 = (int)(tn0/binsize);
  Int_t n1 = (int)(tn1/binsize);
  Int_t n2 = (int)(tn2/binsize);

  Int_t f = 1 - radius;
  Int_t ddF_x = 1;
  Int_t ddF_y = -2 * radius;
  Int_t x = 0;
  Int_t y = radius;
 
  cout<<"i am here"<<endl;
  /*
  circlehist->SetBinContent(n0, n1 + radius,n2,1);
  circlehist->SetBinContent(n0, n1 - radius,n2,1);
  circlehist->SetBinContent(n0 + radius,n1,n2,1);
  circlehist->SetBinContent(n0 - radius,n1,n2,1);
  */
  cout<<"now here"<<endl;

  while(x <= y)
    {
      // ddF_x == 2 * x + 1;
      // ddF_y == -2 * y;
      // f == x*x + y*y - radius*radius + 2*x - y + 1;
      if(f >= 0) 
	{
	  y--;
	  ddF_y += 2;
	  f += ddF_y;
	}
      x++;
      ddF_x += 2;
      f += ddF_x; 

      cout<<x<<" "<<y<<endl;
      /*
      circlehist->SetBinContent(n0 + x, n1 + y,n2,1);
      circlehist->SetBinContent(n0 - x, n1 + y,n2,1);
      circlehist->SetBinContent(n0 + x, n1 - y,n2,1);
      circlehist->SetBinContent(n0 - x, n1 - y,n2,1);
      
      circlehist->SetBinContent(n0 + y, n1 + x,n2,1);
      circlehist->SetBinContent(n0 - y, n1 + x,n2,1);
      circlehist->SetBinContent(n0 + y, n1 - x,n2,1);
      circlehist->SetBinContent(n0 - y, n1 - x,n2,1);
      */
    }

  cout<<"done"<<endl;

  return 1;
}





TH2D* RasterizeCircle::Rasterize2D(Int_t radius, Int_t x0, Int_t y0)
{
  

  TH2D* circle2hist = new TH2D("tt","tt",100,-100.5,100.5,100,-100.5,100.5);

  Int_t f = 1 - radius;
  Int_t ddF_x = 1;
  Int_t ddF_y = -2 * radius;
  Int_t x = 0;
  Int_t y = radius;
 
  cout<<"i am here"<<endl;

  circlehist->SetBinContent(x0, y0 + radius,1);
  circlehist->SetBinContent(x0, y0 - radius,1);
  circlehist->SetBinContent(x0 + radius,y0,1);
  circlehist->SetBinContent(x0 - radius,y0,1);
 
  cout<<"now here"<<endl;

  while(x <= y)
    {
      // ddF_x == 2 * x + 1;
      // ddF_y == -2 * y;
      // f == x*x + y*y - radius*radius + 2*x - y + 1;
      if(f >= 0) 
	{
	  y--;
	  ddF_y += 2;
	  f += ddF_y;
	}
      x++;
      ddF_x += 2;
      f += ddF_x; 

      cout<<x<<" "<<y<<endl;

      circlehist->SetBinContent(x0 + x, y0 + y,1);
      circlehist->SetBinContent(x0 - x, y0 + y,1);
      circlehist->SetBinContent(x0 + x, y0 - y,1);
      circlehist->SetBinContent(x0 - x, y0 - y,1);
      circlehist->SetBinContent(x0 + y, y0 + x,1);
      circlehist->SetBinContent(x0 - y, y0 + x,1);
      circlehist->SetBinContent(x0 + y, y0 - x,1);
      circlehist->SetBinContent(x0 - y, y0 - x,1);
    }
  
  cout<<"done"<<endl;

  return circle2hist;
}






