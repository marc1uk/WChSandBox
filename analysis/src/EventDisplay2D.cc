//////////////////////////////////////////////////////////////////////
///  EventDisplay2D.cpp
///  Sept 4, 2010 Matt Wetstein
///
///  EventDisplay stores a data structure consisting of WC detector hits, as read in from a text file. 
///  These hits can be weighted or unweighted.One can "smear" truth-level events and even apply
///  acceptance cuts, with the option of keeping or dropping the full truth level information.
///  One can also apply quality cuts, with accepted event being rewritten to the original array
///  and rejected events save in a "reject" array. This class also alows one to draw on a variety
///  of summary information, both from the hits themselves, and from extrapolations about the hits
///  for a given track hypothesis.

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMinuit.h>
#include <TStyle.h>
#include <TRandom.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TBox.h>
#include <TEllipse.h>
#include <TMarker.h>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <TApplication.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TVector3.h>

#include "EventDisplay2D.hh"
#include <string.h>


using namespace std;

EventDisplay2D::EventDisplay2D(TString iname)
{
  _isweighted=false;
  _initiatedDetectorGeo=false;
  _initiatedPMTGeo=false;
  _spatialresolving=true;
  _setThePlots=false;
  _openedoutfile=false;
  _settruetrackparams=false;
  _setazimuth=false;

  _iname=iname;

  theHits.clear();
}


EventDisplay2D::~EventDisplay2D()
{

  theHits.clear();
  /*
  delete azimuthvstresid;
  delete azimuthmeanprofile;
  delete azimuthrmsprofile;
  delete azimuthpeakprofile;
  delete azimuthFWHMprofile;
  */
}

void EventDisplay2D::SetHits(vector< vector<double> > ihits)
{
  theHits.clear();

  theHits = ihits;
}

void EventDisplay2D::SetAzimuth(vector< vector<double> > tresid, vector< vector<double> > azimuth, vector< vector<double> > tracktheta){
  _tresid = tresid;
  _azimuth = azimuth;
  _ttheta = tracktheta;
  _setazimuth = true;
}

void EventDisplay2D::SetGeometry(int type, double d1, double d2, double d3, double clearance)
{
  //cylindrical geometry only...right now

  _initiatedDetectorGeo=true;

  _detectortype=type;
  _d1=d1; //diameter, if cylindrical
  _d2=d2;
  _d3=d3;
  _clearance=clearance;
}


void EventDisplay2D::SetPMTGeometry(int pmttype, double diameter, double coverage, bool spatialresolving)
{

  _pmt_type=pmttype;
  
  // if circular pmts
  if(_pmt_type==0){
    _PMTdiameter= diameter;
  }
  else if(_pmt_type==1){
    _PMTdiameter=(sqrt(2)*diameter); //diameter is the diagonal of the sqare for square pmts 
  }
  
  _spatialcoverage=coverage;
  _spatialresolving=spatialresolving;
}

void EventDisplay2D::SetPlots(TH1D* trD,TH1D* tpD,TH1D* ttD)
{
  trackresidualDist=trD;
  trackposDist=tpD;
  trackthetaDist=ttD;

  _setThePlots=true;
}

void EventDisplay2D::SetTrueTrackParams(vector<double> tparams)
{
  truetrack=tparams;
  _settruetrackparams=true;
}


void EventDisplay2D::InitializeStyle(){

  // create new TStyle
  rootStyle = new  TStyle("My Style", "MINOS Style");

  // set the background color to white
  rootStyle->SetFillColor(10);
  rootStyle->SetFrameFillColor(10);
  rootStyle->SetCanvasColor(10);
  rootStyle->SetPadColor(10);
  rootStyle->SetTitleFillColor(0);
  rootStyle->SetStatColor(10);

  // don't put a colored frame around the plots
  rootStyle->SetFrameBorderMode(0);
  rootStyle->SetCanvasBorderMode(0);
  rootStyle->SetPadBorderMode(0);
  rootStyle->SetLegendBorderSize(0);

  // use the primary color palette
  rootStyle->SetPalette(1,0);

// set the default line color for a histogram to be black
  rootStyle->SetHistLineColor(kBlack);
  // set the default line color for a fit function to be red
  rootStyle->SetFuncColor(kRed);
  // make the axis labels black
  rootStyle->SetLabelColor(kBlack,"xyz");
  // set the default title color to be black
  rootStyle->SetTitleColor(kBlack);
  // set the margins
  rootStyle->SetPadBottomMargin(0.18);
  rootStyle->SetPadTopMargin(0.08);
  rootStyle->SetPadRightMargin(0.08);
  rootStyle->SetPadLeftMargin(0.17);
  // set axis label and title text sizes
  rootStyle->SetLabelFont(42,"xyz");
  rootStyle->SetLabelSize(0.06,"xyz");
  rootStyle->SetLabelOffset(0.015,"xyz");
  rootStyle->SetTitleFont(42,"xyz");
  rootStyle->SetTitleSize(0.07,"xyz");
  rootStyle->SetTitleOffset(1.1,"yz");
  rootStyle->SetTitleOffset(1.0,"x");
  rootStyle->SetStatFont(42);
  rootStyle->SetStatFontSize(0.07);
  rootStyle->SetTitleBorderSize(0);
  rootStyle->SetStatBorderSize(0);
  rootStyle->SetTextFont(42);
  // set line widths
  rootStyle->SetFrameLineWidth(2);
  rootStyle->SetFuncWidth(2);
  // set the number of divisions to show
  rootStyle->SetNdivisions(506, "xy");
  // turn off xy grids
  rootStyle->SetPadGridX(0);
  rootStyle->SetPadGridY(0);
  // set the tick mark style
  rootStyle->SetPadTickX(1);
  // turn off stats
  rootStyle->SetOptStat(0);
  // marker settings
  rootStyle->SetMarkerStyle(20);
  rootStyle->SetMarkerSize(1.2);
  rootStyle->SetLineWidth(1);
  // done
  rootStyle->cd();
  gROOT->ForceStyle();
}


TCanvas* EventDisplay2D::MakeCanvas(int ctype, TString namestring, TH2D* &iwcDisplay)
{
  TCanvas* iwcCanvas;
  if(ctype==0) iwcCanvas = this->MakeCanvasCylinder(namestring, iwcDisplay);
  if(ctype==1) iwcCanvas = this->MakeCanvasCube(namestring, iwcDisplay);
  if(ctype==2) iwcCanvas = this->MakeCanvasThetaPhi(namestring, iwcDisplay);
  if(ctype==3) iwcCanvas = this->MakeCanvasTime(namestring, iwcDisplay);


  return iwcCanvas;
}


TCanvas* EventDisplay2D::MakeCanvasCylinder(TString namestring, TH2D* &iwcDisplay){

  int canvasU;
  int canvasV;

  TCanvas* iwcCanvas;

  TString thecanvasname=namestring;
  TString histname=thecanvasname + "_hist";

  if(_initiatedDetectorGeo){

      double fCylRadius = _d1/2;
      double fCylLength = _d2;
      double fCylDiagonal = sqrt( fCylLength*fCylLength + 4.0*fCylRadius*fCylRadius );
      
      // calculate dimesions of event display
      double fU = 2.0*TMath::Pi()*fCylRadius;
      double fV = 2.0*fCylRadius + fCylLength + 2.0*fCylRadius;
      
      // calulate bin width for histogram
      double binsWidth = 0.005*(fU+fV);
      int binsU = (Int_t)(fU/binsWidth);
      int binsV = (Int_t)(fV/binsWidth);
      
      // calculate dimensions of canvas
      double canvasWidth = 800.0;
      double canvasHeight = (fV/fU)*canvasWidth;
      canvasU = (Int_t)(canvasWidth);
      canvasV = (Int_t)(canvasHeight);
      
      // make canvas

      iwcCanvas = new TCanvas(thecanvasname,"Event Display",
			     canvasU, canvasV);

      // make histogram
      iwcDisplay = new TH2D(histname,histname,binsU, -0.5*fU*fScale, 0.5*fU*fScale, binsV, -0.5*fV*fScale, 0.5*fV*fScale);
      iwcDisplay->GetXaxis()->SetTitle("r * #theta (m)");
      iwcDisplay->GetYaxis()->SetTitle("z (m)");  
      
      iwcDisplay->GetXaxis()->SetTitleSize(0.06);
      iwcDisplay->GetYaxis()->SetTitleSize(0.06);  
      
      iwcDisplay->GetXaxis()->SetLabelSize(0.05);
      iwcDisplay->GetYaxis()->SetLabelSize(0.05);
      
      iwcDisplay->GetXaxis()->SetNdivisions(1007);
      iwcDisplay->GetYaxis()->SetNdivisions(1007);
      
      // outline for side of detector
      TBox* wcCylEdgeSide = new TBox( -TMath::Pi()*fCylRadius*fScale, -0.5*fCylLength*fScale,
				      +TMath::Pi()*fCylRadius*fScale, +0.5*fCylLength*fScale );
      wcCylEdgeSide->SetFillStyle(0);
      wcCylEdgeSide->SetLineColor(1);
      wcCylEdgeSide->SetLineWidth(2);
      
      // outline for top face of detector
      TEllipse *wcCylEdgeTop = new TEllipse(0.0, +0.5*fCylLength*fScale+fCylRadius*fScale, fCylRadius*fScale);
      wcCylEdgeTop->SetFillStyle(0);  
      wcCylEdgeTop->SetLineColor(1);
      wcCylEdgeTop->SetLineWidth(2);
      
      // outline for bottom face of detector
      TEllipse* wcCylEdgeBottom = new TEllipse(0.0, -0.5*fCylLength*fScale-fCylRadius*fScale, fCylRadius*fScale);
      wcCylEdgeBottom->SetFillStyle(0);
      wcCylEdgeBottom->SetLineColor(1);
      wcCylEdgeBottom->SetLineWidth(2);  
      
      iwcCanvas->cd();
      iwcDisplay->Draw();
      wcCylEdgeSide->Draw();
      wcCylEdgeTop->Draw();
      wcCylEdgeBottom->Draw();    
      iwcCanvas->Update();
  }
  
}


TCanvas* EventDisplay2D::MakeCanvasCube(TString namestring, TH2D* &iwcDisplay){

  int canvasU;
  int canvasV;

  TCanvas* iwcCanvas;

  TString thecanvasname=namestring;
  TString histname=thecanvasname + "_hist";

    // cubic detector

  if(_initiatedDetectorGeo){
    
      // calculate dimesions of event display
      double fU = 2*_d3 + 2*_d1;
      double fV = _d2 + 2*_d1;
      
      // calulate bin width for histogram
      double binsWidth = 0.005*(fU+fV);
      int binsU = (Int_t)(fU/binsWidth);
      int binsV = (Int_t)(fV/binsWidth);
      
      // calculate dimensions of canvas
      double canvasWidth = 800.0;
      double canvasHeight = (fV/fU)*canvasWidth;
      canvasU = (Int_t)(canvasWidth);
      canvasV = (Int_t)(canvasHeight);
      
      // make canvas
      
      iwcCanvas = new TCanvas(thecanvasname,"Event Display",
			     canvasU, canvasV);
      

      // make histogram
      iwcDisplay = new TH2D(histname,histname,binsU, 0, fU*fScale, binsV, -0.5*fV*fScale, 0.5*fV*fScale);
      
      iwcDisplay->GetXaxis()->SetTitleSize(0.06);
      iwcDisplay->GetYaxis()->SetTitleSize(0.06);  
      
      iwcDisplay->GetXaxis()->SetLabelSize(0.05);
      iwcDisplay->GetYaxis()->SetLabelSize(0.05);
      
      iwcDisplay->GetXaxis()->SetNdivisions(1007);
      iwcDisplay->GetYaxis()->SetNdivisions(1007);

      // outline for back of detector
      TBox* wcBack = new TBox( 0., -0.5*_d2*fScale,
			       _d3*fScale, 0.5*_d2*fScale );
      wcBack->SetFillStyle(0);
      wcBack->SetLineColor(1);
      wcBack->SetLineWidth(2);
      
      // outline for bottom of detector
      TBox* wcBot = new TBox( _d3*fScale, -0.5*_d2*fScale,
			      (_d3+_d1)*fScale, 0.5*_d2*fScale );
      wcBot->SetFillStyle(0);
      wcBot->SetLineColor(1);
      wcBot->SetLineWidth(2);
      
      
      // outline for front of detector
      TBox* wcFront = new TBox( (_d3+_d1)*fScale, -0.5*_d2*fScale,
				(_d3+_d1+_d3)*fScale, 0.5*_d2*fScale );
      wcFront->SetFillStyle(0);
      wcFront->SetLineColor(1);
      wcFront->SetLineWidth(2);
      
      // outline for top of detector
      TBox* wcTop = new TBox( (2*_d3+_d1)*fScale, -0.5*_d2*fScale,
			      (2*_d3+2*_d1)*fScale, 0.5*_d2*fScale );
      wcTop->SetFillStyle(0);
      wcTop->SetLineColor(1);
      wcTop->SetLineWidth(2);
      
      // outline for left of detector
      TBox* wcLeft = new TBox( (_d3+_d1)*fScale, 0.5*_d2*fScale,
			       (2*_d3+_d1)*fScale, (0.5*_d2 + _d1)*fScale );
      wcLeft->SetFillStyle(0);
      wcLeft->SetLineColor(1);
      wcLeft->SetLineWidth(2);
      
      // outline for right of detector
      TBox* wcRight = new TBox( (_d3+_d1)*fScale, -(0.5*_d2 + _d1)*fScale,
				(2*_d3+_d1)*fScale, -0.5*_d2*fScale );
      wcRight->SetFillStyle(0);
      wcRight->SetLineColor(1);
      wcRight->SetLineWidth(2);
      
      iwcCanvas->cd();
      iwcDisplay->Draw();
      wcTop->Draw();
      wcBot->Draw();
      wcLeft->Draw();
      wcRight->Draw();
      wcBack->Draw();
      wcFront->Draw();
            
      iwcCanvas->Update();
   }

 return iwcCanvas;
}





TCanvas* EventDisplay2D::MakeCanvasThetaPhi(TString namestring, TH2D* &iwcDisplay){

  int canvasU;
  int canvasV;

  TCanvas* iwcCanvas;

  TString thecanvasname=namestring;
  TString histname=thecanvasname + "_hist";

  // cylindrical detector
  if(_initiatedDetectorGeo){

    // theta phi plots

    // calculate dimesions of event display
    double fU = 2.0*TMath::Pi();
    double fV = TMath::Pi();
    
    // calulate bin width for histogram
    double binsWidth = 0.005*(fU+fV);
    int binsU = (Int_t)(fU/binsWidth);
    int binsV = (Int_t)(fV/binsWidth);
    
    // calculate dimensions of canvas
    double canvasWidth = 800.0;
    double canvasHeight = (0.5)*canvasWidth;
    canvasU = (Int_t)(canvasWidth);
    canvasV = (Int_t)(canvasHeight);
    
    // make canvas
     
    iwcCanvas = new TCanvas(thecanvasname,"ThetaPhi Display",
				    canvasU, canvasV);
    // make histogram
    iwcDisplay = new TH2D(histname,histname,binsU, -0.5*fU, 0.5*fU, binsV,-0.5*fV,0.5*fV); 
    iwcDisplay->GetXaxis()->SetTitle("phi (rad)");
    iwcDisplay->GetYaxis()->SetTitle("theta (rad)");  
    
    iwcDisplay->GetXaxis()->SetTitleSize(0.06);
    iwcDisplay->GetYaxis()->SetTitleSize(0.06);  
    
    iwcDisplay->GetXaxis()->SetLabelSize(0.05);
    iwcDisplay->GetYaxis()->SetLabelSize(0.05);
    
    iwcCanvas->cd();
    iwcDisplay->Draw();
    iwcCanvas->Update();
  }

  return iwcCanvas;
}




TCanvas* EventDisplay2D::MakeCanvasTime(TString namestring, TH2D* &iwcDisplay){

  int canvasU;
  int canvasV;

  TCanvas* iwcCanvas;

  TString thecanvasname=namestring;
  TString histname=thecanvasname + "_hist";

  if(_initiatedDetectorGeo){

    //arrival time plots
    canvasU = 800;    canvasV = 800;
    iwcCanvas = new TCanvas(thecanvasname,"Time Display",
			    canvasU, canvasV);
  }

 return iwcCanvas;
}


void EventDisplay2D::PlotEventBarrel(string savename, bool savecanvas, double lowt, double hight, bool dowrite){

  InitializeStyle();

  fScale = 0.01;
  TH1D* plainTdist;
  TH1D* NHitdist;
  TH1D** mtdist;
  TCanvas* wcCanvas;
  TH2D* wcDisplay;

  // Make initial Tdist and NHitdist to determine color scales
  plainTdist = new TH1D("plainT","plainT",10000,0.,1000.);
  NHitdist = new TH1D("NHits","NHits",2000,0.,2000.);
  
  for(int i=0; i<theHits.size(); i++){
    double tt = (theHits.at(i)).at(3);
    double tw = (theHits.at(i)).at(4);
    plainTdist->Fill(tt);
    NHitdist->Fill(tw);
  }
  
  double maxT,maxhits,RMST,RMSm;
  // establish color scale
  if(_spatialresolving){
    int mbin = plainTdist->GetMaximumBin();
    maxT = plainTdist->GetBinCenter(mbin);
    RMST = plainTdist->GetRMS();
  }
  else{
    int mbin = NHitdist->GetMaximumBin();
    maxhits = NHitdist->GetBinCenter(mbin);
    RMSm = NHitdist->GetRMS();
  }
  
  double tbin[6]; double wbin[6];
  if(_spatialresolving){
    tbin[0]=maxT-0.5*RMST;
    tbin[1]=maxT;
    tbin[2]=maxT+0.5*RMST;
    tbin[3]=maxT+RMST;
    tbin[4]=maxT+1.5*RMST;
    tbin[5]=maxT+2*RMST;
  }
  else{
    wbin[0]=maxhits-0.5*RMSm;
    wbin[1]=maxhits;
    wbin[2]=maxhits+0.5*RMSm;
    wbin[3]=maxhits+RMSm;
    wbin[4]=maxhits+1.5*RMSm;
    wbin[5]=maxhits+2*RMSm;
  }
  
  
  if(_setazimuth){
    doAzimuthalPlots(10000,12,savename,lowt,hight);
  }
  
  if(_initiatedDetectorGeo){

    double fCylRadius = _d1/2;
    double fCylLength = _d2;
    double fCylDiagonal = sqrt( fCylLength*fCylLength + 4.0*fCylRadius*fCylRadius );
    
    // make canvas
    TString mWCcanvasname;
    mWCcanvasname+="EventView_";
    mWCcanvasname+=savename;

    wcCanvas = MakeCanvas(0,mWCcanvasname,wcDisplay);
    
    double xEnd,yEnd,zEnd,tEnd,wEnd;
    double zbot=100000;
    double ztop=-100000;
    for(int i=0; i<theHits.size(); i++){
      
      xEnd = (theHits.at(i)).at(0); yEnd = (theHits.at(i)).at(1);
      zEnd = (theHits.at(i)).at(2); tEnd = (theHits.at(i)).at(3); wEnd = (theHits.at(i)).at(4);
            
      double phiEnd=0;
      phiEnd = atan2(yEnd, xEnd);
      
      rootStyle->SetPalette(1,0);
      TMarker *m1;
      TBox* mb;
      TEllipse* mc;
      double c1,c2; 

      if(fabs(zEnd)<(0.5*fCylLength)){
	m1 = new TMarker((phiEnd*fCylRadius*fScale),zEnd*fScale,20);
	c1 = (phiEnd*fCylRadius*fScale);
	c2 = zEnd*fScale;
      }
      else if(zEnd==(0.5*fCylLength)){
	m1 = new TMarker(yEnd*fScale,(-fScale*xEnd+fScale*((fCylLength/2.)+fCylRadius)),20);

	c1 = yEnd*fScale;
	c2 = (-fScale*xEnd+fScale*((fCylLength/2.)+fCylRadius));
      } 
      else if(zEnd==-(0.5*fCylLength)){
	m1 = new TMarker(yEnd*fScale,(xEnd*fScale-fScale*((fCylLength/2.)+fCylRadius)),20);	
	c1 = yEnd*fScale;
	c2 = (xEnd*fScale-fScale*((fCylLength/2.)+fCylRadius));
      }

       if(_spatialresolving){
	m1->SetMarkerSize(0.1);
	//m1->SetMarkerSize(0.4);
	if(tEnd<tbin[0]) m1->SetMarkerColor(2);
	if((tEnd>tbin[0])&&(tEnd<tbin[1])) m1->SetMarkerColor(5);
	if((tEnd>tbin[1])&&(tEnd<tbin[2])) m1->SetMarkerColor(7);
	if((tEnd>tbin[2])&&(tEnd<tbin[3])) m1->SetMarkerColor(4);
	if((tEnd>tbin[3])&&(tEnd<tbin[4])) m1->SetMarkerColor(9);
	if((tEnd>tbin[4])&&(tEnd<tbin[5])) m1->SetMarkerColor(12);
	if((tEnd>tbin[5])) m1->SetMarkerColor(14);
      }
      if(!_spatialresolving){
	if(_pmt_type==0){
 	  mc = new TEllipse(c1,c2,0.5*_PMTdiameter*fScale,0.5*_PMTdiameter*fScale);

	  if(wEnd>=wbin[0]) mc->SetFillColor(2);
	  if((wEnd<wbin[0]) && (wEnd>=wbin[1])) mc->SetFillColor(5);
	  if((wEnd<wbin[1]) && (wEnd>=wbin[2])) mc->SetFillColor(7);
	  if((wEnd<wbin[2]) && (wEnd>=wbin[3])) mc->SetFillColor(4);
	  if((wEnd<wbin[3]) && (wEnd>=wbin[4])) mc->SetFillColor(9);
	  if((wEnd<wbin[4]) && (wEnd>=wbin[5])) mc->SetFillColor(12);
	  if(wEnd<wbin[5]) mc->SetFillColor(14);
	}
	else if(_pmt_type==1){
	  double pmtside = fScale*(_PMTdiameter/sqrt(2));
	  mb = new TBox( (c1-0.5*pmtside),(c2-0.5*pmtside),(c1+0.5*pmtside),(c2+0.5*pmtside) );

	  if(wEnd>=wbin[0]) mb->SetFillColor(2);
	  if((wEnd<wbin[0]) && (wEnd>=wbin[1])) mb->SetFillColor(5);
	  if((wEnd<wbin[1]) && (wEnd>=wbin[2])) mb->SetFillColor(7);
	  if((wEnd<wbin[2]) && (wEnd>=wbin[3])) mb->SetFillColor(4);
	  if((wEnd<wbin[3]) && (wEnd>=wbin[4])) mb->SetFillColor(9);
	  if((wEnd<wbin[4]) && (wEnd>=wbin[5])) mb->SetFillColor(12);
	  if(wEnd<wbin[5]) mb->SetFillColor(14);
	}
      }
      
      if( (tEnd>lowt) && (tEnd<hight) ){
	if(_spatialresolving) m1->Draw();
	else if( (!_spatialresolving) && (_pmt_type==0) ) mc->Draw();
	else if( (!_spatialresolving) && (_pmt_type==1) ) mb->Draw();
      }

      double thetEnd=0;
    
    }

  }

  //  Save Relevant Plots

  //  if(!_openedoutfile){


  if(dowrite){
    TString ofname;
    ofname+=_iname;
    ofname+="_plots.root";
  
    TFile *plotfile = new TFile(ofname,"UPDATE");
    _openedoutfile=true;
  
    std::cout<<"about to write"<<std::endl;
     
    if(savecanvas){
     wcCanvas->Write();
    }

    TString ptname;
    ptname+="arrivaltimedist_";
    ptname+=savename;

    if(!_spatialresolving){
      TString nhname;
      nhname+="numPMThits_";
      nhname+=savename;
      NHitdist->Write(nhname);
    }

    plainTdist->Write(ptname);

    std::cout<<"done writing"<<std::endl;

    plotfile->Close();
    std::cout<<"closed"<<std::endl;
    delete wcCanvas;
    delete wcDisplay;
    delete rootStyle;
    delete plainTdist;
    delete NHitdist;

  }

}


void EventDisplay2D::PlotEventThetaPhi(string savename, bool savecanvas, double lowt, double hight, bool dowrite)
{
  // Make initial Tdist and NHitdist to determine color scales
  TH1D* plainTdist = new TH1D("plainT","plainT",10000,0.,1000.);
  TH1D* NHitdist;
  NHitdist = new TH1D("NHits","NHits",2000,0.,2000.);
  TH2D* TPDisplay;

  for(int i=0; i<theHits.size(); i++){
    double tt = (theHits.at(i)).at(3);
    double tw = (theHits.at(i)).at(4);
    plainTdist->Fill(tt);
    NHitdist->Fill(tw);
  }
  
  double maxT,maxhits,RMST,RMSm;
  // establish color scale
  if(_spatialresolving){
    int mbin = plainTdist->GetMaximumBin();
    maxT = plainTdist->GetBinCenter(mbin);
    RMST = plainTdist->GetRMS();
  }
  else{
    int mbin = NHitdist->GetMaximumBin();
    maxhits = NHitdist->GetBinCenter(mbin);
    RMSm = NHitdist->GetRMS();
  }
  
  double tbin[6]; double wbin[6];
  if(_spatialresolving){
    tbin[0]=maxT-0.5*RMST;
    tbin[1]=maxT;
    tbin[2]=maxT+0.5*RMST;
    tbin[3]=maxT+RMST;
    tbin[4]=maxT+1.5*RMST;
    tbin[5]=maxT+2*RMST;
  }
  else{
    wbin[0]=maxhits-0.5*RMSm;
    wbin[1]=maxhits;
    wbin[2]=maxhits+0.5*RMSm;
    wbin[3]=maxhits+RMSm;
    wbin[4]=maxhits+1.5*RMSm;
    wbin[5]=maxhits+2*RMSm;
  }

  TString mTPcanvasname;
  mTPcanvasname+="ThetaPhi_";
  mTPcanvasname+=savename;
  
  TCanvas *TPCanvas = MakeCanvasThetaPhi(mTPcanvasname,TPDisplay);

  for(int i=0; i<theHits.size(); i++){
 
    double xEnd = (theHits.at(i)).at(0);
    double yEnd = (theHits.at(i)).at(1);
    double zEnd = (theHits.at(i)).at(2);
    double tEnd = (theHits.at(i)).at(3);
    double wEnd = (theHits.at(i)).at(4); 
    
    //double phiEnd = kinem::phi(xEnd,yEnd);
    //double thetaEnd = kinem::theta(xEnd,yEnd,zEnd);
    
    double phiEnd = (_azimuth.at(i)).at(0);
    double thetaEnd = (_ttheta.at(i)).at(0);
    double p1 = thetaEnd*cos(phiEnd);
    double p2 = thetaEnd*sin(phiEnd);

    //    TMarker *m1 = new TMarker((phiEnd-3.14159265),(thetaEnd-(0.5*3.14159265)),20);
    TMarker *m1 = new TMarker(p1,p2,20);
    m1->SetMarkerSize(0.1);
    if(!_spatialresolving) m1->SetMarkerSize(0.45);
    if(tEnd<tbin[0]) m1->SetMarkerColor(2);
    if((tEnd>tbin[0])&&(tEnd<tbin[1])) m1->SetMarkerColor(5);
    if((tEnd>tbin[1])&&(tEnd<tbin[2])) m1->SetMarkerColor(7);
    if((tEnd>tbin[2])&&(tEnd<tbin[3])) m1->SetMarkerColor(4);
    if((tEnd>tbin[3])&&(tEnd<tbin[4])) m1->SetMarkerColor(9);
    if((tEnd>tbin[4])&&(tEnd<tbin[5])) m1->SetMarkerColor(12);
    if((tEnd>tbin[5])) m1->SetMarkerColor(14);
    
    if(!_spatialresolving){
      if(wEnd>=wbin[0]) m1->SetMarkerColor(2);
      if((wEnd<wbin[0]) && (wEnd>=wbin[1])) m1->SetMarkerColor(5);
      if((wEnd<wbin[1]) && (wEnd>=wbin[2])) m1->SetMarkerColor(7);
      if((wEnd<wbin[2]) && (wEnd>=wbin[3])) m1->SetMarkerColor(4);
      if((wEnd<wbin[3]) && (wEnd>=wbin[4])) m1->SetMarkerColor(9);
      if((wEnd<wbin[4]) && (wEnd>=wbin[5])) m1->SetMarkerColor(12);
      if(wEnd<wbin[5]) m1->SetMarkerColor(14);
    }
    if( (tEnd>lowt) && (tEnd<hight) ){
      m1->Draw();   
    }   
  }
  
  // Draw the true track center
  if(_settruetrackparams){
    double phiT = truetrack.at(5);
    double thetaT = truetrack.at(4);
  }


  //  Save Relevant Plots
  if(dowrite){
	  TString ofname;
	  ofname+=_iname;
	  ofname+="_plots.root";
  
	  TFile *plotfile = new TFile(ofname,"UPDATE");
	  _openedoutfile=true;
  
	  if(savecanvas){
	    TPCanvas->Write();
	  }

	/*
	  if(_setThePlots){
	    TString trkname;
 	   trkname+="trackTimeResid_";
	    trkname+=savename;
	    trackresidualDist->Write(trkname);
    
	    TString trkname2;
	    trkname2+="trackPosition_";
	    trkname2+=savename;
	    trackposDist->Write(trkname2);
    
	    TString trkname3;
	    trkname3+="trackTheta_";
	    trkname3+=savename;
	    trackthetaDist->Write(trkname3);  
	  }

	  if(_setazimuth) {  
	    //  azimuthvstresid->Write();
	    // azimuthvstheta->Write();
    
	    azimuthmeanprofile->Write();
	    azimuthrmsprofile->Write();
	    azimuthpeakprofile->Write();
	    azimuthFWHMprofile->Write();
    
 	 }
	*/

  	std::cout<<"done writing"<<std::endl;
  	plotfile->Close();

	/*
	  if(_setazimuth) {  
	    delete azimuthvstresid;
	    delete azimuthvstheta;
	    delete azimuthmeanprofile;
	    delete azimuthrmsprofile;
	    delete azimuthpeakprofile;
	    delete azimuthFWHMprofile;
	  }
	*/

	  std::cout<<"closed"<<std::endl;

	  //  Delete;
	  std::cout<<"1st deletes"<<std::endl;

	  delete TPCanvas;
	  delete TPDisplay;
	  delete rootStyle;
	  delete plainTdist;
    }
}


void EventDisplay2D::PlotEventTime(string savename, bool savecanvas, double lowt, double hight, bool dowrite)
{
   // Make initial Tdist and NHitdist to determine color scales
  TH1D* plainTdist = new TH1D("plainT","plainT",10000,0.,1000.);
  TH2D* tDisplay;  

  for(int i=0; i<theHits.size(); i++){
    double tt = (theHits.at(i)).at(3);
    plainTdist->Fill(tt);
  }
  
  double maxT,maxhits,RMST,RMSm;
  // establish color scale
  if(_spatialresolving){
    int mbin = plainTdist->GetMaximumBin();
    maxT = plainTdist->GetBinCenter(mbin);
    RMST = plainTdist->GetRMS();
  }
  
  double tbin[6]; double wbin[6];
  tbin[0]=maxT-0.5*RMST;
  tbin[1]=maxT;
  tbin[2]=maxT+0.5*RMST;
  tbin[3]=maxT+RMST;
  tbin[4]=maxT+1.5*RMST;
  tbin[5]=maxT+2*RMST;
  
  TH1D** mtdist = new TH1D*[7];
  for(int i=0; i<7; i++){
    TString hn;
    hn+="mtdist_";
    hn+=i;
    mtdist[i] = new TH1D(hn,hn,10000,0.,1000.);
  }

  for(int i=0; i<theHits.size(); i++){
      
      double xEnd = (theHits.at(i)).at(0); double yEnd = (theHits.at(i)).at(1);
      double zEnd = (theHits.at(i)).at(2); double tEnd = (theHits.at(i)).at(3); 
      double wEnd = (theHits.at(i)).at(4);
      
      if(tEnd<tbin[0]) mtdist[0]->Fill(tEnd);
      if((tEnd>tbin[0])&&(tEnd<tbin[1])) mtdist[1]->Fill(tEnd);
      if((tEnd>tbin[1])&&(tEnd<tbin[2])) mtdist[2]->Fill(tEnd);
      if((tEnd>tbin[2])&&(tEnd<135)) mtdist[3]->Fill(tEnd);
      if((tEnd>135)&&(tEnd<145)) mtdist[4]->Fill(tEnd);
      if((tEnd>145)&&(tEnd<170)) mtdist[5]->Fill(tEnd);
      if((tEnd>170)) mtdist[6]->Fill(tEnd);
  }

  TString mTcanvasname;
  mTcanvasname+="TimeDisplay_";
  mTcanvasname+=savename;
 
  tDisplay = new TH2D("td","td",100,0.,1.,100,0.,1.);
  TCanvas* tc = MakeCanvas(3,mTcanvasname,tDisplay);

  mtdist[0]->SetLineColor(2);
  mtdist[0]->SetFillColor(2);
  mtdist[0]->Draw();

  mtdist[1]->SetLineColor(5);
  mtdist[1]->SetFillColor(5);
  mtdist[1]->Draw("SAME");

  mtdist[2]->SetLineColor(7);
  mtdist[2]->SetFillColor(7);
  mtdist[2]->Draw("SAME");

  mtdist[3]->SetLineColor(4);
  mtdist[3]->SetFillColor(4);
  mtdist[3]->Draw("SAME");

  mtdist[4]->SetLineColor(9);
  mtdist[4]->SetFillColor(9);
  mtdist[4]->Draw("SAME");

  mtdist[5]->SetLineColor(12);
  mtdist[5]->SetFillColor(12);
  mtdist[5]->Draw("SAME");

  mtdist[6]->SetLineColor(14);
  mtdist[6]->SetFillColor(14);
  mtdist[6]->Draw("SAME");

  //  Save Relevant Plots

  TString ofname;
  ofname+=_iname;
  ofname+="_plots.root";
  
  TFile *plotfile = new TFile(ofname,"UPDATE");
  _openedoutfile=true;
  
  if(savecanvas){
    tc->Write();
  }

  delete tc;
  delete tDisplay;

  for(int j=0; j<7; j++){
    delete mtdist[j];
  }

}



void EventDisplay2D::PlotEventCube(string savename, bool savecanvas, double lowt, double hight, bool dowrite)
{

  InitializeStyle();

  fScale = 0.01;
  int canvasU;
  int canvasV;
  TH1D* plainTdist;
  TH1D* NHitdist;
  TH1D** mtdist;
  TCanvas* wcCanvas;
  TH2D* wcDisplay;
  TH2D* TPDisplay;
  TH2D* tDisplay;

  // Make initial Tdist and NHitdist to determine color scales
  
  plainTdist = new TH1D("plainT","plainT",2400,0.,1000.);
  NHitdist = new TH1D("NHits","NHits",2000,0.,2000.);
  
  for(int i=0; i<theHits.size(); i++){
    
    double tt = (theHits.at(i)).at(3);
    double tw = (theHits.at(i)).at(4);
    
    plainTdist->Fill(tt);
    NHitdist->Fill(tw);
  }
  
  double maxT,maxhits,RMST,RMSm;
  // establish color scale
  if(_spatialresolving){
    int mbin = plainTdist->GetMaximumBin();
    maxT = plainTdist->GetBinCenter(mbin);
    RMST = plainTdist->GetRMS();
  }
  else{
    int mbin = NHitdist->GetMaximumBin();
    maxhits = NHitdist->GetBinCenter(mbin);
    RMSm = NHitdist->GetRMS();
  }
  
  double tbin[6]; double wbin[6];
  if(_spatialresolving){
    tbin[0]=maxT-0.5*RMST;
    tbin[1]=maxT;
    tbin[2]=maxT+0.5*RMST;
    tbin[3]=maxT+RMST;
    tbin[4]=maxT+1.5*RMST;
    tbin[5]=maxT+2*RMST;
  }
  else{
    wbin[0]=maxhits-0.5*RMSm;
    wbin[1]=maxhits;
    wbin[2]=maxhits+0.5*RMSm;
    wbin[3]=maxhits+RMSm;
    wbin[4]=maxhits+1.5*RMSm;
    wbin[5]=maxhits+2*RMSm;
  }
  
  
  mtdist = new TH1D*[7];
  for(int i=0; i<7; i++){
    TString hn;
    hn+="mtdist_";
    hn+=i;
    mtdist[i] = new TH1D(hn,hn,10000,0.,1000.);
  }

  if(_setazimuth){
    doAzimuthalPlots(2400,12,savename,lowt,hight);
  }
  
  if(_initiatedDetectorGeo){

    // make canvas
    TString mWCcanvasname;
    mWCcanvasname+="EventView_";
    mWCcanvasname+=savename;
    
    wcCanvas = MakeCanvas(1,mWCcanvasname,wcDisplay);    

    double xEnd,yEnd,zEnd,tEnd,wEnd;
    double zbot=100000;
    double ztop=-100000;
    for(int i=0; i<theHits.size(); i++){
      
      xEnd = (theHits.at(i)).at(0); yEnd = (theHits.at(i)).at(1);
      zEnd = (theHits.at(i)).at(2); tEnd = (theHits.at(i)).at(3); wEnd = (theHits.at(i)).at(4);
      
      if(tEnd<tbin[0]) mtdist[0]->Fill(tEnd);
      if((tEnd>tbin[0])&&(tEnd<tbin[1])) mtdist[1]->Fill(tEnd);
      if((tEnd>tbin[1])&&(tEnd<tbin[2])) mtdist[2]->Fill(tEnd);
      if((tEnd>tbin[2])&&(tEnd<135)) mtdist[3]->Fill(tEnd);
      if((tEnd>135)&&(tEnd<145)) mtdist[4]->Fill(tEnd);
      if((tEnd>145)&&(tEnd<170)) mtdist[5]->Fill(tEnd);
      if((tEnd>170)) mtdist[6]->Fill(tEnd);
      
      double phiEnd=0;
      phiEnd = atan2(yEnd, xEnd);
      
      rootStyle->SetPalette(1,0);
      
      TMarker *m1;
      TBox* mb;
      TEllipse* mc;
      double c1,c2;

      if(xEnd==(_d1/2.)){
	//front
	m1 = new TMarker((1.5*_d3+_d1+zEnd)*fScale,yEnd*fScale,20);
	c1 = (1.5*_d3+_d1+zEnd)*fScale;
	c2 = yEnd*fScale;
      }
      else if(xEnd==(-_d1/2.)){
	//back
	m1 = new TMarker((0.5*_d3-zEnd)*fScale,yEnd*fScale,20);
	c1 = (0.5*_d3-zEnd)*fScale;
	c2 = yEnd*fScale;
      }
      else if(yEnd==(_d2/2.)){
	//left
	m1 = new TMarker((1.5*_d3+_d1+zEnd)*fScale,(0.5*_d2 + 0.5*_d1 - xEnd)*fScale,20);
	c1 = (1.5*_d3+_d1+zEnd)*fScale;
	c2 = (0.5*_d2 + 0.5*_d1 - xEnd)*fScale;
      }
      else if(yEnd==(-_d2/2.)){
	//right
	m1 = new TMarker((1.5*_d3+_d1+zEnd)*fScale,(-0.5*_d2 - 0.5*_d1 + xEnd)*fScale,20);
	c1 = (1.5*_d3+_d1+zEnd)*fScale;
	c2 = (-0.5*_d2 - 0.5*_d1 + xEnd)*fScale;
      }
      else if(zEnd==(_d3/2.)){
	//top
	m1 = new TMarker((-xEnd + 1.5*_d1 + 2*_d3)*fScale,yEnd*fScale,20);
	c1 = (-xEnd + 1.5*_d1 + 2*_d3)*fScale;
	c2 = yEnd*fScale;
      }
      else if(zEnd==(-_d3/2.)){
	//bottom
	m1 = new TMarker((xEnd + 0.5*_d1 + _d3)*fScale,yEnd*fScale,20);
	c1 = (xEnd + 0.5*_d1 + _d3)*fScale;
	c2 = yEnd*fScale;
      }

      if(_spatialresolving){
	m1->SetMarkerSize(0.1);
	if(tEnd<tbin[0]) m1->SetMarkerColor(2);
	if((tEnd>tbin[0])&&(tEnd<tbin[1])) m1->SetMarkerColor(5);
	if((tEnd>tbin[1])&&(tEnd<tbin[2])) m1->SetMarkerColor(7);
	if((tEnd>tbin[2])&&(tEnd<tbin[3])) m1->SetMarkerColor(4);
	if((tEnd>tbin[3])&&(tEnd<tbin[4])) m1->SetMarkerColor(9);
	if((tEnd>tbin[4])&&(tEnd<tbin[5])) m1->SetMarkerColor(12);
	if((tEnd>tbin[5])) m1->SetMarkerColor(14);
      }
      if(!_spatialresolving){
	if(_pmt_type==0){
	  mc = new TEllipse(c1,c2,0.5*_PMTdiameter*fScale,0.5*_PMTdiameter*fScale);

	  if(wEnd>=wbin[0]) mc->SetFillColor(2);
	  if((wEnd<wbin[0]) && (wEnd>=wbin[1])) mc->SetFillColor(5);
	  if((wEnd<wbin[1]) && (wEnd>=wbin[2])) mc->SetFillColor(7);
	  if((wEnd<wbin[2]) && (wEnd>=wbin[3])) mc->SetFillColor(4);
	  if((wEnd<wbin[3]) && (wEnd>=wbin[4])) mc->SetFillColor(9);
	  if((wEnd<wbin[4]) && (wEnd>=wbin[5])) mc->SetFillColor(12);
	  if(wEnd<wbin[5]) mc->SetFillColor(14);
	}
	else if(_pmt_type==1){
	  double pmtside = fScale*(_PMTdiameter/sqrt(2));
	  mb = new TBox( (c1-0.5*pmtside),(c2-0.5*pmtside),(c1+0.5*pmtside),(c2+0.5*pmtside) );

	  if(wEnd>=wbin[0]) mb->SetFillColor(2);
	  if((wEnd<wbin[0]) && (wEnd>=wbin[1])) mb->SetFillColor(5);
	  if((wEnd<wbin[1]) && (wEnd>=wbin[2])) mb->SetFillColor(7);
	  if((wEnd<wbin[2]) && (wEnd>=wbin[3])) mb->SetFillColor(4);
	  if((wEnd<wbin[3]) && (wEnd>=wbin[4])) mb->SetFillColor(9);
	  if((wEnd<wbin[4]) && (wEnd>=wbin[5])) mb->SetFillColor(12);
	  if(wEnd<wbin[5]) mb->SetFillColor(14);
	}
      }
    
      if( (tEnd>lowt) && (tEnd<hight) ){
	if(_spatialresolving) m1->Draw();
	else if( (!_spatialresolving) && (_pmt_type==0) ) mc->Draw();
	else if( (!_spatialresolving) && (_pmt_type==1) ) mb->Draw();
      }
      double thetEnd=0;
     }
   
    //   cout<<ztop<<" "<<zbot<<endl;
  }

  // Plot Theta/Phi

  TString mTPcanvasname;
  mTPcanvasname+="ThetaPhi_";
  mTPcanvasname+=savename;

  TCanvas *TPCanvas = MakeCanvas(2,mTPcanvasname,TPDisplay);

  for(int i=0; i<theHits.size(); i++){
 
    double xEnd = (theHits.at(i)).at(0);
    double yEnd = (theHits.at(i)).at(1);
    double zEnd = (theHits.at(i)).at(2);
    double tEnd = (theHits.at(i)).at(3);
    double wEnd = (theHits.at(i)).at(4); 

    /*
    double phiEnd = kinem::phi(xEnd,yEnd);
    double thetaEnd = kinem::theta(xEnd,yEnd,zEnd);
    */

    double phiEnd = (_azimuth.at(i)).at(0);
    double thetaEnd = (_ttheta.at(i)).at(0);
    double p1 = thetaEnd*cos(phiEnd);
    double p2 = thetaEnd*sin(phiEnd);

    //    TMarker *m1 = new TMarker((phiEnd-3.14159265),(thetaEnd-(0.5*3.14159265)),20);
    TMarker *m1 = new TMarker(p1,p2,20);
    m1->SetMarkerSize(0.1);
    if(!_spatialresolving) m1->SetMarkerSize(0.45);

    if(tEnd<tbin[0]) m1->SetMarkerColor(2);
    if((tEnd>tbin[0])&&(tEnd<tbin[1])) m1->SetMarkerColor(5);
    if((tEnd>tbin[1])&&(tEnd<tbin[2])) m1->SetMarkerColor(7);
    if((tEnd>tbin[2])&&(tEnd<tbin[3])) m1->SetMarkerColor(4);
    if((tEnd>tbin[3])&&(tEnd<tbin[4])) m1->SetMarkerColor(9);
    if((tEnd>tbin[4])&&(tEnd<tbin[5])) m1->SetMarkerColor(12);
    if((tEnd>tbin[5])) m1->SetMarkerColor(14);
    
    if(!_spatialresolving){
      if(wEnd>=wbin[0]) m1->SetMarkerColor(2);
      if((wEnd<wbin[0]) && (wEnd>=wbin[1])) m1->SetMarkerColor(5);
      if((wEnd<wbin[1]) && (wEnd>=wbin[2])) m1->SetMarkerColor(7);
      if((wEnd<wbin[2]) && (wEnd>=wbin[3])) m1->SetMarkerColor(4);
      if((wEnd<wbin[3]) && (wEnd>=wbin[4])) m1->SetMarkerColor(9);
      if((wEnd<wbin[4]) && (wEnd>=wbin[5])) m1->SetMarkerColor(12);
      if(wEnd<wbin[5]) m1->SetMarkerColor(14);
    }
    if( (tEnd>lowt) && (tEnd<hight) ){	
      m1->Draw();    
    }  
  }
  

  // Draw the true track center
  if(_settruetrackparams){

    double phiT = truetrack.at(5);
    double thetaT = truetrack.at(4);

  }


 // Plot Arrival Time Distribution
  TString mTcanvasname;
  mTcanvasname+="TimeDisplay_";
  mTcanvasname+=savename;

  tDisplay = new TH2D("td","td",100,0.,1.,100,0.,1.); 

  TCanvas* tc = MakeCanvas(3,mTcanvasname,tDisplay);

  mtdist[0]->SetLineColor(2);
  mtdist[0]->SetFillColor(2);
  mtdist[0]->Draw();

  mtdist[1]->SetLineColor(5);
  mtdist[1]->SetFillColor(5);
  mtdist[1]->Draw("SAME");

  mtdist[2]->SetLineColor(7);
  mtdist[2]->SetFillColor(7);
  mtdist[2]->Draw("SAME");

  mtdist[3]->SetLineColor(4);
  mtdist[3]->SetFillColor(4);
  mtdist[3]->Draw("SAME");

  mtdist[4]->SetLineColor(9);
  mtdist[4]->SetFillColor(9);
  mtdist[4]->Draw("SAME");

  mtdist[5]->SetLineColor(12);
  mtdist[5]->SetFillColor(12);
  mtdist[5]->Draw("SAME");

  mtdist[6]->SetLineColor(14);
  mtdist[6]->SetFillColor(14);
  mtdist[6]->Draw("SAME");

  
  //  Save Relevant Plots
 
  //  if(!_openedoutfile){
  TString ofname;
  ofname+=_iname;
  ofname+="_plots.root";
  
  TFile *plotfile = new TFile(ofname,"UPDATE");
  _openedoutfile=true;
  //  }
 
  if(savecanvas){
    wcCanvas->Write();
    TPCanvas->Write();
    tc->Write();
  }
  if(_setThePlots){
    TString trkname;
    trkname+="trackTimeResid_";
    trkname+=savename;
    trackresidualDist->Write(trkname);
    
    TString trkname2;
    trkname2+="trackPosition_";
    trkname2+=savename;
    trackposDist->Write(trkname2);
    
    TString trkname3;
    trkname3+="trackTheta_";
    trkname3+=savename;
    trackthetaDist->Write(trkname3);  
  }

  TString ptname;
  ptname+="arrivaltimedist_";
  ptname+=savename;
  plainTdist->Write(ptname);
  if(!_spatialresolving){
    TString nhname;
    nhname+="numPMThits_";
    nhname+=savename;
    NHitdist->Write(nhname);
  }

  if(_setazimuth) {  
    azimuthvstresid->Write();
    azimuthvstheta->Write();
    /*
    azimuthmeanprofile->Write();
    azimuthrmsprofile->Write();
    azimuthpeakprofile->Write();
    azimuthFWHMprofile->Write();
    */
  }

  plotfile->Close();

  if(_setazimuth) {  
    delete azimuthvstresid;
    delete azimuthvstheta;
    delete azimuthmeanprofile;
    delete azimuthrmsprofile;
    delete azimuthpeakprofile;
    delete azimuthFWHMprofile;
  }

  delete wcCanvas;
  delete wcDisplay;
  delete TPDisplay;
  delete tDisplay;
  delete TPCanvas;
  delete tc;
  delete rootStyle;
  for(int j=0; j<7; j++){
    delete mtdist[j];
  }
  
  delete plainTdist;
  delete NHitdist;
 
  //  cout<<"DONE saving."<<endl;
}







void EventDisplay2D::doAzimuthalPlots(int tbins, int phibins, string savename, double lowt, double hight){


    TString ahistname;
    ahistname+="azimuthvstresid_";
    ahistname+= savename;

    TString athistname;
    athistname+="azimuthvstheta_";
    athistname+= savename;

    TString apmhistname;
    apmhistname+="azimuthmeanprofile_";
    apmhistname+= savename;

    TString aprhistname;
    aprhistname+="azimuthrmsprofile_";
    aprhistname+= savename;

    TString apphistname;
    apphistname+="azimuthpeakprofile_";
    apphistname+= savename;

    TString apwhistname;
    apwhistname+="azimuthFWHMprofile_";
    apwhistname+= savename;

    TH1D** TrHist;
    TrHist = new TH1D*[phibins];

    for(int nn=0; nn<phibins; nn++){
      TString trhname;
      trhname+="trhist";
      trhname+=nn;
      TrHist[nn]= new TH1D(trhname,trhname,tbins,-100,500);
    }

    azimuthmeanprofile = new TH1D(apmhistname,apmhistname,phibins,-3.2,3.2);
    azimuthpeakprofile = new TH1D(apphistname,apphistname,phibins,-3.2,3.2);
    azimuthrmsprofile = new TH1D(aprhistname,aprhistname,phibins,-3.2,3.2);
    azimuthFWHMprofile = new TH1D(apwhistname,apwhistname,phibins,-3.2,3.2);
    azimuthvstresid = new TH2D(ahistname,ahistname,(phibins*10),-3.14159654,3.141592654,(tbins/2),-15.,40.);
    azimuthvstheta = new TH2D(athistname,athistname,(phibins*10),-3.141592654,3.141592654,(tbins/2),0.,1.5708);

    for(int mm=0; mm<(_tresid.size()); mm++){
      
      vector<double> thisresid = _tresid.at(mm);
      vector<double> thisazimuth = _azimuth.at(mm);
      vector<double> thisthet = _ttheta.at(mm);
      double hittime = (theHits.at(mm)).at(3);
      double myazimuth = thisazimuth.at(0);
      double mytresid = thisresid.at(0);
      double myweight = thisresid.at(1);

      int whichbin = azimuthmeanprofile->FindBin(myazimuth);

      if( (hittime>lowt) && (hittime<hight) ){
	TrHist[(whichbin-1)]->Fill(mytresid,myweight);
	azimuthvstresid->Fill((thisazimuth.at(0)),(thisresid.at(0)),(thisresid.at(1)));
	azimuthvstheta->Fill((thisazimuth.at(0)),(thisthet.at(0)),(thisresid.at(1)));
      }
    }

    for(int nn=0; nn<phibins; nn++){
      
      double nevents = TrHist[nn]->GetEntries();
      double rmsres = TrHist[nn]->GetRMS();
      double meanres = TrHist[nn]->GetMean();
      double reserror;
      if(nevents>0){
	reserror = (rmsres/sqrt(nevents));
      } else{
	reserror=0;
      }

      // get peak value
      int maxbin = TrHist[nn]->GetMaximumBin();
      double maxvalue = ((TrHist[nn]->GetBinCenter(maxbin)) + (TrHist[nn]->GetBinCenter(maxbin-2)) + (TrHist[nn]->GetBinCenter(maxbin-1)) 
			 + (TrHist[nn]->GetBinCenter(maxbin+1)) + (TrHist[nn]->GetBinCenter(maxbin+2)))/5.;

      double maxcontent = TrHist[nn]->GetBinContent(maxbin);
      double halfmax = maxcontent*0.5;
      int ntrbins = TrHist[nn]->GetNbinsX();

      int lowbin,highbin;
      int crossedonce =0;
      for(int ll=1; ll<(ntrbins+1); ll++){
	double bc = TrHist[nn]->GetBinContent(ll);
	if( (bc>halfmax) && (crossedonce==0) ){
	  lowbin = ll;
	  crossedonce++;
	}
	if( (bc<halfmax) && (crossedonce==1) ){
	  highbin = ll;
	  crossedonce++;
	}
      }

      double fwhm = ((TrHist[nn]->GetBinCenter(highbin)) - (TrHist[nn]->GetBinCenter(lowbin)));


      azimuthrmsprofile->SetBinContent((nn+1),rmsres);
      azimuthrmsprofile->SetBinError((nn+1),(reserror/sqrt(2)));

      azimuthFWHMprofile->SetBinContent((nn+1),fwhm);
      azimuthFWHMprofile->SetBinError((nn+1),(reserror/sqrt(2)));

      azimuthmeanprofile->SetBinContent((nn+1),meanres);
      azimuthmeanprofile->SetBinError((nn+1),reserror);
      azimuthpeakprofile->SetBinContent((nn+1),maxvalue);
      azimuthpeakprofile->SetBinError((nn+1),reserror); 
      
    }
 
    for(int nn=0; nn<phibins; nn++){
      delete TrHist[nn];    
    }
    
}







// ClassImp(genGauss)
