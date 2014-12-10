#include "EventDisplay3D.hh"
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"

#include "TGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveViewer.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TEveStraightLineSet.h"
#include "TEveArrow.h"
#include "TEveText.h"
#include "TStyle.h"


#include <iostream>
#include <cmath>


ClassImp(EventDisplay3D)

EventDisplay3D::EventDisplay3D()
{
  TEveManager::Create();

  // as a default, do not draw gen event
  _drawgen = false;	  

  // create new TStyle
  TStyle *fRootStyle = new  TStyle("My Style", "");

  // set the background color to white
  fRootStyle->SetFillColor(10);
  fRootStyle->SetFrameFillColor(10);
  fRootStyle->SetCanvasColor(10);
  fRootStyle->SetPadColor(10);
  fRootStyle->SetTitleFillColor(0);
  fRootStyle->SetStatColor(10);

  // don't put a colored frame around the plots
  fRootStyle->SetFrameBorderMode(0);
  fRootStyle->SetCanvasBorderMode(0);
  fRootStyle->SetPadBorderMode(0);
  fRootStyle->SetLegendBorderSize(0);

  // use the primary color palette
  fRootStyle->SetPalette(1,0);

  // set the default line color for a histogram to be black
  fRootStyle->SetHistLineColor(kBlack);

  // set the default line color for a fit function to be red
  fRootStyle->SetFuncColor(kRed);

  // make the axis labels black
  fRootStyle->SetLabelColor(kBlack,"xyz");

  // set the default title color to be black
  fRootStyle->SetTitleColor(kBlack);

  // set the margins
  fRootStyle->SetPadBottomMargin(0.18);
  fRootStyle->SetPadTopMargin(0.08);
  fRootStyle->SetPadRightMargin(0.08);
  fRootStyle->SetPadLeftMargin(0.17);

  // set axis label and title text sizes
  fRootStyle->SetLabelFont(42,"xyz");
  fRootStyle->SetLabelSize(0.06,"xyz");
  fRootStyle->SetLabelOffset(0.015,"xyz");
  fRootStyle->SetTitleFont(42,"xyz");
  fRootStyle->SetTitleSize(0.07,"xyz");
  fRootStyle->SetTitleOffset(1.1,"yz");
  fRootStyle->SetTitleOffset(1.0,"x");
  fRootStyle->SetStatFont(42);
  fRootStyle->SetStatFontSize(0.07);
  fRootStyle->SetTitleBorderSize(0);
  fRootStyle->SetStatBorderSize(0);
  fRootStyle->SetTextFont(42);

  // set line widths
  fRootStyle->SetFrameLineWidth(2);
  fRootStyle->SetFuncWidth(2);

  // set the number of divisions to show
  fRootStyle->SetNdivisions(506, "xy");

  // turn off xy grids
  fRootStyle->SetPadGridX(0);
  fRootStyle->SetPadGridY(0);

  // set the tick mark style
  fRootStyle->SetPadTickX(1);
  fRootStyle->SetPadTickY(1);

  // turn off stats
  fRootStyle->SetOptStat(0);

  // marker settings
  fRootStyle->SetMarkerStyle(20);
  fRootStyle->SetMarkerSize(1.2);
  fRootStyle->SetLineWidth(1);

  // done
  fRootStyle->cd();
//  gROOT->ForceStyle();
//  gStyle->ls();


  // set default time ranges
  _tb[0] = 5.0;
  _tb[1] = 10.0;
  _tb[2] = 15.0;
  _tb[3] = 20.0;
  _tb[4] = 25.0;
  _tb[5] = 30.0;
  _tb[6] = 40.0;
  _tb[7] = 50.0;
         

}

EventDisplay3D::~EventDisplay3D()
{

}

void EventDisplay3D::DrawDetector(){

  TGeoManager *geom = new TGeoManager("DetectorGeometry", "detector geometry");

  // materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0.0,0.0,0.0);
  TGeoMaterial *matWater = new TGeoMaterial("Water",18.0,8.0,1.0);

  // media
  TGeoMedium* Vacuum = new TGeoMedium("Vacuum",1,matVacuum);
  TGeoMedium* Water = new TGeoMedium("Water",2,matWater);

  // top volume
  TGeoVolume *top = geom->MakeBox("Detector", Vacuum, 10000.0, 10000.0, 10000.0);
  geom->SetTopVolume(top);

  Double_t fCylRadius = 1250.;
  Double_t fCylLength = 3500.;

  TGeoVolume* myCylinder = geom->MakeTube("Cylinder",Water,0.0,fCylRadius,0.5*fCylLength);
  myCylinder->SetLineColor(kCyan);
  myCylinder->SetTransparency(70);  // percentage transparency [0-100]
  myCylinder->SetVisibility(1);     
  top->AddNode( myCylinder, 0, new TGeoTranslation(0, 0, 0));

  // close geometry
  geom->CloseGeometry();

  std::cout<<"Made is this far!"<<std::endl;
	

  // Create Geometry in Event Display
  // ================================

  TGeoNode* node = gGeoManager->GetTopNode();
  TEveGeoTopNode* eveNode = new TEveGeoTopNode(gGeoManager, node);
  eveNode->SetVisLevel(1);
  gEve->AddGlobalElement(eveNode);



  // Draw Display
  // ============
  gEve->Redraw3D(kTRUE);



}



void EventDisplay3D::PlotEvent(){





  // Containers for Hits
  // ===================
  Int_t markerStyle = 4;   
  Double_t markerSize = 0.25;

  Int_t colourCode1 = kBlue+1;
  Int_t colourCode2 = kCyan+1;
  Int_t colourCode3 = kGreen;
  Int_t colourCode4 = kYellow;
  Int_t colourCode5 = kOrange;
  Int_t colourCode6 = kOrange+7;
  Int_t colourCode7 = kRed;
  Int_t colourCode8 = kRed;
  Int_t colourCode9 = kRed;

  TEvePointSet* fHitList1 = new TEvePointSet(); 
  fHitList1->SetMarkerColor(colourCode1);
  fHitList1->SetMarkerSize(markerSize);
  fHitList1->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList1);

  TEvePointSet* fHitList2 = new TEvePointSet(); 
  fHitList2->SetMarkerColor(colourCode2);
  fHitList2->SetMarkerSize(markerSize);
  fHitList2->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList2);

  TEvePointSet* fHitList3 = new TEvePointSet();  
  fHitList3->SetMarkerColor(colourCode3);
  fHitList3->SetMarkerSize(markerSize);
  fHitList3->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList3);

  TEvePointSet* fHitList4 = new TEvePointSet();  
  fHitList4->SetMarkerColor(colourCode4);
  fHitList4->SetMarkerSize(markerSize);
  fHitList4->SetMarkerStyle(markerStyle);  
  gEve->AddElement(fHitList4);

  TEvePointSet* fHitList5 = new TEvePointSet();  
  fHitList5->SetMarkerColor(colourCode5);
  fHitList5->SetMarkerSize(markerSize);
  fHitList5->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList5);

  TEvePointSet* fHitList6 = new TEvePointSet();  
  fHitList6->SetMarkerColor(colourCode6);
  fHitList6->SetMarkerSize(markerSize);
  fHitList6->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList6);

  TEvePointSet* fHitList7 = new TEvePointSet();  
  fHitList7->SetMarkerColor(colourCode7);
  fHitList7->SetMarkerSize(markerSize);
  fHitList7->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList7);

  TEvePointSet* fHitList8 = new TEvePointSet(); 
  fHitList8->SetMarkerColor(colourCode8);
  fHitList8->SetMarkerSize(markerSize);
  fHitList8->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList8);

  TEvePointSet* fHitList9 = new TEvePointSet();
  fHitList9->SetMarkerColor(colourCode9);
  fHitList9->SetMarkerSize(markerSize);
  fHitList9->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitList9);

  TEvePointSet* fHitListGen = new TEvePointSet();
  fHitListGen->SetMarkerColor(kWhite);
  fHitListGen->SetMarkerSize(markerSize);
  fHitListGen->SetMarkerStyle(markerStyle);
  gEve->AddElement(fHitListGen);


  // Loop over digits
  // ================
  for( Int_t ihit=0; ihit<_hits.size(); ihit++ ){

    std::vector<double> vhit = _hits.at(ihit);
     
    double x = vhit.at(0);
    double y = vhit.at(1);
    double z = vhit.at(2);
    double t = vhit.at(3);
    double wl = vhit.at(4);
    double Q = vhit.at(5);

    // in future
    // if( Q<GetPulseHeightCut() ) continue;

    Int_t listNumber = 0;

    // color code by charge
    /*
    if( Q<0.8 )             listNumber = 1;
    if( Q>=0.8 && Q<1.5 )   listNumber = 2;
    if( Q>=1.5 && Q<2.5 )   listNumber = 3;
    if( Q>=2.5 && Q<5.0 )   listNumber = 4;
    if( Q>=5.0 && Q<10.0 )  listNumber = 5;
    if( Q>=10.0 && Q<15.0 ) listNumber = 6;
    if( Q>=15.0 && Q<20.0 ) listNumber = 7;
    if( Q>=20.0 && Q<30.0 ) listNumber = 8;
    if( Q>=30.0 )           listNumber = 9;
    */

    if( t<_tb[0] )               listNumber = 1;
    if( t>=_tb[0] && t<_tb[1] )  listNumber = 2;
    if( t>=_tb[1] && t<_tb[2] )  listNumber = 3;
    if( t>=_tb[2] && t<_tb[3] )  listNumber = 4;
    if( t>=_tb[3] && t<_tb[4] )  listNumber = 5;
    if( t>=_tb[4] && t<_tb[5] )  listNumber = 6;
    if( t>=_tb[5] && t<_tb[6] )  listNumber = 7;
    if( t>=_tb[6] && t<_tb[7] )  listNumber = 8;
    if( t>=_tb[7] )              listNumber = 9;



    switch( listNumber ){
      case 1: fHitList1->SetNextPoint(x,y,z); break;
      case 2: fHitList2->SetNextPoint(x,y,z); break;
      case 3: fHitList3->SetNextPoint(x,y,z); break;
      case 4: fHitList4->SetNextPoint(x,y,z); break;
      case 5: fHitList5->SetNextPoint(x,y,z); break;
      case 6: fHitList6->SetNextPoint(x,y,z); break;
      case 7: fHitList7->SetNextPoint(x,y,z); break;
      case 8: fHitList8->SetNextPoint(x,y,z); break;
      case 9: fHitList9->SetNextPoint(x,y,z); break;
      default: break;
    }
  }

  // Loop over true photon starting positions
  // ================


  if(_drawgen){
  	for( Int_t ipos=0; ipos<_genpositions.size(); ipos++ ){

   		 std::vector<double> vpos = _genpositions.at(ipos);
     
  		  double gx = vpos.at(0);
  		  double gy = vpos.at(1);
   	  	  double gz = vpos.at(2);
  		  double gt = vpos.at(3);
  		  //double gwl = vpos.at(4);

                 fHitListGen->SetNextPoint(gx,gy,gz); 
        }
   }


  // Re-draw Event Display
  // =====================
  gEve->Redraw3D();

  return;
}

