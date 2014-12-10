//////////////////////////////////////////////////////////////////////
///  EventDisplay2D.hh
///  Sept 4, 2010 Matt Wetstein
///



#ifndef EventDisplay2D_hh
#define EventDisplay2D_hh

#include <TROOT.h>
#include <vector>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TString.h>
#include <TApplication.h>
#include <vector>
#include <TRandom.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TVector3.h>
#include <string.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>


#include <iostream>
#include <fstream>

using namespace std;

class EventDisplay2D
{

 public:
  
  //Do Something

  EventDisplay2D(TString iname);
  ~EventDisplay2D();
  void SetHits(vector< vector<double> > ihits);
  void SetAzimuth(vector< vector<double> > tresid, vector< vector<double> > azimuth, vector< vector<double> > tracktheta);
  void SetGeometry(int type, double d1, double d2, double d3, double clearance);
  void SetPMTGeometry(int pmttype,double diameter, double coverage, bool spatialresolving);
  void SetTrueTrackParams(vector<double> tparams);
  void SetPlots(TH1D* trD,TH1D* tpD, TH1D* ttD);
  void PlotEventBarrel(string savename, bool savecanvas, double lowt, double hight, bool dowrite);
  void PlotEventCube(string savename, bool savecanvas, double lowt, double hight, bool dowrite);
  void PlotEventThetaPhi(string savename, bool savecanvas, double lowt, double hight, bool dowrite);
  void PlotEventTime(string savename, bool savecanvas, double lowt, double hight, bool dowrite);


//  void SavePlotsBarrel(string savename, bool savecanvas, double lowt, double hight);
//  void SavePlotsCube(string savename, bool savecanvas, double lowt, double hight);
  void doAzimuthalPlots(int tbins, int phibins, string savename, double lowt, double hight);
  void InitializeStyle();
  TCanvas* MakeCanvas(int ctype, TString namestring, TH2D* &iwcDisplay);
  TCanvas* MakeCanvasCylinder(TString namestring, TH2D* &iwcDisplay);
  TCanvas* MakeCanvasCube(TString namestring, TH2D* &iwcDisplay);
  TCanvas* MakeCanvasThetaPhi(TString namestring, TH2D* &iwcDisplay);
  TCanvas* MakeCanvasTime(TString namestring, TH2D* &iwcDisplay);



  //  ClassDef(EventDisplay, 0);
  
 private:
  
  // name of EventDisplay instance

  TString _iname;

  // scale of the plots

  double fScale;

  //booleans to keep track of what's been done
  bool _isweighted;
  bool _smeared;
  bool _initiatedPMTGeo;
  bool _initiatedDetectorGeo;
  bool _setThePlots;
  bool _openedoutfile;
  bool _settruetrackparams;
  bool _setazimuth;

  //the main data structures
  vector< vector<double> > theHits;
  vector< vector<double> > _tresid;
  vector< vector<double> > _azimuth;
  vector< vector<double> > _ttheta;

  //track parameters
  vector<double> truetrack;
  vector<double> fittrack;

  //parameters to model individual PMTs
  int _pmt_type;
  double _spatialcoverage;
  double _PMTdiameter;
  bool _spatialresolving;

  //detector geometry parameters
  int _detectortype;
  double _d1; // diameter if cylindrical
  double _d2;
  double _d3; // only used for mailbox geometry
  double _clearance;

  //style class
  TStyle *rootStyle; 

  //data-hypothesis comparison histograms
  TH1D* pointresidualDist;
  TH1D* trackresidualDist;
  TH1D* trackposDist;
  TH1D* trackthetaDist;

  TH2D* azimuthvstresid;
  TH2D* azimuthvstheta;
  TH1D* azimuthmeanprofile;
  TH1D* azimuthrmsprofile;
  TH1D* azimuthpeakprofile;
  TH1D* azimuthFWHMprofile;


  //  TFile* plotfile;
};

#endif
