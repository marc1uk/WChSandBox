{

  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");
  gSystem->Load("../lib/libWCLAnalysis.so");

  #include "TRandom3.h"
  #include <iomanip>
  #include "../include/WCLTreeReader.hh"
  #include "../include/SandBoxPMTcoverage.hh"

  TRandom3 numberp(42484);


  // Path to WCSim ROOT file
  // =======================
  TString filename("testcoverage.root");

  TFile *tf = new TFile(filename,"READ");
  TTree* mt = (TTree*) tf->Get("SmearedEventTree");
  
  TString whattodraw;
  TString conditions;

  whattodraw+="phot_xEnd:phot_yEnd:phot_zEnd:phot_tEnd";
  conditions+="evt==1&&phot_tEnd<30";

  //whattodraw+="phot_xEnd:phot_yEnd:phot_tEnd";
  //conditions+="evt==1&&phot_tEnd<30&&phot_zEnd==1500.";

  mt->Draw(whattodraw,conditions,"colz");

}

