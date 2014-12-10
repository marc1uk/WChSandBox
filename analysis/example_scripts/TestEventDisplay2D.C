{
  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");
  gSystem->Load("lib/EventDisplay.so");

  #include "TRandom3.h"
  #include <iomanip>
  #include "include/EventDisplay3D.hh"

  TRandom3 numberp(42484);


  // Path to WCSim ROOT file
  // =======================
  TString filename("FullEvent.root");
  TString treename("EventTree");
//  TString filename("SmearedEventTree.root");
//  TString treename("SmearedEventTree");

  int ievent=66;

  // Load Data
  // =========

  vector<vector<double> > theHits;
  vector<vector<double> > theGens;


  int knphotmax=10000000;

  int *phot_hit = new double[knphotmax];

  double *phot_xStart = new double[knphotmax];
  double *phot_yStart = new double[knphotmax];
  double *phot_zStart = new double[knphotmax];
  double *phot_tStart = new double[knphotmax];


  double *phot_xEnd = new double[knphotmax];
  double *phot_yEnd = new double[knphotmax];
  double *phot_zEnd = new double[knphotmax];
  double *phot_tEnd = new double[knphotmax];
  double *phot_wavelength = new double[knphotmax];
  int nphot;

  TFile* mtf = new TFile(filename);
  TTree* fChain = (TTree*) mtf->Get(treename);
  cout<<fChain->GetEntries()<<endl;

//  fChain->LoadTree();
  fChain->SetBranchAddress("nphot",&nphot);
  fChain->SetBranchAddress("phot_hit",phot_hit);

  fChain->SetBranchAddress("phot_xEnd",phot_xEnd);
  fChain->SetBranchAddress("phot_yEnd",phot_yEnd);
  fChain->SetBranchAddress("phot_zEnd",phot_zEnd);
  fChain->SetBranchAddress("phot_tEnd",phot_tEnd);
  fChain->SetBranchAddress("phot_wavelength",phot_wavelength);

  fChain->SetBranchAddress("phot_xStart",phot_xStart);
  fChain->SetBranchAddress("phot_yStart",phot_yStart);
  fChain->SetBranchAddress("phot_zStart",phot_zStart);
  fChain->SetBranchAddress("phot_tStart",phot_tStart);


  fChain->GetEntry(ievent); 
  cout<<ievent<<" "<<nphot<<endl;


  TH1D* timehist = new TH1D("timehist","timehist",100,0.,50.);

  for(int j=0; j<nphot; j++){

    timehist->Fill(phot_tEnd[j]);
    vector<double> ahit;
    ahit.push_back(phot_xEnd[j]);
    ahit.push_back(phot_yEnd[j]);
    ahit.push_back(phot_zEnd[j]);
    ahit.push_back(phot_tEnd[j]);
    ahit.push_back(phot_wavelength[j]);



    if(phot_hit[j]>0){
      theHits.push_back(ahit);
    }
  }

    TString dname;
    dname+="Event_";
    dname+=ievent;
    cout<<"Hits array size "<<theHits.size()<<endl;
    EventDisplay2D* mED = new EventDisplay2D("covgdisplayB");
    mED->SetGeometry(1,2500,3500,3500,10.);
    mED->SetPMTGeometry(1,(203/1.414),0,1);
    mED->SetHits(theHits);

    mED->PlotEventBarrel("testteset",1,0,1000, false);


  TCanvas* timedist = new TCanvas();
  timehist->Draw();

}

