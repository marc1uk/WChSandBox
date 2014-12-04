{

  //************************************************************//
  //								//
  // 	Set the value of "opt" to choose between 		//
  //	TITUS or a ANNIE detector.				//
  //								//
  // 	Set the value of "config" to choose the coverage	//
  //	of the PMTs.						//
  //								//
  // 	Set the value of "LAPPDs" to choose if there is		//
  //	LAPPDs or not.						//
  //								//
  //	And set other values to choose the time resolution,	//
  //	quantum efficiency, size, shape of the PMTs or LAPPDs	//
  //								//
  //************************************************************//


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
  //  TString filename("../../FullEvent.root");
  //  TString genfilename("../../generatorcardfile.root");
    TString filename("../../SmearedEventTree.root");
  //TString genfilename("../../generatorcardfile.root");

  // Load Data
  // =========
  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree
    //  WCLTreeReader *mTR = new WCLTreeReader(0,1);
      WCLTreeReader *mTR = new WCLTreeReader(1,0);



  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)
  //  mTR->LoadData(filename,genfilename);  
  mTR->LoadData(filename);


  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

   vector< vector<double> > hitsarray;

  // Loop over number of entries
  for(int i=0; i<nentries; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

    hitsarray.clear();

    // Load the event variables into TreeReader
    mTR->LoadEvent(i);

    // loop over photons
    int nphot = mTR->get_nphot();
    int npart = mTR->get_npart();

    /// grab the relevant variables from the TreeReader event
    Int_t* phot_hit = mTR->get_phot_hit();
    Int_t* phot_isScat = mTR->get_phot_isScat();
    Int_t* phot_processStart = mTR->get_phot_processStart();
    Int_t* phot_trackid = mTR->get_phot_trackid(); 
    Int_t* phot_parentid = mTR->get_phot_parentid();
    Int_t* phot_capnum = mTR->get_phot_capnum();
    Int_t* part_processEnd = mTR->get_part_processEnd();
    Int_t* part_pid = mTR->get_part_pid();
    Double_t* phot_xStart = mTR->get_phot_xStart();
    Double_t* phot_yStart = mTR->get_phot_yStart();
    Double_t* phot_zStart = mTR->get_phot_zStart();
    Double_t* phot_tStart = mTR->get_phot_tStart();
    Double_t* phot_xEnd = mTR->get_phot_xEnd();
    Double_t* phot_yEnd = mTR->get_phot_yEnd();
    Double_t* phot_zEnd = mTR->get_phot_zEnd();
    Double_t* phot_tEnd = mTR->get_phot_tEnd();

    Double_t* phot_wavelength = mTR->get_phot_wavelength();

 

    //Loop over photons to add those which hit a PMT or LAPPD
    for(int k=0; k<nphot; k++){
      if( (k%10000==0) && (i%10==0) ) cout<<k<<" photons out of "<<nphot<<endl;  
      int PMTid=0;
      int ishit=0;
      int isscat = phot_isScat[k];
      int process = phot_processStart[k];
      int trackid = phot_trackid[k];
      int parentid = phot_parentid[k];
      int capnum = phot_capnum[k];
      int Nrow;
      int Ncol;
      int mcase;
      double xS = phot_xStart[k]; double yS = phot_yStart[k]; double zS = phot_zStart[k]; 
      double tS = phot_tStart[k];
      double xE = phot_xEnd[k]; double yE = phot_yEnd[k]; double zE = phot_zEnd[k];
      double tE = phot_tEnd[k];
      double wl = phot_wavelength[k];

      if(tE<45.0){
      vector<double> hitcoordinates;
      hitcoordinates.push_back(xE);
      hitcoordinates.push_back(yE);      
      hitcoordinates.push_back(zE);
      hitcoordinates.push_back(tE);
      hitcoordinates.push_back(wl);

      hitsarray.push_back(hitcoordinates);
      }

    }
   
    TString dname;
    dname+="Event_";
    dname+=i;
    cout<<"Hits array size "<<hitsarray.size()<<endl;
    EventDisplay* mED = new EventDisplay("covgdisplayB");
    mED->SetGeometry(1,2500,3500,3500,10.);
    mED->SetPMTGeometry(1,(203/1.414),0,1);
    mED->SetHits(hitsarray);



     mED->SavePlotsBarrel("testteset",1,0,1000);

  }



}  
