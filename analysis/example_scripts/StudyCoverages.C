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
  TString filename("../../FullEvent_TITUS500.root");
  TString genfilename("../../generatorcardfile_TITUS500.root");

  // Load Data
  // =========
  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree
  WCLTreeReader *mTR = new WCLTreeReader(0,1);



  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)
  mTR->LoadData(filename,genfilename);


  // WCLTREEWRITER
  // *********************************************************
  // specify the output filename and the2nd parameter is the  mode
  // mode=1 is the only relevant mode for most people
  WCLTreeWriter *mTW = new WCLTreeWriter("testcoverage.root",1);

  // SANDBOXPMTCOVERAGE
  //**********************************************************
  // specify NbyN = the number of LAPPDs in a single row (or column) 
  // in a square grid. 13 is largest you can get before they no longer all fit
  // right now this only works for even arrays. I just need to fix the 
  // formular for odd arrays.
  SandBoxPMTcoverage* sbPMT = new SandBoxPMTcoverage(10);
  // set the dimensions of the detector -x,+x,-y,+y,-z,+z
  sbPMT->SetBoxDimensions(-1500.,1500.,-1500.,1500.,-1500.,1500.);
  //configure each wall
  //                front-wall, PMTs, square-shaped, 203mm, 4x4 grid 
  sbPMT->SetWallConfiguration(0,1,1,203.2,4);
   //                back-wall, totally active, ..irrelevant parameters 
  sbPMT->SetWallConfiguration(1,2,0,0,0);

  //                top-wall, PMTs, circular, 203mm, 6x6 grid 
  sbPMT->SetWallConfiguration(2,1,0,203.2,6)
  //                bottom-wall, PMTs, circular, 203mm, 6x6 grid 
  sbPMT->SetWallConfiguration(3,1,0,203.2,6);

  //                left-wall, inactive, ..irrelevant parameters 
  sbPMT->SetWallConfiguration(4,0,0,0,0);
  //                right-wall, inactive,..irrelevant parameters 
  sbPMT->SetWallConfiguration(5,0,0,0,0);
  
  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  // Loop over number of entries
  for(int i=0; i<nentries; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

    // Load the event variables into TreeReader
    mTR->LoadEvent(i);

    // loop over photons
    int nphot = mTR->get_nphot();

    /// grab the relevant variables from the TreeReader event
    Int_t* phot_hit = mTR->get_phot_hit();
    Int_t* phot_isScat = mTR->get_phot_isScat();
    Int_t* phot_processStart = mTR->get_phot_processStart();
    Int_t* phot_trackid = mTR->get_phot_trackid(); 
    Int_t* phot_parentid = mTR->get_phot_parentid();
    Int_t* phot_capnum = mTR->get_phot_capnum();
    Double_t* phot_xStart = mTR->get_phot_xStart();
    Double_t* phot_yStart = mTR->get_phot_yStart();
    Double_t* phot_zStart = mTR->get_phot_zStart();
    Double_t* phot_tStart = mTR->get_phot_tStart();
    Double_t* phot_xEnd = mTR->get_phot_xEnd();
    Double_t* phot_yEnd = mTR->get_phot_yEnd();
    Double_t* phot_zEnd = mTR->get_phot_zEnd();
    Double_t* phot_tEnd = mTR->get_phot_tEnd();
    Double_t* phot_wavelength = mTR->get_phot_wavelength();

    
    //Intialize a new event
    mTW->InitializeEvent();
    //Add all of the information from the event, except for the photon hits
    mTW->AddWholeBranches(mTR,0,1,1,1,1);
    
    //Loop over photons to add only those which it an LAPPD
    for(int k=0; k<nphot; k++){
     
      if( (k%10000==0) && (i%10==0) ) cout<<k<<" photons out of "<<nphot<<endl;
 
      int PMTid=0;
      int ishit=0;
      int isscat = phot_isScat[k];
      int process = phot_processStart[k];
      int trackid = phot_trackid[k];
      int parentid = phot_parentid[k];
      int capnum = phot_capnum[k];
      double xS = phot_xStart[k]; double yS = phot_yStart[k]; double zS = phot_zStart[k]; 
      double tS = phot_tStart[k];
      double xE = phot_xEnd[k]; double yE = phot_yEnd[k]; double zE = phot_zEnd[k];
      double tE = phot_tEnd[k];
      double wl = phot_wavelength[k];

      //if the photon lands on an LAPPD, according to the PMTcoverage class
      if( sbPMT->isActiveHit(xE, yE, zE, PMTid) ){

	//	cout<<"is a hit"<<endl;

	ishit=2;
	WCSimTrueLight* mTL = new WCSimTrueLight(xS,yS,zS,tS,xE,yE,zE,tE,wl,
						 isscat,process,parentid,trackid,ishit,capnum,PMTid);
	mTW->AddPhoton(mTL);
      }
    }

    // Fill the event
    mTW->FillEvent();
    
  }
  mTW->WriteTreeToFile();
}  

