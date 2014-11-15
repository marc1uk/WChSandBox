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

  TRandom3 numberp(42484);


  // Path to WCSim ROOT file
  // =======================
  TString filename("../../FullEvent.root");
  TString genfilename("../../generatorcardfile.root");

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

  // Setup Sample Histograms

  TH1D* Lmu = new TH1D("lm","lm",200,-10.,600.);

  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  // Loop over number of entries
  //  for(int i=0; i<nentries; i++){
  for(int i=0; i<40; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

    // Load the event variables into TreeReader
    mTR->LoadEvent(i);

    // loop over photons
    int nphot = mTR->get_nphot();
    
    for(int k=0; k<nphot; k++){
      
     //       cout<<i<<" "<<k<<endl;
      
      /// grab the relevant variables from the TreeReader event
      Int_t* phot_hit = mTR->get_phot_hit();
      Int_t* phot_isScat = mTR->get_phot_isScat();
      Double_t* phot_xStart = mTR->get_phot_xStart();
      Double_t* phot_yStart = mTR->get_phot_yStart();
      Double_t* phot_zStart = mTR->get_phot_zStart();
      Double_t* phot_tEnd = mTR->get_phot_tEnd();
      Double_t mvtxx = mTR->get_genvtxx();
      Double_t mvtxy = mTR->get_genvtxy();
      Double_t mvtxz = mTR->get_genvtxz();
      
      
      // Determine distance between true vertex and the starting positions
      // of each photon...fill to a hist
      
      if(phot_hit[k]==1 && phot_tEnd[k]<30. && phot_isScat[k]==0){
	double dxs = (phot_xStart[k] - 10*mvtxx)*(phot_xStart[k] - 10*mvtxx);
	double dys = (phot_yStart[k] - 10*mvtxy)*(phot_yStart[k] - 10*mvtxy);
	double dzs = (phot_zStart[k] - 10*mvtxz)*(phot_zStart[k] - 10*mvtxz);
	double vtxdistance = sqrt(dxs+dys+dzs);
	
	//if(vtxdistance<hibH) cout<<hibH<<" "<<vtxdistance<<endl;
	Lmu->Fill(vtxdistance);
      }
    }   
  }

  Lmu->Draw();
}  

