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
  TString filename("../../FullEvent.root.root");
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

  cout<<"initializing treewriter"<<endl;


  // WCLTREEWRITER
  // *********************************************************
  // specify the output filename and the2nd parameter is the  mode
  // mode=1 is the only relevant mode for most people
  WCLTreeWriter *mTW = new WCLTreeWriter("testout.root",1);

  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  
  cout<<nentries<<endl;

  for(int i=0; i<nentries; i++){

    if(i%10==0) cout<<"event: "<<i<<endl;

    mTR->LoadEvent(i);
    int nphot = mTR->get_nphot();
    Double_t* photxstart = mTR->get_phot_xStart();
    if((nphot>5)&&(i%10==0)){
      cout<<"x-start of the fifth photon: "<<photxstart[4]<<endl;
    }

    mTW->InputEvent(mTR);
  }

  mTW->WriteTreeToFile();

}

