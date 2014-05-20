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

  cout<<"initializing treewriter"<<endl;

  // WCLTREEWRITER
  // *********************************************************
  // specify the output filename and the2nd parameter is the  mode
  // mode=1 is the only relevant mode for most people
  WCLTreeWriter *mTW = new WCLTreeWriter("testout_Titus_200.root",1);

  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  
  cout<<nentries<<endl;

  for(int i=0; i<nentries; i++){

    if(i%100==0) cout<<"event: "<<i<<endl;

    mTR->LoadEvent(i);
    int nphot = mTR->get_nphot();
    Double_t* photxstart = mTR->get_phot_xStart();

    mTW->InitializeEvent();

    // AddWholeBranches(treereader,copyphot,copypart,copycapt,copygen,copymrd)
    //                 (WCLTreeReader*, int, int, int, int, int)
    //**************************
    // the Root trees store data for several categories of arrays:
    // photons, non-photon particles, neutron captures, gen level, and MRD
    // AddWholeBranches let's you pass the TreeReader in to the Tree Writer
    // and copy the entire contents of one or more sets of arrays, without
    // necessarily copying over another. In this example, we copy over
    // all of the non-photon particle information but no other data.
    // We leave the photon information blank, because we are going to
    // manually input the photons.

    mTW->AddWholeBranches(mTR,1,1,1,1,1);

    mTW->FillEvent();

  }

  mTW->WriteTreeToFile();

}

