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
  WCLTreeWriter *mTW = new WCLTreeWriter("testout.root",1);

  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  
  cout<<nentries<<endl;

  for(int i=0; i<nentries; i++){

    if(i%10==0) cout<<"event: "<<i<<endl;

    mTR->LoadEvent(i);
    int nphot = mTR->get_nphot();
    Double_t* photxstart = mTR->get_phot_xStart();
    if(nphot>5){
      cout<<"x-start of the fifth photon: "<<photxstart[4]<<endl;
    }

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
       mTW->AddWholeBranches(mTR,0,1,0,0,0);
    
    double xs;
      // manually filling instants of WCSimTrueLight and adding them to each event in TreeWriter Using the AddPhoton class
    if(i==0) {
      xs=5.2;
      WCSimTrueLight* mTL1 = new WCSimTrueLight(0,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL1);
      WCSimTrueLight* mTL2 = new WCSimTrueLight(xs+0.2,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL2);
      WCSimTrueLight* mTL3 = new WCSimTrueLight(xs-0.2,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL3);
      WCSimTrueLight* mTL4 = new WCSimTrueLight(xs+3.4,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL4);
    }

    if(i==1){
      xs=3.7;
      WCSimTrueLight* mTL5 = new WCSimTrueLight(0,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL5);
      WCSimTrueLight* mTL6 = new WCSimTrueLight(xs+0.2,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL6);
      WCSimTrueLight* mTL7 = new WCSimTrueLight(xs-0.2,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL7);
      WCSimTrueLight* mTL8 = new WCSimTrueLight(xs+3.4,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL8);
      WCSimTrueLight* mTL9 = new WCSimTrueLight(xs-22.4,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL9);
    }
    if(i==2){
      xs=-1.7;
      WCSimTrueLight* mTL10 = new WCSimTrueLight(0,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL10);
      WCSimTrueLight* mTL11 = new WCSimTrueLight(xs+3.4,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL11);
      WCSimTrueLight* mTL12 = new WCSimTrueLight(xs-22.4,0,0,0,1,1,1,1,432,4*i+3,0,i,i+2,1,-666,12);
      mTW->AddPhoton(mTL12);
    }

    mTW->FillEvent();

  }

  mTW->WriteTreeToFile();

}

