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
  TString filename("../../FullEvent_3k_-50.root");
  TString genfilename("../../generatorcardfile_3k_-50.root");

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

  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  
  cout<<nentries<<endl;

  //  for(int i=0; i<nentries; i++){
  for(int i=0; i<200; i++){

    if(i%10==0) cout<<"event: "<<i<<endl;
    cout<<"**************************************************************"<<endl;
    mTR->LoadEvent(i);
    int nneut = mTR->get_neutroncount();
    int nneut_gen = mTR->get_gennneutrons();
    int npart = mTR->get_npart();
    int ncapcount = mTR->get_ncapturecount();
    
    Int_t* partpid = mTR->get_part_pid();

    Double_t* partxstart = mTR->get_part_xStart();
    Double_t* partystart = mTR->get_part_yStart();
    Double_t* partzstart = mTR->get_part_zStart();
    Int_t* partprocessend = mTR->get_part_processEnd();

    Double_t* partxend = mTR->get_part_xEnd();
    Double_t* partyend = mTR->get_part_yEnd();
    Double_t* partzend = mTR->get_part_zEnd();
    Int_t* partprocessstart = mTR->get_part_processStart();

    vector<double> neut_xstart;    
    vector<double> neut_ystart;    
    vector<double> neut_zstart;    
    vector<double> neut_xend;    
    vector<double> neut_yend;    
    vector<double> neut_zend;    

    vector<int> neut_pstart;
    vector<int> neut_pend;


    for(int j=0; j<npart; j++){

      if(partpid[j]==2112){
	
	neut_xstart.push_back(partxstart[j]);
	neut_ystart.push_back(partystart[j]);
	neut_zstart.push_back(partzstart[j]);

	neut_xend.push_back(partxend[j]);
	neut_yend.push_back(partyend[j]);
	neut_zend.push_back(partzend[j]);

	neut_pstart.push_back(partprocessstart[j]);
	neut_pend.push_back(partprocessend[j]);
      }
    }

    int newnneut=0;
    
    for(int k=0; k<nneut; k++){
      
      if(neut_pend!=11){
	
	int ismatched=0;
	
	for( int mm=0; mm<nneut; mm++){
	  
	  if(neut_xend[k]==neut_xstart[mm]){
	    
	    cout<<"I think we have a match"<<endl;
	    cout<<neut_xend[k]<<" "<<neut_yend[k]<<" "<<neut_zend[k]<<" "<<neut_pend[k]<<endl;
	    cout<<neut_xstart[mm]<<" "<<neut_ystart[mm]<<" "<<neut_zstart[mm]<<" "<<neut_pstart[mm]<<endl;
	    ismatched = 1;
	  }
	}
	if(ismatched==0){cout<<"UNMATCHED!!"<<endl; newnneut++;}
      } else{newnneut++;}
      
    }
    
    cout<<"original nneut: "<<nneut<<" corrected: "<<newnneut<<" capture count: "<<ncapcount<<endl;
    cout<<"**************************************************************"<<endl;

  }

}

