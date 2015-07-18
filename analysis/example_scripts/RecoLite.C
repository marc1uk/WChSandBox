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
  std::vector<double> fDigitX;
  std::vector<double> fDigitY;
  std::vector<double> fDigitZ;
  std::vector<double> fDigitT;
  std::vector<double> fDigitQ;
  std::vector<double> fDigitPE;
  std::vector<double> fDigitW;
  std::vector<double> fDigitV;
  std::vector<double> fDelta; // time residual

  std::vector<double> vSeedVtxX;
  std::vector<double> vSeedVtxY;
  std::vector<double> vSeedVtxZ;
  std::vector<double> vSeedVtxTime;
  std::vector<int> vSeedDigitList;

  int NMAX_PHOT=100000; //
  int NPHI=12;
  int NTHETA=12;

  //  int EVT_NUM=100; //controls maximum number of events to be processed
  double R_SPHERE=650; //sphere diameter [cm]
  double N_REF=1.53; //average index of refraction
  double C_VAC=29.9792458; //speed of light in vacuum [cm/ns]
  int NSeedsTarget=400; //number of quadruplets
  double TSIGMA=0.5; //total time spread (including detector TTS chromatic dispersions)

  // Path to WCSim ROOT file
  // =======================
  TString filename("../../FullEvent_e.root");
  TString genfilename("../../generatorcardfile_e.root");

  // Load Data
  // =========
  
  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree

  WCLTreeReader *mTR = new WCLTreeReader(0,1);

  // Initialize WChRecoLite
  WChRecoLite* mRec = WChRecoLite::Instance();

  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)

  mTR->LoadData(filename,genfilename);

  cout<<"initializing treewriter"<<endl;
  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  //define reconstrution output here
  TString fOutputName;
  fOutputName+="testout.root";

  TFile f_out(fOutputName,"recreate");
  double recoVtxX;
  double recoVtxY;
  double recoVtxZ;
  double recoVtxTime;

  Double_t trueVtxX; 
  Double_t trueVtxY; 
  Double_t trueVtxZ;
  Double_t trueVtxTime;

  double DrecoVtxX;
  double DrecoVtxY;
  double DrecoVtxZ;
  double DrecoVtxTime;

    
  TTree* reco_out_ntuple = new TTree("ntuple","ntuple");
  reco_out_ntuple->Branch("recoVtxX",&recoVtxX,"recoVtxX/D");
  reco_out_ntuple->Branch("recoVtxY",&recoVtxY,"recoVtxY/D");
  reco_out_ntuple->Branch("recoVtxZ",&recoVtxZ,"recoVtxZ/D");
  reco_out_ntuple->Branch("recoVtxTime",&recoVtxTime,"recoVtxTime/D");

  reco_out_ntuple->Branch("trueVtxX",&trueVtxX,"trueVtxX/D");
  reco_out_ntuple->Branch("trueVtxY",&trueVtxY,"trueVtxY/D");
  reco_out_ntuple->Branch("trueVtxZ",&trueVtxZ,"trueVtxZ/D");
  reco_out_ntuple->Branch("trueVtxTime",&trueVtxTime,"trueVtxTime/D");

  reco_out_ntuple->Branch("DrecoVtxX",&DrecoVtxX,"DrecoVtxX/D");
  reco_out_ntuple->Branch("DrecoVtxY",&DrecoVtxY,"DrecoVtxY/D");
  reco_out_ntuple->Branch("DrecoVtxZ",&DrecoVtxZ,"DrecoVtxZ/D");
  reco_out_ntuple->Branch("DrecoVtxTime",&DrecoVtxTime,"DrecoVtxTime/D")

  //====================
  //results of the first itteration of reco is defined here
  double X0, Y0, Z0;
  //  int NNN=1;
  TTree* recoTree;
  /*
  if(fRecoIt==2)
  {
    TFile* recoFile = new TFile(fFirstRecoName);
    recoTree = (TTree*)recoFile->Get("ntuple");
    recoTree->SetBranchAddress("recoVtxX",&X0);
    recoTree->SetBranchAddress("recoVtxY",&Y0);
    recoTree->SetBranchAddress("recoVtxZ",&Z0);
  }
  */
  //=====================


  for(int i=0; i<nentries; i++){
  //  for(int i=0; i<20; i++){

    int currentEvent=i;

    fThisDigit=0;
    fDigitX.clear();
    fDigitY.clear();
    fDigitZ.clear();
    fDigitT.clear();
    fDigitW.clear();
    fDigitV.clear();
    fDigitQ.clear();
    fDigitPE.clear();
    vSeedDigitList.clear();

    if(i%10==0) cout<<"event: "<<i<<endl;

    mTR->LoadEvent(i);

    int N_phot_v = mTR->get_nphot();
    Double_t* x_hit_v = mTR->get_phot_xEnd();
    Double_t* y_hit_v = mTR->get_phot_yEnd();
    Double_t* z_hit_v = mTR->get_phot_zEnd();
    Double_t* true_time_corrected_v = mTR->get_phot_tEnd();
    Double_t* true_time_v = mTR->get_phot_tEnd();
    Double_t* PE_time_v = mTR->get_phot_tEnd();
    Double_t* photon_wavelength_v = mTR->get_phot_wavelength();
    Double_t* process_v = mTR->get_phot_processStart();

    trueVtxX = mTR->get_genvtxx();
    trueVtxY = mTR->get_genvtxy();
    trueVtxZ = mTR->get_genvtxz();
    trueVtxTime = 0;

    cout<<"TV!! :"<<trueVtxX<<" "<<trueVtxY<<" "<<trueVtxZ<<endl;

    for(int iphot=0;iphot!=N_phot_v;iphot++)
    {
        double distL = TMath::Sqrt(TMath::Power((x_hit_v[iphot] - X0*10), 2) + (y_hit_v[iphot]-Y0*10)*(y_hit_v[iphot]-Y0*10) + (z_hit_v[iphot]-Z0*10)*(z_hit_v[iphot]-Z0*10));
        double light_vel = 300./1.53;
        double TPredicted = distL/light_vel;
        if((PE_time_v[iphot] - TPredicted) > 3.0) continue;
	
        fDigitX.push_back(x_hit_v[iphot]/10.); //!Sphere1
        fDigitY.push_back(y_hit_v[iphot]/10.);//!Sphere1
        fDigitZ.push_back(z_hit_v[iphot]/10.);//!Sphere1
	fDigitT.push_back(PE_time_v[iphot]);// - min_PE_time; //!Sphere1
	fDigitQ.push_back(1);
        fDigitW.push_back(photon_wavelength_v[iphot]);
//	int lambda = (int)
//	fDigitV.push_back(C_VAC/INDEX[(int)photon_wavelength_v[iphot]]);
	fDigitV.push_back(C_VAC/1.33);
	vSeedDigitList.push_back(fThisDigit);
        fThisDigit++;
    }

    std::cout<<"Photon filtering for event #"<<i<<" has just finished. fDigits are ready."<<std::endl;
    
    //input digits
    mRec->SetDigits(fDigitX, fDigitY, fDigitZ, fDigitT, fDigitQ, fDigitPE, fDigitW, fDigitV, fDelta, vSeedDigitList);

    
    //calculate seed vertices
    mRec->CalcVertexSeeds();
    
    vSeedVtxX = mRec->GetSeedVtx(0);
    vSeedVtxY = mRec->GetSeedVtx(1);
    vSeedVtxZ = mRec->GetSeedVtx(2);
    vSeedVtxTime = mRec->GetSeedVtx(3);
    
    //now select the best vertex
    int best_seed = mRec->SelectBestSeed(currentEvent);
    std::cout<<"best_seed = "<<best_seed<<"  vSeedVtxX.size() = "<<vSeedVtxX.size()<<std::endl;
    //and save reco data
    
    //   cout<<"UUUUU "<<vSeedVtxTime.size()<<endl;

    recoVtxX = vSeedVtxX[best_seed];
    recoVtxY = vSeedVtxY[best_seed];
    recoVtxZ = vSeedVtxZ[best_seed];
    recoVtxTime = vSeedVtxTime[best_seed];

    DrecoVtxX = vSeedVtxX[best_seed]-trueVtxX;
    DrecoVtxY = vSeedVtxY[best_seed]-trueVtxY;
    DrecoVtxZ = vSeedVtxZ[best_seed]-trueVtxZ;
    DrecoVtxTime = vSeedVtxTime[best_seed]-trueVtxTime;
    
    reco_out_ntuple->Fill();
    
  } //end i-loop over Hits_Tree entries
  
  f_out.cd();
  reco_out_ntuple->Write();
  return 0;
 
}

