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
  TString filename("/vols/t2k03/hyperk/output/sim/neut/FullEvent_nu_0.root");
  TString genfilename("/vols/t2k03/hyperk/output/sim/neut/generatorcardfile_nu_0.root");

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
  WCLTreeWriter *mTW = new WCLTreeWriter("/vols/t2k03/hyperk/output/tg1814/testcoverage.root",1);

  // set the option: 0=TITUS, 1=ANNIE
  int opt = 0;

  // set the configuration: 
  // 0 to set the values by yourself
  // 1=40% coverage by 20 inch PMTs
  // 2=40% coverage by 12 inch PMTs
  // 3=40% coverage by  8 inch PMTs
  // 4=30% coverage by 20 inch PMTs
  // 5=30% coverage by 12 inch PMTs
  // 6=30% coverage by  8 inch PMTs
  // 7=20% coverage by 20 inch PMTs
  // 8=20% coverage by 12 inch PMTs
  // 9=20% coverage by  8 inch PMTs
  int config = 1;

  // set the presence of LAPPDs
  bool LAPPDs = true;

  // set the dimensions of the detector.
  double Rdet = 5500.;
  double LdetCylindar = 22000.; 
  double LdetCube = 3000.;
  double Ldet;

  // set the time resolution
  double tResPMT = 2.;
  double tResLAPPD = 0.1;
  // set the quantum efficiency
  double QE_PMT = 22.;
  double QE_LAPPD = 22.;
  // set the size of the detectors
  double sizePMT;
  double sizeLAPPD = 203.2;
  double Sdet;
  // set the PMT coverage
  double cover = 40.;
  // set the shape of the detectors: 0=circular, 1=square
  int shapePMT = 0;
  int shapeLAPPD = 1;

  // Number of PMTs in each wall for the option=0.
  int NFront;
  int NBack;
  int NrowCurve;
  int NcolCurve;
  // set the number of PMTs in each wall for the option=1.
  int NbyN = 5;

  if(config==1) {
    sizePMT = 508.;
    cover = 40.;
  }
  if(config==2) {
    sizePMT = 304.8;
    cover = 40.;
  }
  if(config==3) {
    sizePMT = 203.2;
    cover = 40.;
  }
  if(config==4) {
    sizePMT = 508.;
    cover = 30.;
  }
  if(config==5) {
    sizePMT = 304.8;
    cover = 30.;
  }
  if(config==6) {
    sizePMT = 203.2;
    cover = 30.;
  }
  if(config==7) {
    sizePMT = 508.;
    cover = 20.;
  }
  if(config==8) {
    sizePMT = 304.8;
    cover = 20.;
  }
  if(config==9) {
    sizePMT = 203.2;
    cover = 20.;
  }

  if(shapePMT==0) Sdet=3.14*sizePMT*sizePMT/4;
  if(shapePMT==1) Sdet=sizePMT*sizePMT;
  double Ntot = (Rdet*Rdet*3.14+LdetCylindar*2*Rdet*3.14)/Sdet*cover/100;
  if(opt==0) {
    NFront = sqrt(0.8*Ntot/3.14)-fmod(sqrt(0.8*Ntot/3.14),1);
    NBack = sqrt(0.8*Ntot/3.14)-fmod(sqrt(0.8*Ntot/3.14),1);
    NcolCurve = sqrt(0.4*3.14*Ntot)-fmod(sqrt(0.4*3.14*Ntot),1);
    NrowCurve = (0.8*Ntot/NcolCurve)-fmod(0.8*Ntot/NcolCurve,1);
  }
  if(opt==1) {
    NbyN = sqrt(Ntot/6)-fmod(sqrt(Ntot/6),1);
  }


  // SANDBOXPMTCOVERAGE
  //**********************************************************
  // specify NbyN = the number of PMTs in a single row (or column) 
  // in a square grid.
  SandBoxPMTcoverage* sbPMT = new SandBoxPMTcoverage();

  // set the dimensions of the detector L/2, R
  if( opt==0 ){
  Ldet=LdetCylindar;
  }
  else {
  Ldet=LdetCube;
  }
  sbPMT->SetBoxDimensions(opt,Ldet/2.,Rdet);

  //configure each wall for cylindar option
  if(opt==0){
  //                 front-wall, detectors
  sbPMT->SetWallConfiguration(5,1,shapePMT,sizePMT,NFront,NFront,QE_PMT,shapeLAPPD,sizeLAPPD,NFront,NFront,QE_LAPPD);
   //                back-wall, detectors
  sbPMT->SetWallConfiguration(6,1,shapePMT,sizePMT,NBack,NBack,QE_PMT,shapeLAPPD,sizeLAPPD,NBack,NBack,QE_LAPPD);
  //                 curve-wall, detectors 
  sbPMT->SetWallConfiguration(0,1,shapePMT,sizePMT,NrowCurve,NcolCurve,QE_PMT,shapeLAPPD,sizeLAPPD,NrowCurve,NcolCurve,QE_LAPPD);
  }

  //configure each wall for cube option
  if(opt==1){
  // walls are all the same.
    for(int ii=1;ii<7;ii++){
      sbPMT->SetWallConfiguration(ii,1,shapePMT,sizePMT,NbyN,NbyN,QE_PMT,shapeLAPPD,sizeLAPPD,NbyN,NbyN,QE_LAPPD);
    }
  }

  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  double zEndBarrel[5]={-5500.,-5500.,0.,5500.,11000.};
  int NpartLeaving = 0;
  int NpartThruMRD[5] = {0,0,0,0,0};
  int ModeConf[4][4]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  TH1D *CCOth1 = new TH1D("CCOth1","Energy spectrum of muons that reach MRD for configuration 1",20,0,3000);
  TH1D *CCMEC1 = new TH1D("CCMEC1","CCMEC",20,0,3000);
  TH1D *CCQE1 = new TH1D("CCQE1","CCQE",20,0,3000);
  TH1D *NC1 = new TH1D("NC1","NC",20,0,3000);
  TH1D *CCOth2 = new TH1D("CCOth2","Energy spectrum of muons that reach MRD for configuration 2",20,0,3000);
  TH1D *CCMEC2 = new TH1D("CCMEC2","CCMEC",20,0,3000);
  TH1D *CCQE2 = new TH1D("CCQE2","CCQE",20,0,3000);
  TH1D *NC2 = new TH1D("NC2","NC",20,0,3000);
  TH1D *CCOth3 = new TH1D("CCOth3","Energy spectrum of muons that reach MRD for configuration 3",20,0,3000);
  TH1D *CCMEC3 = new TH1D("CCMEC3","CCMEC",20,0,3000);
  TH1D *CCQE3 = new TH1D("CCQE3","CCQE",20,0,3000);
  TH1D *NC3 = new TH1D("NC3","NC",20,0,3000);
  TH1D *CCOth4 = new TH1D("CCOth4","Energy spectrum of muons that reach MRD for configuration 4",20,0,3000);
  TH1D *CCMEC4 = new TH1D("CCMEC4","CCMEC",20,0,3000);
  TH1D *CCQE4 = new TH1D("CCQE4","CCQE",20,0,3000);
  TH1D *NC4 = new TH1D("NC4","NC",20,0,3000);


  // Loop over number of entries
  for(int i=0; i<nentries; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

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
    Double_t* part_xEnd = mTR->get_part_xEnd();
    Double_t* part_yEnd = mTR->get_part_yEnd();
    Double_t* part_zEnd = mTR->get_part_zEnd();
    Double_t* part_pxEnd = mTR->get_part_pxEnd();
    Double_t* part_pyEnd = mTR->get_part_pyEnd();
    Double_t* part_pzEnd = mTR->get_part_pzEnd();
    Double_t* part_KEend = mTR->get_part_KEend();
    Double_t neutrino_E = mTR->get_genE();
    Int_t* mmode = mTR->get_genmode();

    //Intialize a new event
    mTW->InitializeEvent();
    //Add all of the information from the event, except for the photon hits
    mTW->AddWholeBranches(mTR,0,1,1,1,1);


    //Loop over muons to know how many pass thru MRD and to plot their energy spectrum for each configuration
    for(int l=0; l<npart; l++){
      if(part_processEnd[l]==0 && fabs(part_pid[l])==13){
	NpartLeaving = NpartLeaving+1;

	for(int barrelConfig=1; barrelConfig<5; barrelConfig++){

	  double partEnergy = 105.7+part_KEend[l];
	  double ratio = (zEndBarrel[barrelConfig]-part_zEnd[l])/fabs(part_pzEnd[l]);
	  double posfTheta = atan2(part_yEnd[l]+ratio*(part_pyEnd[l]),part_xEnd[l]+ratio*(part_pxEnd[l]));
          double posfR = sqrt((part_xEnd[l]+ratio*(part_pxEnd[l]))*(part_xEnd[l]+ratio*(part_pxEnd[l]))
                        +(part_yEnd[l]+ratio*(part_pyEnd[l]))*((part_yEnd[l])+ratio*(part_pyEnd[l])));

	  if(fabs(posfTheta)>1.0472 && fabs(posfTheta)<2.0944){
	    posfTheta = posfTheta-1.0472;
	  }
	  if(fabs(posfTheta)>2.0944 && fabs(posfTheta)<3.1416){
	    posfTheta = posfTheta-2.0944;
	  }

	  if((part_zEnd[l]<zEndBarrel[barrelConfig] && posfR<-649.*sin(fabs(posfTheta))+Rdet*(1.118)+500. 
		&& part_pzEnd[l]>0) || (part_zEnd[l]>=zEndBarrel[barrelConfig] 
		&& ((posfR>-649.*sin(fabs(posfTheta))+Rdet*(1.118) && part_pzEnd[l]<0) || part_pzEnd[l]>=0))){
	    NpartThruMRD[barrelConfig] = NpartThruMRD[barrelConfig]+1;

	    if(barrelConfig==1){
              if(mmode==1) {
		CCQE1->Fill(partEnergy);
		ModeConf[0][0]++;
	      }
	      if(mmode==2) {
		CCMEC1->Fill(partEnergy);
                ModeConf[1][0]++;
              }
	      if(mmode>=3 && mmode<=30) {
                CCOth1->Fill(partEnergy);
                ModeConf[2][0]++;
              }
	      if(mmode>=31 && mmode<=60) {
                NC1->Fill(partEnergy);
                ModeConf[3][0]++;
              }
	    }
	    if(barrelConfig==2){
              if(mmode==1) {
                CCQE2->Fill(partEnergy);
                ModeConf[0][1]++;
              }
              if(mmode==2) {
                CCMEC2->Fill(partEnergy);
                ModeConf[1][1]++;
              }
              if(mmode>=11 && mmode<=30) {
                CCOth2->Fill(partEnergy);
                ModeConf[2][1]++;
              }
              if(mmode>=31 && mmode<=60) {
                NC2->Fill(partEnergy);
                ModeConf[3][1]++;
              }
	    }
	    if(barrelConfig==3){
              if(mmode==1) {
                CCQE3->Fill(partEnergy);
                ModeConf[0][2]++;
              }
              if(mmode==2) {
                CCMEC3->Fill(partEnergy);
                ModeConf[1][2]++;
              }
              if(mmode>=3 && mmode<=30) {
                CCOth3->Fill(partEnergy);
                ModeConf[2][2]++;
              }
              if(mmode>=31 && mmode<=60) {
                NC3->Fill(partEnergy);
                ModeConf[3][2]++;
              }
	    }
	    if(barrelConfig==4){
              if(mmode==1) {
                CCQE4->Fill(partEnergy);
                ModeConf[0][3]++;
              }
              if(mmode==2) {
                CCMEC4->Fill(partEnergy);
                ModeConf[1][3]++;
              }
              if(mmode>=3 && mmode<=30) {
                CCOth4->Fill(partEnergy);
                ModeConf[2][3]++;
              }
              if(mmode>=31 && mmode<=60) {
                NC4->Fill(partEnergy);
                ModeConf[3][3]++;
              }
	    }
	  }
        }
      }
    }

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


      if(opt==1){
	Nrow = NbyN;
	Ncol = NbyN;
	if(xE==Ldet/2.) mcase=1;
	if(xE==-Ldet/2.) mcase=2;
	if(yE==Ldet/2.) mcase=3;
	if(yE==-Ldet/2.) mcase=4;
	if(zE==Ldet/2.) mcase=5;
	if(zE==-Ldet/2.) mcase=6;
      }
      else {
        if(zE==Ldet/2.){
	  Nrow = NFront;
	  Ncol = NFront;
	  mcase=5;
        }
        if(zE==-Ldet/2.){
          Nrow = NBack;
          Ncol = NBack;
	  mcase=6;
        }
        if((xE*xE + yE*yE)==(Rdet*Rdet)){
	  Nrow = NrowCurve;
          Ncol = NcolCurve;
	  mcase=0;
        }
      }

      int WhichDet = sbPMT->isActiveHit(xE, yE, zE, PMTid, LAPPDs); // 0=none, 1=PMT, 2=LAPPD

      //if the photon lands on a PMT
      if( WhichDet==1 ){
        double tRes = tResPMT;
      }
      //if the photon lands on a LAPPD
      if( WhichDet==2){
	double tRes = tResLAPPD;
      }
        //      cout<<"is a hit"<<endl;
     //if the photon lands on a PMT or LAPPD
      if( WhichDet!=0){
        ishit=2;
        WCSimTrueLight* mTL = new WCSimTrueLight(xS,yS,zS,tS,xE,yE,zE,tE,wl,tRes
					,Rdet,Ldet,isscat,process,parentid,trackid
					,ishit,capnum,PMTid,Nrow,Ncol,mcase,WhichDet,opt);
        mTW->AddPhoton(mTL);
      }
    }

    // Fill the event
    mTW->FillEvent();

  }

  std::cout<<"The number of muons which leave the tank is "<<NpartLeaving<<std::endl;
  std::cout<<"The number of muons which pass thru the MRD in configuration 0 and 1 is "<<NpartThruMRD[1]<<" -> "<<(NpartThruMRD[1]*1.0)/(NpartLeaving*1.0)*100.<<"%"<<std::endl;
  std::cout<<"The number of muons which pass thru the MRD in configuration 2 is "<<NpartThruMRD[2]<<" -> "<<(NpartThruMRD[2]*1.0)/(NpartLeaving*1.0)*100.<<"%"<<std::endl;
  std::cout<<"The number of muons which pass thru the MRD in configuration 3 is "<<NpartThruMRD[3]<<" -> "<<(NpartThruMRD[3]*1.0)/(NpartLeaving*1.0)*100.<<"%"<<std::endl;
  std::cout<<"The number of muons which pass thru the MRD in configuration 4 is "<<NpartThruMRD[4]<<" -> "<<(NpartThruMRD[4]*1.0)/(NpartLeaving*1.0)*100.<<"%"<<std::endl;

  std::cout<<"CCQE & conf 1 ->"<<ModeConf[0][0]<<std::endl;
  std::cout<<"CCQE & conf 2 ->"<<ModeConf[0][1]<<std::endl;
  std::cout<<"CCQE & conf 3 ->"<<ModeConf[0][2]<<std::endl;
  std::cout<<"CCQE & conf 4 ->"<<ModeConf[0][3]<<std::endl;
  std::cout<<"CCMEC & conf 1 ->"<<ModeConf[1][0]<<std::endl;
  std::cout<<"CCMEC & conf 2 ->"<<ModeConf[1][1]<<std::endl;
  std::cout<<"CCMEC & conf 3 ->"<<ModeConf[1][2]<<std::endl;
  std::cout<<"CCMEC & conf 4 ->"<<ModeConf[1][3]<<std::endl;
  std::cout<<"CCOth & conf 1 ->"<<ModeConf[2][0]<<std::endl;
  std::cout<<"CCOth & conf 2 ->"<<ModeConf[2][1]<<std::endl;
  std::cout<<"CCOth & conf 3 ->"<<ModeConf[2][2]<<std::endl;
  std::cout<<"CCOth & conf 4 ->"<<ModeConf[2][3]<<std::endl;
  std::cout<<"NC & conf 1 ->"<<ModeConf[3][0]<<std::endl;
  std::cout<<"NC & conf 2 ->"<<ModeConf[3][1]<<std::endl;
  std::cout<<"NC & conf 3 ->"<<ModeConf[3][2]<<std::endl;
  std::cout<<"NC & conf 4 ->"<<ModeConf[3][3]<<std::endl;

  new TCanvas("conf1","conf1");
  CCOth1->SetLineColor(1);
  CCOth1->Draw();
  CCOth1->GetXaxis()->SetTitle("Energy of muons");
  CCMEC1->SetLineColor(2);
  CCMEC1->Draw("same");
  CCQE1->SetLineColor(3);
  CCQE1->Draw("same");
  NC1->SetLineColor(4);
  NC1->Draw("same");

  TLegend *legend = new TLegend(0.45,0.55,0.76,0.82);
  legend->AddEntry(CCQE1,"CCQE","l");
  legend->AddEntry(CCMEC1,"CCMEC","l");
  legend->AddEntry(CCOth1,"CCOth","l");
  legend->AddEntry(NC1,"NC","l");
  legend->Draw();

  new TCanvas("conf2","conf2");
  CCOth2->SetLineColor(1);
  CCOth2->Draw();
  CCOth2->GetXaxis()->SetTitle("Energy of muons");
  CCMEC2->SetLineColor(2);
  CCMEC2->Draw("same");
  CCQE2->SetLineColor(3);
  CCQE2->Draw("same");
  NC2->SetLineColor(4);
  NC2->Draw("same");
  legend->Draw();

  new TCanvas("conf3","conf3");
  CCOth3->SetLineColor(1);
  CCOth3->Draw();
  CCOth3->GetXaxis()->SetTitle("Energy of muons");
  CCMEC3->SetLineColor(2);
  CCMEC3->Draw("same");
  CCQE3->SetLineColor(3);
  CCQE3->Draw("same");
  NC3->SetLineColor(4);
  NC3->Draw("same");
  legend->Draw();

  new TCanvas("conf4","conf4");
  CCOth4->SetLineColor(1);
  CCOth4->Draw();
  CCOth4->GetXaxis()->SetTitle("Energy of muons");
  CCMEC4->SetLineColor(2);
  CCMEC4->Draw("same");
  CCQE4->SetLineColor(3);
  CCQE4->Draw("same");
  NC4->SetLineColor(4);
  NC4->Draw("same");
  legend->Draw();

  mTW->WriteTreeToFile();
}  
