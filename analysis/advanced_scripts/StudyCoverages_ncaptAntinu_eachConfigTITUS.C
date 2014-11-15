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
  TString filename("/vols/t2k03/hyperk/output/sim/genie/genie/FullEvent_antinu_*.root");
  TString genfilename("/vols/t2k03/hyperk/output/sim/genie/genie/generatorcardfile_antinu_*.root");

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
  int config = 2;

  TH1D *nLightAntinu1 = new TH1D("nLightAntinu1","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu2 = new TH1D("nLightAntinu2","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu3 = new TH1D("nLightAntinu3","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu4 = new TH1D("nLightAntinu4","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu5 = new TH1D("nLightAntinu5","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu6 = new TH1D("nLightAntinu6","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu7 = new TH1D("nLightAntinu7","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu8 = new TH1D("nLightAntinu8","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);
  TH1D *nLightAntinu9 = new TH1D("nLightAntinu9","Number of photons which reach a photosensor for a single Gd-captured neutron.",50,0,500);

  for(int config = 1; config<10; config++){  
    std::cout<<"Configuration number "<<config<<std::endl;

    // set if you use NEUT (whichFile=0) or GENIE (whichFile=1) file
    int whichFile = 1;

    // set the presence of LAPPDs
    bool LAPPDs = false;

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
    // double the QE if true because WChSandbox divide by 2 the number of photons
    bool QEset = true;
    // set the size of the detectors
    double sizePMT = 304.8;
    double sizeLAPPD = 203.2;
    double Sdet;
    // set the PMT coverage
    double cover = 40.;
    // set the shape of the detectors: 0=circular, 1=square
    int shapePMT = 0;
    int shapeLAPPD = 1;

    // Number of PMTs in each wall for the option=0.
    int NFront = 34;
    int NBack = 34;
    int NrowCurve = 49;
    int NcolCurve = 78;
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

    if(QEset){
      QE_PMT = QE_PMT*2.;
      QE_LAPPD = QE_LAPPD*2.;
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
      Int_t ncapturecount = mTR->get_ncapturecount();
      Int_t mmode = mTR->get_genmode();
      Double_t* phot_xStart = mTR->get_phot_xStart();
      Double_t* phot_yStart = mTR->get_phot_yStart();
      Double_t* phot_zStart = mTR->get_phot_zStart();
      Double_t* phot_tStart = mTR->get_phot_tStart();
      Double_t* phot_xEnd = mTR->get_phot_xEnd();
      Double_t* phot_yEnd = mTR->get_phot_yEnd();
      Double_t* phot_zEnd = mTR->get_phot_zEnd();
      Double_t* phot_tEnd = mTR->get_phot_tEnd();
      Double_t* phot_wavelength = mTR->get_phot_wavelength();
      Double_t* capt_x = mTR->get_capt_x();
      Double_t* capt_y = mTR->get_capt_y();
      Double_t* capt_z = mTR->get_capt_z();
      Double_t* capt_t0 = mTR->get_capt_t0();
      Double_t* capt_E = mTR->get_capt_E();

      //Intialize a new event
      mTW->InitializeEvent();
      //Add all of the information from the event, except for the photon hits
      mTW->AddWholeBranches(mTR,0,1,1,1,1);

      int nphotbyNeutron = 0;
      if(ncapturecount==1 && capt_E[ncapturecount-1]>7.7 && capt_E[ncapturecount-1]<8.7){

        //Loop over photons to add those which hit a PMT or LAPPD
        for(int k=0; k<nphot; k++){
//          if( (k%10000==0) && (i%10==0) ) cout<<k<<" photons out of "<<nphot<<endl;  
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

          double t_flight = tE-capt_t0[ncapturecount-1];
          double distSquare = (xE-capt_x[ncapturecount-1])*(xE-capt_x[ncapturecount-1])
			 + (yE-capt_y[ncapturecount-1])*(yE-capt_y[ncapturecount-1])
			 + (zE-capt_z[ncapturecount-1])*(zE-capt_z[ncapturecount-1]);

          if( t_flight < (sqrt(distSquare)/(300./1.33))+50
		 && t_flight > (sqrt(distSquare)/(300./1.33))-5 && t_flight > 0. ){


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
            if(WhichDet==1){
              double tRes = tResPMT;
              nphotbyNeutron++;
            }
            //if the photon lands on a LAPPD
            if(WhichDet==2){
	      double tRes = tResLAPPD;
              nphotbyNeutron++;
            }
            //      cout<<"is a hit"<<endl;
            //if the photon lands on a PMT or LAPPD
            if(WhichDet!=0){
              ishit=2;
              WCSimTrueLight* mTL = new WCSimTrueLight(xS,yS,zS,tS,xE,yE,zE,tE,wl,tRes,Rdet,Ldet,isscat
				,process,parentid,trackid,ishit,capnum,PMTid,Nrow,Ncol,mcase,WhichDet,opt);
              mTW->AddPhoton(mTL);
            }
          }
        }
      }


      if (nphotbyNeutron != 0){
        if (config==1) nLightAntinu1->Fill(nphotbyNeutron);
        if (config==2) nLightAntinu2->Fill(nphotbyNeutron);
        if (config==3) nLightAntinu3->Fill(nphotbyNeutron);
        if (config==4) nLightAntinu4->Fill(nphotbyNeutron);
        if (config==5) nLightAntinu5->Fill(nphotbyNeutron);
        if (config==6) nLightAntinu6->Fill(nphotbyNeutron);
        if (config==7) nLightAntinu7->Fill(nphotbyNeutron);
        if (config==8) nLightAntinu8->Fill(nphotbyNeutron);
        if (config==9) nLightAntinu9->Fill(nphotbyNeutron);
      }

      // Fill the event
      mTW->FillEvent();
    }
  } 

  new TCanvas("canvas12","canvas12");
  nLightAntinu1->SetLineWidth(2);
  nLightAntinu2->SetLineWidth(2);
  nLightAntinu2->SetLineStyle(5);
  nLightAntinu3->SetLineWidth(2);
  nLightAntinu3->SetLineStyle(3);
  nLightAntinu4->SetLineWidth(2);
  nLightAntinu4->SetLineColor(2);
  nLightAntinu5->SetLineWidth(2);
  nLightAntinu5->SetLineColor(2);
  nLightAntinu5->SetLineStyle(5);
  nLightAntinu6->SetLineWidth(2);
  nLightAntinu6->SetLineColor(2);
  nLightAntinu6->SetLineStyle(3);
  nLightAntinu7->SetLineWidth(2);
  nLightAntinu7->SetLineColor(4);
  nLightAntinu8->SetLineWidth(2);
  nLightAntinu8->SetLineColor(4);
  nLightAntinu8->SetLineStyle(5);
  nLightAntinu9->SetLineWidth(2);
  nLightAntinu9->SetLineColor(4);
  nLightAntinu9->SetLineStyle(3);

  nLightAntinu7->GetXaxis()->SetTitle("Number of photons");
  nLightAntinu7->GetYaxis()->SetTitle("Number of events");
  nLightAntinu7->Draw();
  nLightAntinu1->Draw("same");
  nLightAntinu2->Draw("same");
  nLightAntinu3->Draw("same");
  nLightAntinu4->Draw("same");
  nLightAntinu5->Draw("same");
  nLightAntinu6->Draw("same");
  nLightAntinu8->Draw("same");
  nLightAntinu9->Draw("same");

  TLegend *legend = new TLegend(0.70,0.40,0.86,0.86);
  legend->AddEntry(nLightAntinu1,"40% 20\"","l");
  legend->AddEntry(nLightAntinu2,"40% 12\"","l");
  legend->AddEntry(nLightAntinu3,"40% 8\"","l");
  legend->AddEntry(nLightAntinu4,"30% 20\"","l");
  legend->AddEntry(nLightAntinu5,"30% 12\"","l");
  legend->AddEntry(nLightAntinu6,"30% 8\"","l");
  legend->AddEntry(nLightAntinu7,"20% 20\"","l");
  legend->AddEntry(nLightAntinu8,"20% 12\"","l");
  legend->AddEntry(nLightAntinu9,"20% 8\"","l");
  legend->Draw();

  mTW->WriteTreeToFile();
}
