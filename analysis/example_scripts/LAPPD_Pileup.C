{
  // This script studies how many photons hit each LAPPD 
  // in different coverage schemes
  // Matt W


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
  TString filename("/lbne/data2/users/wetstein/LiteSimulations/LAPPDcovgstudy/FullEvent_-50cm.root");
  TString genfilename("/lbne/data2/users/wetstein/LiteSimulations/LAPPDcovgstudy/generatorcardfile_-50cm.root");

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
  SandBoxPMTcoverage* sbPMT = new SandBoxPMTcoverage();
  // set the dimensions of the detector -x,+x,-y,+y,-z,+z
  sbPMT->SetBoxDimensions(-1500.,1500.,-1500.,1500.,-1500.,1500.);
  //configure each wall
 

  // front-wall, PMTs, square-shaped, 203mm, 4x4 grid 
  sbPMT->SetWallConfiguration(0,1,1,203.2,6);
  sbPMT->SetWallConfiguration(1,1,1,203.2,6);
  sbPMT->SetWallConfiguration(2,1,1,203.2,6);
  sbPMT->SetWallConfiguration(3,1,1,203.2,6);
  sbPMT->SetWallConfiguration(4,1,1,203.2,6);
  sbPMT->SetWallConfiguration(5,1,1,203.2,6);

  // create relevant histograms
    
  int nbinsH=160;
  double lowbH=-2.5;
  double hibH=897.5;

  // nearest point of origin of light for all light and with LAPPD cvg
  TH1D* nearestlight_cvg = new TH1D("nl_c","nl_c",nbinsH,lowbH,hibH);
  TH1D* nearestlight_all = new TH1D("nl_a","nl_a",nbinsH,lowbH,hibH);
  // how likely to find light from the track within the specified x
  TH1D* confidence_cvg = new TH1D("cfl_c","cfl_c",nbinsH,lowbH,hibH);
  TH1D* confidence_all = new TH1D("cfl_a","cfl_a",nbinsH,lowbH,hibH);
  // locations of the vertices in x, y, and z
  TH1D* hvtxx = new TH1D("hvtxx","hvtxx",300,-150,150);
  TH1D* hvtxy = new TH1D("hvtxy","hvtxy",300,-150,150);
  TH1D* hvtxz = new TH1D("hvtxz","hvtxz",300,-150,150);
  // fraction of "saturated LAPPDs" per event. Histogram of overall hits-per-LAPPD (many events)
  TH1D* fracsat = new TH1D("fracsat","fracsat",201,-0.25,100.25);
  TH1D* pileuphist = new TH1D("multi-photon-pileup","multi-photon-pileup",101,-0.5,100.5);
  // pile-up per event
  TH1D** pu_evt = new TH1D*[1000];

  for(int pp=0; pp<1000; pp++){
    TString pun;
    pun+="pileup";
    pun+=pp;
    pu_evt[pp] = new TH1D(pun,pun,101,-0.5,100.5);
  }

  int nentries = mTR->GetEntries();
  cout<<nentries<<endl;

  // Loop over number of entries
  for(int i=0; i<nentries; i++){
    if(i%10==0) cout<<"event: "<<i<<endl;

    // event-by-event histograms
    TH1D* hvdist_all = new TH1D("hvdist_all","hvdist_all",nbinsH,lowbH,hibH);
    TH1D* hvdist_cvg = new TH1D("hvdist_cvg","hvdist_cvg",nbinsH,lowbH,hibH);
    TH1D* PMTidhist = new TH1D("pmtidh","pmtidh",601,-0.5,600.5);

    // Load the event variables into TreeReader
    mTR->LoadEvent(i);

    // generator-level vertices
    double mvtxx = mTR->get_genvtxx();
    double mvtxy = mTR->get_genvtxy();
    double mvtxz = mTR->get_genvtxz();
    hvtxx->Fill(mvtxx);
    hvtxy->Fill(mvtxy);
    hvtxz->Fill(mvtxz);

    //calculate distance in z direction from the center of the MRD
    double MRDthickness = 30.;
    double zdist = 150-mvtxz;
    double zp = zdist + (MRDthickness/2.);

    // generator variables
    bool ismuon=false;
    int ngentrks = mTR->get_genntrks();
    int* genpid = mTR->get_genpid();
    double mKE,mpx,mpy,mpz,mthet;
    double slopexz,slopeyz;
    
    // loop over gen tracks
    for(int gi=0; gi<ngentrks; gi++){

      // if there is a muon in the event, make note
      if(genpid[gi]==13){
	ismuon=true;

	// what is the direction and energy of the muon?
	double* genpx = mTR->get_genpx();
	double* genpy = mTR->get_genpy();
	double* genpz = mTR->get_genpz();
	double* genKE = mTR->get_genKE();
	
	if(genpz[gi]!=0) slopexz = genpx[gi]/genpz[gi];
	else slopexz=100.;

	if(genpz[gi]!=0) slopeyz = genpy[gi]/genpz[gi];
	else slopeyz=100.;
	
	// what is the angle off of the beam axis?
	TVector3 vbeam(0,0,1);
	TVector3 vmuon(mpx,mpy,mpz);
	mthet = vbeam.Angle(vmuon); 

      }
    }

    if(ismuon){

      // does the muon fall comfortably in the MRD?
      double mdx = slopexz*zp;
      double mdy = slopeyz*zp;
      bool inMRDangle = false;
      if( (mvtxx+mdx<145)&&(mvtxx+mdx>-145)
	  &&(mvtxy+mdy<145)&&(mvtxy+mdy>-145) ) inMRDangle=true;

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
	
	//calculate the distance of the phot_start from the vertex
	double dxs = (phot_xStart[k] - 10*mvtxx)*(phot_xStart[k] - 10*mvtxx);
	double dys = (phot_yStart[k] - 10*mvtxy)*(phot_yStart[k] - 10*mvtxy);
	double dzs = (phot_zStart[k] - 10*mvtxz)*(phot_zStart[k] - 10*mvtxz);
	double vtxdistance = sqrt(dxs+dys+dzs);

	if(phot_isScat[k]==0) hvdist_all->Fill(vtxdistance);
	
	//if the photon lands on an LAPPD, according to the PMTcoverage class
	if( sbPMT->isActiveHitLAPPD(xE, yE, zE, PMTid) ){
	  
	  ishit=2;
	  WCSimTrueLight* mTL = new WCSimTrueLight(xS,yS,zS,tS,xE,yE,zE,tE,wl,
						   process,isscat,parentid,trackid,ishit,capnum,PMTid);

	  // prompt light (<50 nseconds after event start)
	  if(tS<50.){
	    PMTidhist->Fill(PMTid);
	  }

	  mTW->AddPhoton(mTL);

	  
	  if((phot_isScat[k]==0) && inMRDangle) hvdist_cvg->Fill(vtxdistance);

	}
      }
      
      // Fill pile-up hist
      for(int gg=1; gg<602; gg++){

	double mbc = PMTidhist->GetBinContent(gg);
	if(mbc>0){

	  pileuphist->Fill(mbc/2.);
	  pu_evt[i]->Fill(mbc/2.);
	}

      }
      
      // fill nearest light hist
      bool firstbin_all = true;
      bool firstbin_cvg = true;
      for(int mm=1; mm<nbinsH+1; mm++){
	
	double bc_all = hvdist_all->GetBinContent(mm);
	double bc_cvg = hvdist_cvg->GetBinContent(mm);

	if( (bc_cvg>10) && (firstbin_cvg) ){
	  double bcenter = hvdist_cvg->GetBinCenter(mm);
	  firstbin_cvg=false;
	  nearestlight_cvg->Fill(bcenter);
	}

	if( (bc_all>10) && (firstbin_all) ){
	  double bcenter = hvdist_all->GetBinCenter(mm);
	  firstbin_all=false;
	  nearestlight_all->Fill(bcenter);
	}

	
      }
            
      // Fill the event
      mTW->FillEvent();
    }
    delete PMTidhist;
    delete hvdist_all;
    delete hvdist_cvg;
  }
  mTW->WriteTreeToFile();

  /*
  nearestlight_all->Draw();
  nearestlight_cvg->SetLineColor(2);
  nearestlight_cvg->Draw("SAME");
  */

  // Normalize and write histograms

  int nbnl = nearestlight_all->GetNbinsX();
  double bsum=0;
  for(int j=0; j<nbnl; j++){
    bsum+=(nearestlight_all->GetBinContent(j));
    confidence_all->SetBinContent(j,bsum);	     
  }
  if(bsum>0) confidence_all->Scale(1/bsum);

  nbnl = nearestlight_cvg->GetNbinsX();
  bsum=0;
  for(int j=0; j<nbnl; j++){
    bsum+=(nearestlight_cvg->GetBinContent(j));
    confidence_cvg->SetBinContent(j,bsum);	     
  }
  if(bsum>0) confidence_cvg->Scale(1/bsum);

  TCanvas* cvgcan = new TCanvas();
  confidence_all->Draw();
  confidence_cvg->SetLineColor(2);
  confidence_cvg->Draw("SAME");

  TFile* mmm = new TFile("pileupana.root","RECREATE");
  nearestlight_cvg->Write("nearestlight_cvg");
  nearestlight_all->Write("nearestlight_all");
  confidence_cvg->Write("confidence_cvg");
  confidence_all->Write("confidence_all");
  pileuphist->Write("pileuphist");

  double pdsum=0;
  double pdsat=0;
  for(int qq=0; qq<1000; qq++){
    pu_evt[qq]->Write();

    for(int bb=1; bb<101; bb++){
      pdsum+=(pu_evt[qq]->GetBinContent(bb));
      double pdbcent = pu_evt[qq]->GetBinCenter(bb);
      if(pdbcent>30) pdsat+=(pu_evt[qq]->GetBinContent(bb));
    }
  
    double fsat;
    if(pdsum<=0) fsat=0;
    else fsat = pdsat/pdsum;
    fracsat->Fill(fsat*100);

  }
  fracsat->Write();
  mmm->Close();

  //The End
}  

