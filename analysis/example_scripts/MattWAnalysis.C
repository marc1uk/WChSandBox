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
  TString filename("../../SmearedEvents_test_6.root");
  //  TString filename("../../SmearedEvents_3k_-50_6by6.root");

  // Load Data
  // =========
  
  // WCLTREEREADER
  //**********************************************************
  // first parameter in the construct: 0=direct output of WCLite
  //                                   1=a "smeared" file made from WCL files
  // second parameter in the constructor: 1 means "include gen level tree"
  //                                      0 means no gen level tree

  WCLTreeReader *mTR = new WCLTreeReader(1,0);

  // LoadData is an overloaded function. If you specify only one name,
  // only the main tree is loaded. If you specify two files, it loads
  // both the geant output file and the file containing the gen tree
  // (in that order)

  cout<<"loading "<<filename<<endl;
  
  mTR->LoadData(filename);

  cout<<"initializing treewriter"<<endl;

  cout<<"getting nentries"<<endl;

  int nentries = mTR->GetEntries();
  
  cout<<nentries<<endl;

  int nbinsH=160;
  double lowbH=-2.5;
  double hibH=897.5;
  TH1D* nearestlight = new TH1D("nl","nl",nbinsH,lowbH,hibH);
  TH1D* confidence = new TH1D("cfl","cfl",nbinsH,lowbH,hibH);

  TH1D* hq = new TH1D("q","q",20,0,2000);
  TH1D* hqmrd = new TH1D("qmrd","qmrd",20,0,2000);

  TH1D* hvtxx = new TH1D("hvtxx","hvtxx",300,-150,150);
  TH1D* hvtxy = new TH1D("hvtxy","hvtxy",300,-150,150);
  TH1D* hvtxz = new TH1D("hvtxz","hvtxz",300,-150,150);


  TH2D** photxystarts = new TH1D*[10];
  TH2D** photxzstarts = new TH1D*[10];
  TH2D** photyzstarts = new TH1D*[10];

  for(int mm=0; mm<10; mm++){

    photxystarts[mm] = new TH1D(nn,nn,300,-150.,150.);
    photxzstarts[mm] = new TH1D(nn,nn,300,-150.,150.);
    photyzstarts[mm] = new TH1D(nn,nn,300,-150.,150.);

  }



  //  for(int i=0; i<nentries; i++){
  nentries=10;
  for(int i=0; i<nentries; i++){

    if(i%100==0) cout<<"event: "<<i<<endl;

    mTR->LoadEvent(i);
    int nphot = mTR->get_nphot();

    double mvtxx = mTR->get_genvtxx();
    double mvtxy = mTR->get_genvtxy();
    double mvtxz = mTR->get_genvtxz();
    
    cout<<mvtxx<<" "<<mvtxy<<" "<<mvtxz<<endl;

    hvtxx->Fill(mvtxx);
    hvtxy->Fill(mvtxy);
    hvtxz->Fill(mvtxz);

    double* phot_xStart = mTR->get_phot_xStart();
    double* phot_yStart = mTR->get_phot_yStart();
    double* phot_zStart = mTR->get_phot_zStart();
    double* phot_isScat = mTR->get_phot_isScat();

    TH1D* hvdist = new TH1D("hvdist","hvdist",nbinsH,lowbH,hibH);
    
    //calculate distance in z direction from the center of the MRD
    double MRDthickness = 30.;
    double zdist = 150-mvtxz;
    double zp = zdist + (MRDthickness/2.);
    // cout<<"distance from approx MRD center: "<<zp<<endl;

    // loop over gen particles. Is there a muon? What's the slope?
    // What's the scattering angle? Momentum transfer?
    bool ismuon=false;
    int ngentrks = mTR->get_genntrks();
    double* genpx = mTR->get_genpx();
    double* genpy = mTR->get_genpy();
    double* genpz = mTR->get_genpz();
    double* genKE = mTR->get_genKE();
    int* genpid = mTR->get_genpid();

    double mKE,mpx,mpy,mpz,mthet;
    for(int gi=0; gi<ngentrks; gi++){
      
      if(genpid[gi]==13){
	ismuon=true;
	mKE = genKE[gi];
	mpx = genpx[gi];
	mpy = genpy[gi];
	mpz = genpz[gi];	  
	
	double slopexz,slopeyz;
	
	if(genpz[gi]!=0) slopexz = genpx[gi]/genpz[gi];
	else slopexz=100.;

	if(genpz[gi]!=0) slopeyz = genpy[gi]/genpz[gi];
	else slopeyz=100.;
	
	TVector3 vbeam(0,0,1);
	TVector3 vmuon(mpx,mpy,mpz);
	mthet = vbeam.Angle(vmuon); 
      }
    }      
    if(ismuon){
      
      double totE = mKE + 105.7;
      double p = sqrt( totE*totE - 105.7*105.7 );
      double qabs = fabs(2*p*sin(mthet/2.));
      hq->Fill(qabs);
      
      // does the muon fall comfortably in the MRD?
      double mdx = slopexz*zp;
      double mdy = slopeyz*zp;
      bool inMRDangle = false;
      if( (mvtxx+mdx<145)&&(mvtxx+mdx>-145)
	  &&(mvtxy+mdy<145)&&(mvtxy+mdy>-145) ) inMRDangle=true;
      
      if(inMRDangle) hqmrd->Fill(qabs);
      
      //cout<<"is muon! "<<qabs<<" "<<mKE<<" "<<mpx<<" "<<mpy<<" "<<mpz<<endl;
      
    } //else cout<<"NO MUON! "<< (mTR->get_genmode()) <<endl;
      
      
    
    for(int j=0; j<nphot; j++){
      
      double dxs = (phot_xStart[j] - 10*mvtxx)*(phot_xStart[j] - 10*mvtxx);
      double dys = (phot_yStart[j] - 10*mvtxy)*(phot_yStart[j] - 10*mvtxy);
      double dzs = (phot_zStart[j] - 10*mvtxz)*(phot_zStart[j] - 10*mvtxz);
      double vtxdistance = sqrt(dxs+dys+dzs);
      
      //      cout<<"vtx distance "<<vtxdistance<<endl;
      if(phot_isScat[j]==0) hvdist->Fill(vtxdistance);
    }
    
    bool firstbin = true;
    for(int mm=1; mm<nbinsH+1; mm++){
      
      double bc = hvdist->GetBinContent(mm);
      nearestlight->SetBinContent(mm,bc);
      /*
      //      if(mm<10) cout<<"mm "<<mm<<" BC "<<bc<<endl;

      if( (bc>2) && (firstbin) ){
	double bcenter = hvdist->GetBinCenter(mm);
	// cout<<mm<<" "<<bcenter<<endl;
	firstbin=false;
	nearestlight->Fill(bcenter);
      }
      */
    }
    //    cout<<"firstbin: "<<firstbin<<endl;

    delete hvdist; 
    
  }


  nearestlight->Draw();
  TCanvas* tc = new TCanvas();
  hq->Draw();
  hqmrd->SetLineColor(2);
  hqmrd->Draw("SAME");

  /*  
  TCanvas* effic = new TCanvas();
  hqeffic = hqmrd->Divide(hq);
  */

  TCanvas* vt = new TCanvas();
  hvtxx->Draw();
  hvtxy->SetLineColor(4);
  hvtxy->Draw("SAME");
  hvtxz->SetLineColor(2);
  hvtxz->Draw("SAME");


  int nbnl = nearestlight->GetNbinsX();
  double bsum=0;
  for(int j=0; j<nbnl; j++){
    bsum+=(nearestlight->GetBinContent(j));
    confidence->SetBinContent(j,bsum);	     
  }
  if(bsum>0) confidence->Scale(1/bsum);
  
  TFile* pp = new TFile("test.root","RECREATE");
  
  hq->Write();
  hqmrd->Write();

  nearestlight->Write();
  confidence->Write();

  pp->Close();
}

