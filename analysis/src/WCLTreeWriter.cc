#include "WCLTreeWriter.hh"

#include "TF1.h"
#include "TRandom.h"
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <cassert>

ClassImp(WCLTreeWriter)


WCLTreeWriter::WCLTreeWriter(TString filename, int mode)
{

  fmode=mode;

   if(mode==1){ 
    outfile = new TFile(filename,"RECREATE");
    fChain = new TTree("SmearedEventTree","outputEventTree");
    this->SetOutBranchAddys();
  }
  if(mode==0){
    outfile = new TFile(filename,"RECREATE");
    fChain = new TTree("Hits_Tree","Hits_Tree");
    truetree = new TTree("Events_Tree","Events_Tree");
    this->SetFlatBranchAddys();
  }

  eventcount=0;
  nphot=0; nhits=0; ncapturecount=0; npart=0; neutroncount=0;
  fNewEvent=false;
}

WCLTreeWriter::~WCLTreeWriter()
{
  outfile->Close();
}

void WCLTreeWriter::InitializeEvent()
{
  fNewEvent=true;
  nphot=0; nhits=0; ncapturecount=0; npart=0; neutroncount=0;
}


void WCLTreeWriter::AddPhoton(WCSimTrueLight* fPhot){
  phot_xStart[nphot] = fPhot->GetXstart();
  phot_yStart[nphot] = fPhot->GetYstart();
  phot_zStart[nphot] = fPhot->GetZstart();
  phot_tStart[nphot] = fPhot->GetTstart();
  phot_xEnd[nphot] = fPhot->GetXend();
  phot_yEnd[nphot] = fPhot->GetYend();
  phot_zEnd[nphot] = fPhot->GetZend();
  phot_tEnd[nphot] = fPhot->GetTend();
  phot_wavelength[nphot] = fPhot->GetWavelength();
  phot_processStart[nphot] = fPhot->GetProcessStart();
  phot_isScat[nphot] = fPhot->GetIsScat();
  phot_parentid[nphot] = fPhot->GetParentID();
  phot_trackid[nphot] = fPhot->GetTrackID();
  phot_hit[nphot] = fPhot->GetIsHit();
  phot_capnum[nphot] = fPhot->GetCaptureNum();
  phot_PMTid[nphot] = fPhot->GetPMTID();
  nphot++;

}


void WCLTreeWriter::InputEvent(WCLTreeReader* fEvent)
{
  // this function automatically copies the entire contents of an event
  // and then writes that full event to the output tree
  this->InputEventContents(fEvent,1,1,1,1,1,1);
}


void WCLTreeWriter::AddWholeBranches(WCLTreeReader* fEvent, int addphot, int addpart, int addcapt, int addmrd, int addgen)
{
  // this function copies the contents of entire sections of an event
  // but does not write the event to the tree yet
  this->InputEventContents(fEvent,0,addphot,addpart,addcapt,addmrd,addgen);
}
 

  
void WCLTreeWriter::InputEventContents(WCLTreeReader* fEvent, int whole_event, int addphot, int addpart, int addcapt, int addmrd, int addgen)
{
  // this private function is called by InputEvent and AddWholeBranches
  // this copies over the contents of selected parts of the event
  // and write the event to the output tree, if requested.

  if(addphot==1){
    
    nphot = fEvent->get_nphot();
    nhits = fEvent->get_nhits();
    Double_t *pxs = fEvent->get_phot_xStart();
    Double_t *pys = fEvent->get_phot_yStart();
    Double_t *pzs = fEvent->get_phot_zStart();
    Double_t *pts = fEvent->get_phot_tStart();
    Double_t *pxe = fEvent->get_phot_xEnd();
    Double_t *pye = fEvent->get_phot_yEnd();
    Double_t *pze = fEvent->get_phot_zEnd();
    Double_t *pte = fEvent->get_phot_tEnd();
    Double_t* pwl = fEvent->get_phot_wavelength();
    Int_t* pps = fEvent->get_phot_processStart();
    Int_t* pis = fEvent->get_phot_isScat(); 
    Int_t* ppid  = fEvent->get_phot_parentid();
    Int_t* ptid  = fEvent->get_phot_trackid();
    Int_t* ph  = fEvent->get_phot_hit();
    Int_t* pcn  = fEvent->get_phot_capnum();
    Int_t* ppmt  = fEvent->get_phot_PMTid();

    for(int i=0; i<1000000; i++){
      if(i<nphot){
	phot_xStart[i] = pxs[i];
	phot_yStart[i] = pys[i];
	phot_zStart[i] = pzs[i];
	phot_tStart[i] = pts[i];
	phot_xEnd[i] = pxe[i];
	phot_yEnd[i] = pye[i];
	phot_zEnd[i] = pze[i];
	phot_tEnd[i] = pte[i];
	phot_wavelength[i] = pwl[i];
	phot_processStart[i] = pps[i];
	phot_isScat[i] = pis[i];
	phot_parentid[i] = ppid[i];
	phot_trackid[i] = ptid[i];
	phot_hit[i] = ph[i];
	phot_capnum[i] = pcn[i];
	phot_PMTid[i] = ppmt[i];
      } else{
	phot_xStart[i]=0;
	phot_yStart[i] = 0;
	phot_zStart[i] = 0;
	phot_tStart[i] = 0;
	phot_xEnd[i] = 0;
	phot_yEnd[i] = 0;
	phot_zEnd[i] = 0;
	phot_tEnd[i] = 0;
	phot_wavelength[i] = 0;
	phot_processStart[i] = 0;
	phot_isScat[i] = 0;
	phot_parentid[i] = 0;
	phot_trackid[i] = 0;
	phot_hit[i] = 0;
	phot_capnum[i] = 0;
	phot_PMTid[i] = 0;
      }
    
    }
  }


  if(addpart==1){
    npart = fEvent->get_npart();
    Double_t *paxs = fEvent->get_part_xStart();
    Double_t *pays = fEvent->get_part_yStart();
    Double_t *pazs = fEvent->get_part_zStart();
    Double_t *pats = fEvent->get_part_tStart();
    Double_t *paxe = fEvent->get_part_xEnd();
    Double_t *paye = fEvent->get_part_yEnd();
    Double_t *paze = fEvent->get_part_zEnd();
    Double_t *pate = fEvent->get_part_tEnd();
    Double_t *papxs = fEvent->get_part_pxStart();
    Double_t *papys = fEvent->get_part_pyStart();
    Double_t *papzs = fEvent->get_part_pzStart();
    Double_t *pakes = fEvent->get_part_KEstart();
    Double_t *papxe = fEvent->get_part_pxEnd();
    Double_t *papye = fEvent->get_part_pyEnd();
    Double_t *papze = fEvent->get_part_pzEnd();
    Double_t *pakee = fEvent->get_part_KEend();
    Int_t* paps = fEvent->get_part_processStart();
    Int_t* pape = fEvent->get_part_processEnd();
    Int_t* papid  = fEvent->get_part_parentid();
    Int_t* patid  = fEvent->get_part_trackid();
    Int_t* pappid  = fEvent->get_part_pid();

    for(int i=0; i<10000; i++){
      if(i<npart){
	part_xStart[i] = paxs[i];
	part_yStart[i] = pays[i];
	part_zStart[i] = pazs[i];
	part_tStart[i] = pats[i];
	part_xEnd[i] = paxe[i];
	part_yEnd[i] = paye[i];
	part_zEnd[i] = paze[i];
	part_tEnd[i] = pate[i];
	part_processStart[i] = paps[i];
	part_processEnd[i] = pape[i];
	part_parentid[i] = papid[i];
	part_trackid[i] = patid[i];
	part_pid[i] = pappid[i];
	part_pxStart[i] = papxs[i];
	part_pyStart[i] = papys[i];
	part_pzStart[i] = papzs[i];
	part_KEstart[i] = pakes[i];
	part_pxEnd[i] = papxe[i];
	part_pyEnd[i] = papye[i];
	part_pzEnd[i] = papze[i];
	part_KEend[i] = pakee[i];
      } else{
	part_xStart[i]=0;
	part_yStart[i] = 0;
	part_zStart[i] = 0;
	part_tStart[i] = 0;
	part_xEnd[i] = 0;
	part_yEnd[i] = 0;
	part_zEnd[i] = 0;
	part_tEnd[i] = 0;
	part_processStart[i] = 0;
	part_processEnd[i] = 0;
	part_parentid[i] = 0;
	part_trackid[i] = 0;
	part_pid[i] = 0;
	part_pxStart[i] = 0;
	part_pyStart[i] = 0;
	part_pzStart[i] = 0;
	part_KEstart[i] = 0;
	part_pxEnd[i] = 0;
	part_pyEnd[i] = 0;
	part_pzEnd[i] = 0;
	part_KEend[i] = 0;
      }
    }
  }

  if(addcapt==1){
    ncapturecount=fEvent->get_ncapturecount();
    neutroncount=fEvent->get_neutroncount();
    Double_t* cx = fEvent->get_capt_x();
    Double_t* cy = fEvent->get_capt_y();
    Double_t* cz = fEvent->get_capt_z();
    Double_t* ct = fEvent->get_capt_t0();
    Double_t* ce = fEvent->get_capt_E();
    Int_t* cn = fEvent->get_capt_num();
    Int_t* cnc = fEvent->get_capt_nucleus();
    Int_t* cid = fEvent->get_capt_pid();
    Int_t* cnp = fEvent->get_capt_nphot();
    Int_t* cng = fEvent->get_capt_ngamma();
    
    for(int i=0; i<50; i++){
      if(i<ncapturecount){
	capt_x[i] = cx[i];
	capt_y[i] = cy[i];
	capt_z[i] = cz[i];
	capt_t0[i] = ct[i];
	capt_E[i] = ce[i];
	capt_num[i] = cn[i];
	capt_nucleus[i] = cnc[i];
	capt_pid[i] = cid[i];
	capt_nphot[i] = cnp[i];
	capt_ngamma[i] = cng[i];
      } else{
	capt_x[i] = 0;
	capt_y[i] = 0;
	capt_z[i] = 0;
	capt_t0[i] = 0;
	capt_E[i] = 0;
	capt_num[i] = 0;
	capt_nucleus[i] = 0;
	capt_pid[i] = 0;
	capt_nphot[i] = 0;
	capt_ngamma[i] = 0;
      }
      
    }


  }

  if(addmrd==1){

  }

  if(addgen==1){
    
    ntrks = fEvent->get_genntrks();
    mmode = fEvent->get_genmode();
    mbeam_id = fEvent->get_genbeam_id();
    nneutrons = fEvent->get_gennneutrons();
    mevt = fEvent->get_genevt();
    mvtxx = fEvent->get_genvtxx();
    mvtxy = fEvent->get_genvtxy();
    mvtxz = fEvent->get_genvtxz();
    mE = fEvent->get_genE();
    mbeam_px = fEvent->get_genbeam_px();
    mbeam_py = fEvent->get_genbeam_py();
    mbeam_pz = fEvent->get_genbeam_pz();

    Int_t* gid = fEvent->get_genpid();
    Double_t* gpx = fEvent->get_genpx();
    Double_t* gpy = fEvent->get_genpy();
    Double_t* gpz = fEvent->get_genpz();
    Double_t* gpke = fEvent->get_genKE();
    

    for(int i=0; i<100; i++){
      if(i<ntrks){
	//	std::cout<<"gen level "<<nneutrons<<" "<<gid[i]<<std::endl;
	mpid[i] = gid[i];
	mpx[i] = gpx[i];
	mpy[i] = gpy[i];
	mpz[i] = gpz[i];
	mKE[i] = gpke[i];
      } else{
	mpid[i]= 0;
	mpx[i] = 0;
	mpy[i] = 0;
	mpz[i] = 0;
	mKE[i] = 0;
      }
   
    }
    
  }

 
  if(whole_event==1){
    eventcount++;  
    this->FillEvent();
  }
}
  

void WCLTreeWriter::FillEvent()
{
  
  //  std::cout<<"FillEvent(), nphot: "<<nphot<<std::endl;

  fChain->Fill();
  fNewEvent=false;
  eventcount++;
}


void WCLTreeWriter::WriteTreeToFile()
{
  outfile->cd();
  fChain->Write();
  if(fmode==0) truetree->Write();
}


Int_t WCLTreeWriter::GetEntries()
{
 return  fChain->GetEntries();
}

void WCLTreeWriter::ResetForNewSample()
{
  if( fChain->GetEntries()>0 ){
    std::cout << " *** WCLTreeWriter::Reset() *** " << std::endl;  
  }

}

void WCLTreeWriter::SetOutBranchAddys()
{
  std::cout<<"Setting Branch Addresses"<<std::endl;

  const int knphotmax=1000000;

  phot_xStart = new Double_t[knphotmax];
  phot_yStart = new Double_t[knphotmax];
  phot_zStart = new Double_t[knphotmax];
  phot_tStart = new Double_t[knphotmax];
  phot_xEnd = new Double_t[knphotmax];
  phot_yEnd = new Double_t[knphotmax];
  phot_zEnd = new Double_t[knphotmax];
  phot_tEnd = new Double_t[knphotmax];
  phot_wavelength = new Double_t[knphotmax];
  phot_processStart = new Int_t[knphotmax];
  phot_isScat = new Int_t[knphotmax];
  phot_parentid = new Int_t[knphotmax];
  phot_trackid = new Int_t[knphotmax];
  phot_hit = new Int_t[knphotmax];
  phot_capnum = new Int_t[knphotmax];
  phot_PMTid = new int[knphotmax];

  fChain->Branch("evt",&eventcount);
  fChain->Branch("nphot",&nphot);
  fChain->Branch("nhits",&nhits);
  fChain->Branch("npart",&npart);
   
  fChain->Branch("phot_xStart",phot_xStart,"phot_xStart[nphot]/D");
  fChain->Branch("phot_yStart",phot_yStart,"phot_yStart[nphot]/D");
  fChain->Branch("phot_zStart",phot_zStart,"phot_zStart[nphot]/D");
  fChain->Branch("phot_tStart",phot_tStart,"phot_tStart[nphot]/D");
  fChain->Branch("phot_xEnd",phot_xEnd,"phot_xEnd[nphot]/D");
  fChain->Branch("phot_yEnd",phot_yEnd,"phot_yEnd[nphot]/D");
  fChain->Branch("phot_zEnd",phot_zEnd,"phot_zEnd[nphot]/D");
  fChain->Branch("phot_tEnd",phot_tEnd,"phot_tEnd[nphot]/D");
  fChain->Branch("phot_wavelength",phot_wavelength,"phot_wavelength[nphot]/I");
  fChain->Branch("phot_processStart",phot_processStart,"phot_processStart[nphot]/I");
  fChain->Branch("phot_isScat",phot_isScat,"phot_isScat[nphot]/I");
  fChain->Branch("phot_parentid",phot_parentid,"phot_parentid[nphot]/I");
  fChain->Branch("phot_trackid",phot_trackid,"phot_trackid[nphot]/I");
  fChain->Branch("phot_hit",phot_hit,"phot_hit[nphot]/I");
  fChain->Branch("phot_PMTid",phot_PMTid,"phot_PMTid[nphot]/I");
  fChain->Branch("phot_capnum",phot_capnum,"phot_capnum[nphot]/I");

  const int knpartmax = 10000;
  part_xStart = new double[knpartmax];
  part_yStart = new double[knpartmax];
  part_zStart = new double[knpartmax];
  part_tStart = new double[knpartmax];
  part_xEnd = new double[knpartmax];
  part_yEnd = new double[knpartmax];
  part_zEnd = new double[knpartmax];
  part_tEnd = new double[knpartmax];
  part_KEstart = new double[knpartmax];
  part_KEend = new double[knpartmax];
  part_pxStart = new double[knpartmax];
  part_pyStart = new double[knpartmax];
  part_pzStart = new double[knpartmax];
  part_pxEnd = new double[knpartmax];
  part_pyEnd = new double[knpartmax];
  part_pzEnd = new double[knpartmax];
  part_processStart = new int[knpartmax];
  part_processEnd = new int[knpartmax];
  part_parentid = new int[knpartmax];
  part_trackid = new int[knpartmax];
  part_pid = new int[knpartmax];
  
  fChain->Branch("part_xStart",part_xStart,"part_xStart[npart]/D");
  fChain->Branch("part_yStart",part_yStart,"part_yStart[npart]/D");
  fChain->Branch("part_zStart",part_zStart,"part_zStart[npart]/D");
  fChain->Branch("part_tStart",part_tStart,"part_tStart[npart]/D");
  fChain->Branch("part_xEnd",part_xEnd,"part_xEnd[npart]/D");
  fChain->Branch("part_yEnd",part_yEnd,"part_yEnd[npart]/D");
  fChain->Branch("part_zEnd",part_zEnd,"part_zEnd[npart]/D");
  fChain->Branch("part_tEnd",part_tEnd,"part_tEnd[npart]/D");
  fChain->Branch("part_pxStart",part_pxStart,"part_pxStart[npart]/D");
  fChain->Branch("part_pyStart",part_pyStart,"part_pyStart[npart]/D");
  fChain->Branch("part_pzStart",part_pzStart,"part_pzStart[npart]/D");
  fChain->Branch("part_pxEnd",part_pxEnd,"part_pxEnd[npart]/D");
  fChain->Branch("part_pyEnd",part_pyEnd,"part_pyEnd[npart]/D");
  fChain->Branch("part_pzEnd",part_pzEnd,"part_pzEnd[npart]/D");
  fChain->Branch("part_KEstart",part_KEstart,"part_KEstart[npart]/D");
  fChain->Branch("part_KEend",part_KEend,"part_KEend[npart]/D");
  fChain->Branch("part_processStart",part_processStart,"part_processStart[npart]/I");
  fChain->Branch("part_processEnd",part_processEnd,"part_processEnd[npart]/I");
  fChain->Branch("part_parentid",part_parentid,"part_parentid[npart]/I");
  fChain->Branch("part_trackid",part_trackid,"part_trackid[npart]/I");
  fChain->Branch("part_pid",part_pid,"part_pid[npart]/I");

  const int kcapmax = 50;
  capt_num = new int[kcapmax];
  capt_nucleus = new int[kcapmax];
  capt_pid = new int[kcapmax];
  capt_nphot = new int[kcapmax];
  capt_ngamma = new int[kcapmax];
  capt_x = new double[kcapmax];
  capt_y = new double[kcapmax];
  capt_z = new double[kcapmax];
  capt_t0 = new double[kcapmax];
  capt_E = new double[kcapmax];

  fChain->Branch("ncapturecount",&ncapturecount);
  fChain->Branch("neutroncount",&neutroncount);
  fChain->Branch("capt_x",capt_x,"capt_x[ncapturecount]/D");
  fChain->Branch("capt_y",capt_y,"capt_y[ncapturecount]/D");
  fChain->Branch("capt_z",capt_z,"capt_z[ncapturecount]/D");
  fChain->Branch("capt_t0",capt_t0,"capt_t0[ncapturecount]/D");
  fChain->Branch("capt_E",capt_E,"capt_E[ncapturecount]/D");
  fChain->Branch("capt_num",capt_num,"capt_num[ncapturecount]/I");
  fChain->Branch("capt_pid",capt_pid,"capt_pid[ncapturecount]/I");
  fChain->Branch("capt_nucleus",capt_nucleus,"capt_nucleus[ncapturecount]/I");
  fChain->Branch("capt_nphot",capt_nphot,"capt_nphot[ncapturecount]/I");
  fChain->Branch("capt_ngamma",capt_ngamma,"capt_ngamma[ncapturecount]/I");

  const int kmaxtrk=100;
  const int kmaxMRDhits=50;
  mpid = new int[kmaxtrk];
  mpx = new double[kmaxtrk];
  mpy = new double[kmaxtrk];
  mpz = new double[kmaxtrk];
  mKE = new double[kmaxtrk];

  fChain->Branch("mode",&mmode);
  fChain->Branch("neutrino_E",&mE);
  fChain->Branch("neutrino_id",&mbeam_id);
  fChain->Branch("neutrino_px",&mbeam_px);
  fChain->Branch("neutrino_py",&mbeam_py);
  fChain->Branch("neutrino_pz",&mbeam_pz);
  fChain->Branch("ntrks",&ntrks);
  fChain->Branch("nneutrons",&nneutrons);
  fChain->Branch("vtxx",&mvtxx); 
  fChain->Branch("vtxy",&mvtxy); 
  fChain->Branch("vtxz",&mvtxz);
  fChain->Branch("mpid",mpid,"mpid[ntrks]/I");
  fChain->Branch("px",mpx,"px[ntrks]/D"); 
  fChain->Branch("py",mpy,"py[ntrks]/D"); 
  fChain->Branch("pz",mpz,"pz[ntrks]/D");
  fChain->Branch("KE",mKE,"KE[ntrks]/D");

}



void WCLTreeWriter::SetFlatBranchAddys(){


  fChain->Branch("x_hit",&x_hit);
  fChain->Branch("y_hit",&y_hit);
  fChain->Branch("z_hit",&z_hit);
  fChain->Branch("true_time",&true_time);
  fChain->Branch("true_time_corrected",&true_time_corrected);
  fChain->Branch("process",&process);
  fChain->Branch("eventID",&evtID);
  fChain->Branch("eventglobal",&evtID);
  fChain->Branch("photon_wavelength",&photon_wavelength);  

  truetree->Branch("eventID_ET",&evtID);
  truetree->Branch("pos_x",&pos_x);
  truetree->Branch("pos_y",&pos_y);
  truetree->Branch("pos_z",&pos_z);
  truetree->Branch("dir_x1",&dir_x);
  truetree->Branch("dir_y1",&dir_y);
  truetree->Branch("dir_z1",&dir_z);
  truetree->Branch("t",&tt);
  truetree->Branch("Et",&Et);
  truetree->Branch("Ekin1",&Ekin);
  truetree->Branch("pdg",&pidd);
  truetree->Branch("processStep",&mprocessStep);
  truetree->Branch("processStep",&mprocessStep);
  truetree->Branch("parentID",&mparentID);
  truetree->Branch("trackID",&mtrackID);
  truetree->Branch("processStart",&mprocessStart);


}



