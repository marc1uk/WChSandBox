#include "WCLTreeReader.hh"

#include "TF1.h"
#include "TRandom.h"
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <cassert>

#include "WCSimTrueLight.hh"
#include "WCSimTruePart.hh"
#include "WCSimTrueCapture.hh"


ClassImp(WCLTreeReader)



Bool_t WCLTreeReader::TouchData()
{

  return true;
}

void WCLTreeReader::LoadData(const char* file)
{
  std::cout << " *** WCLTreeReader::LoadData(...) *** " << std::endl;
  std::cout << "  adding: " << file << std::endl;

  fChain->Add(file);
  std::cout << "   ... total entries=" << fChain->GetEntries() << std::endl;
}


void WCLTreeReader::LoadData(const char* dfile, const char* gfile)
{
  std::cout << " *** WCLTreeReader::LoadData(...) *** " << std::endl;
  std::cout << "  adding: " << dfile << " and " <<gfile<<std::endl;

  fincgen=1;
  fChain->Add(dfile);
  fChainT->Add(gfile);

  std::cout << "   ... total entries=" << fChainT->GetEntries() << std::endl;
}

void WCLTreeReader::LoadEvent(Int_t ievent)
{
  //  return WCLTreeReader::Instance()->BuildEvent(ievent);
  fChain->GetEntry(ievent);
  if(fincgen==1) fChainT->GetEntry(ievent);
}


WCSimTrueLight* WCLTreeReader::getPhot(Int_t ip)
{

  WCSimTrueLight* thisphot;

  if(ip>nphot){
    thisphot = new WCSimTrueLight(0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0,0,0,0,0,0);
    std::cout<<"asking for more photons than were produced in this event"<<std::endl;
  } else{
    thisphot = new WCSimTrueLight(phot_xStart[ip],phot_yStart[ip],phot_zStart[ip],phot_tStart[ip],
				  phot_xEnd[ip],phot_yEnd[ip],phot_zEnd[ip],phot_tEnd[ip],
				  phot_wavelength[ip],phot_processStart[ip], phot_isScat[ip],
				  phot_parentid[ip],phot_trackid[ip],phot_hit[ip],phot_capnum[ip],0);
  }
  return thisphot;
}


WCSimTruePart* WCLTreeReader::getPart(Int_t ip){

  WCSimTruePart* thispart;

  if(ip>npart){
    thispart = new WCSimTruePart(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0,0,0,0);
    std::cout<<"asking for more particles than were produced in this event"<<std::endl;
  } else{
    thispart = new WCSimTruePart(part_xStart[ip],part_yStart[ip],part_zStart[ip],part_tStart[ip],
				 part_xEnd[ip],part_yEnd[ip],part_zEnd[ip],part_tEnd[ip],
				 part_pxStart[ip],part_pyStart[ip],part_pzStart[ip],
				 part_pxEnd[ip],part_pyEnd[ip],part_pzEnd[ip],part_KEstart[ip],part_KEend[ip],
				 part_processStart[ip],part_processEnd[ip],part_parentid[ip],part_trackid[ip],
				 part_pid[ip]);  
  }
  return thispart;
}


WCSimTrueCapture* WCLTreeReader::getCapt(Int_t ic){

  WCSimTrueCapture* thiscapt;

  if(ic>ncapturecount){
    thiscapt = new WCSimTrueCapture(0.,0.,0.,0.,0.,0,0,0,0,0);
    std::cout<<"asking for more particles than were produced in this event"<<std::endl;
  } else{
    thiscapt = new WCSimTrueCapture(capt_x[ic],capt_y[ic],capt_z[ic],capt_t0[ic],capt_E[ic],
				 capt_num[ic],capt_nucleus[ic],capt_pid[ic],capt_ngamma[ic],capt_nphot[ic]);

  }
  return thiscapt;
}



Int_t WCLTreeReader::GetEntries()
{
  return fChain->GetEntries();
}

Bool_t WCLTreeReader::TouchEvent()
{
  //  return WCLTreeReader::Instance()->CheckEvent();
  return true;
}

/*
void WCLTreeReader::Reset()
{
  return WCLTreeReader::Instance()->ResetForNewSample();
}
*/

WCLTreeReader::WCLTreeReader(int filetype, int includegen)
{

  TString tname;
  if(filetype==1) tname+="SmearedEventTree";
  else tname+="EventTree";
 
  // create event chain
  fChain = new TChain(tname,"chain");
  if(includegen==1) fChainT = new TChain("tcardfile","chain");
  //  this->ResetForNewSample();

  this->SetEventBranchAddys(filetype);
  if(includegen==1) this->SetGenBranchAddys();

  std::cout<<fChain->GetEntries()<<std::endl;
  
  fincgen=includegen;

}


WCLTreeReader::~WCLTreeReader()
{


}

void WCLTreeReader::ResetForNewSample()
{
  if( fChain->GetEntries()>0 ){
    std::cout << " *** WCLTreeReader::Reset() *** " << std::endl;  
  }

  // event chain
  /*
  if( fChain ){
    fChain->Reset();
    fChain->SetBranchAddress("wcsimrootevent",&fEvent);
    if( fEvent ) delete fEvent; fEvent = 0;
  }
  */
  //  fCounter = -1;

  return;
}



/*
void WCLTreeReader::BuildEvent(Int_t ievent)
{
  std::cout << " *** WCLTreeReader::BuildEvent(" << ievent << ") *** " << std::endl;


  return;
}

WCSimRootEvent* WCLTreeReader::GetWCSimEvent(Int_t ievent)
{
  this->BuildEvent(ievent);

  //  return fEvent;
}
*/


void WCLTreeReader::SetEventBranchAddys(int iS)
{
  bool isSmeared=false;
  if(iS==1) isSmeared = true;

  const int knphotmax=1000000;
  phot_xStart = new double[knphotmax];
  phot_yStart = new double[knphotmax];
  phot_zStart = new double[knphotmax];
  phot_tStart = new double[knphotmax];
  phot_xEnd = new double[knphotmax];
  phot_yEnd = new double[knphotmax];
  phot_zEnd = new double[knphotmax];
  phot_tEnd = new double[knphotmax];
  phot_wavelength = new double[knphotmax];

  phot_processStart = new int[knphotmax];
  phot_isScat = new int[knphotmax];
  phot_parentid = new int[knphotmax];
  phot_trackid = new int[knphotmax];
  phot_hit = new int[knphotmax];
  phot_capnum = new int[knphotmax]; 
 

  phot_PMTid = new int[knphotmax];


  /*
    phot_pxStart = new double[knphotmax];
    phot_pyStart = new double[knphotmax];
    phot_pyStart = new double[knphotmax];
    phot_pxEnd = new double[knphotmax];
    phot_pyEnd = new double[knphotmax];
    phot_pyEnd = new double[knphotmax];
  */
  
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

  fChain->SetBranchAddress("evt",&eventcount);
  fChain->SetBranchAddress("nphot",&nphot);
  fChain->SetBranchAddress("npart",&npart);
  fChain->SetBranchAddress("ncapturecount",&ncapturecount);
  fChain->SetBranchAddress("neutroncount",&neutroncount);
  fChain->SetBranchAddress("phot_xStart",phot_xStart);
  fChain->SetBranchAddress("phot_yStart",phot_yStart);
  fChain->SetBranchAddress("phot_zStart",phot_zStart);
  fChain->SetBranchAddress("phot_tStart",phot_tStart);
  fChain->SetBranchAddress("phot_xEnd",phot_xEnd);
  fChain->SetBranchAddress("phot_yEnd",phot_yEnd);
  fChain->SetBranchAddress("phot_zEnd",phot_zEnd);
  fChain->SetBranchAddress("phot_tEnd",phot_tEnd);
  fChain->SetBranchAddress("phot_wavelength",phot_wavelength);
  fChain->SetBranchAddress("phot_processStart",phot_processStart);
  fChain->SetBranchAddress("phot_isScat",phot_isScat);
  fChain->SetBranchAddress("phot_parentid",phot_parentid);
  fChain->SetBranchAddress("phot_trackid",phot_trackid);
  fChain->SetBranchAddress("phot_hit",phot_hit);
  fChain->SetBranchAddress("phot_capnum",phot_capnum);
  if(isSmeared)  fChain->SetBranchAddress("phot_PMTid",phot_PMTid);
  
  fChain->SetBranchAddress("part_xStart",part_xStart);
  fChain->SetBranchAddress("part_yStart",part_yStart);
  fChain->SetBranchAddress("part_zStart",part_zStart);
  fChain->SetBranchAddress("part_tStart",part_tStart);
  fChain->SetBranchAddress("part_xEnd",part_xEnd);
  fChain->SetBranchAddress("part_yEnd",part_yEnd);
  fChain->SetBranchAddress("part_zEnd",part_zEnd);
  fChain->SetBranchAddress("part_tEnd",part_tEnd);
  fChain->SetBranchAddress("part_pxStart",part_pxStart);
  fChain->SetBranchAddress("part_pyStart",part_pyStart);
  fChain->SetBranchAddress("part_pzStart",part_pzStart);
  fChain->SetBranchAddress("part_pxEnd",part_pxEnd);
  fChain->SetBranchAddress("part_pyEnd",part_pyEnd);
  fChain->SetBranchAddress("part_pzEnd",part_pzEnd);
  fChain->SetBranchAddress("part_KEstart",part_KEstart);
  fChain->SetBranchAddress("part_KEend",part_KEend);
  fChain->SetBranchAddress("part_processStart",part_processStart);
  fChain->SetBranchAddress("part_processEnd",part_processEnd);
  fChain->SetBranchAddress("part_parentid",part_parentid);
  fChain->SetBranchAddress("part_trackid",part_trackid);
  fChain->SetBranchAddress("part_pid",part_pid);
  
  fChain->SetBranchAddress("capt_x",capt_x);
  fChain->SetBranchAddress("capt_y",capt_y);
  fChain->SetBranchAddress("capt_z",capt_z);
  fChain->SetBranchAddress("capt_t0",capt_t0);
  fChain->SetBranchAddress("capt_E",capt_E);
  fChain->SetBranchAddress("capt_num",capt_num);
  fChain->SetBranchAddress("capt_pid",capt_num);
  fChain->SetBranchAddress("capt_nucleus",capt_nucleus);
  fChain->SetBranchAddress("capt_nphot",capt_nphot);
  fChain->SetBranchAddress("capt_ngamma",capt_ngamma);



  if(isSmeared){
    std::cout<<"We are indeed smeared"<<std::endl;
    const int kmaxtrk=100;
    mpid = new int[kmaxtrk];
    mpx = new double[kmaxtrk];
    mpy = new double[kmaxtrk];
    mpz = new double[kmaxtrk];
    mKE = new double[kmaxtrk];

    fChain->SetBranchAddress("mode",&mmode);
    fChain->SetBranchAddress("neutrino_E",&mE);
    fChain->SetBranchAddress("neutrino_id",&mbeam_id);
    fChain->SetBranchAddress("neutrino_px",&mbeam_px);
    fChain->SetBranchAddress("neutrino_py",&mbeam_py);
    fChain->SetBranchAddress("neutrino_pz",&mbeam_pz);
    fChain->SetBranchAddress("ntrks",&ntrks);
    fChain->SetBranchAddress("nneutrons",&nneutrons);
    fChain->SetBranchAddress("vtxx",&mvtxx); 
    fChain->SetBranchAddress("vtxy",&mvtxy); 
    fChain->SetBranchAddress("vtxz",&mvtxz);
    fChain->SetBranchAddress("mpid",mpid);
    fChain->SetBranchAddress("px",mpx);
    fChain->SetBranchAddress("py",mpy);
    fChain->SetBranchAddress("pz",mpz);
    fChain->SetBranchAddress("KE",mKE);
  }

  
}



void WCLTreeReader::SetGenBranchAddys(){

  std::cout<<"setting generator branch addresses"<<std::endl;

  const int kmaxtrk=100;
  const int kmaxMRDhits=50;
  mpid = new int[kmaxtrk];
  mpx = new double[kmaxtrk];
  mpy = new double[kmaxtrk];
  mpz = new double[kmaxtrk];
  mKE = new double[kmaxtrk];

  mMRDhitlayer = new int[kmaxMRDhits];
  mMRDhitorientation = new int[kmaxMRDhits];
  mMRDhitEdep = new double[kmaxMRDhits];
  mMRDhitEchdep = new double[kmaxMRDhits];
  mMRDhitx = new double[kmaxMRDhits];
  mMRDhity = new double[kmaxMRDhits];
  mMRDhitz = new double[kmaxMRDhits];

  
  fChainT->SetBranchAddress("evt",&mevt);
  fChainT->SetBranchAddress("mode",&mmode);
  fChainT->SetBranchAddress("neutrino_E",&mE);
  fChainT->SetBranchAddress("neutrino_id",&mbeam_id);
  fChainT->SetBranchAddress("neutrino_px",&mbeam_px);
  fChainT->SetBranchAddress("neutrino_py",&mbeam_py);
  fChainT->SetBranchAddress("neutrino_pz",&mbeam_pz);
  fChainT->SetBranchAddress("ntrks",&ntrks);
  fChainT->SetBranchAddress("nneutrons",&nneutrons);
  fChainT->SetBranchAddress("vtxx",&mvtxx); 
  fChainT->SetBranchAddress("vtxy",&mvtxy); 
  fChainT->SetBranchAddress("vtxz",&mvtxz);
  fChainT->SetBranchAddress("mpid",mpid);
  fChainT->SetBranchAddress("px",mpx);
  fChainT->SetBranchAddress("py",mpy);
  fChainT->SetBranchAddress("pz",mpz);
  fChainT->SetBranchAddress("KE",mKE);

  fChainT->SetBranchAddress("nMRDlayers",&nMRDlayers); 
  fChainT->SetBranchAddress("mMRDtotEdep",&mMRDtotEdep); 
  fChainT->SetBranchAddress("mMRDtotEchdep",&mMRDtotEchdep);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitlayer);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitorientation);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitEdep);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitEchdep);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitx);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhity);
  fChainT->SetBranchAddress("mMRDtotEchdep",mMRDhitz);
  
}
