#ifndef WCLTREEWRITER_HH
#define WCLTREEWRITER_HH

#include "TObject.h"
#include "TChain.h"
#include "TFile.h"

#include "WCSimTrueLight.hh"
#include "WCSimTruePart.hh"
#include "WCSimTrueCapture.hh"
#include "WCLTreeReader.hh"

#include <vector>
//#include "TRandom.h"
//#include "TRandom3.h"

class WCLTreeWriter : public TObject {

 public: 

  WCLTreeWriter(TString filename, int mode);
  ~WCLTreeWriter();

  void SetOutBranchAddys();
  void SetFlatBranchAddys();
  
  Int_t GetEntries();

  void InputEvent(WCLTreeReader* fEvent);
  void InitializeEvent();
  void AddWholeBranches(WCLTreeReader* fEvent, int addphot, int addpart, int addcapt, int addmrd, int addgen);
  void AddPhoton(WCSimTrueLight* fPhot);
//  void AddParticle(WCSimTruePart* fPart);
//  void AddCapture(WCSimTrueCapture* fCapt);
  void FillEvent();
  
  void WriteTreeToFile();

 
 private:

  void InputEventContents(WCLTreeReader* fEvent, int whole_event, int addphot, int addpart, int addcapt, int addmrd, int addgen);

  void ResetForNewSample();

  TTree* fChain;
  WCSimTrueLight* thisphot;
  int fmode;
  bool fNewEvent;

  //main file
  TFile* tf;
  //main event tree
  TTree* evttree;
  //main tree varialbles
  Double_t *phot_xStart,*phot_yStart,*phot_zStart,*phot_tStart,*phot_xEnd,*phot_yEnd,*phot_zEnd,*phot_tEnd,*phot_wavelength;
  Int_t *phot_processStart,*phot_isScat,*phot_parentid,*phot_trackid,*phot_hit,*phot_capnum;
  Double_t *part_xStart,*part_yStart,*part_zStart,*part_tStart,*part_xEnd,*part_yEnd,*part_zEnd,*part_tEnd,*part_KEstart,*part_KEend,*part_pxStart,*part_pyStart,*part_pzStart,*part_pxEnd,*part_pyEnd,*part_pzEnd;
  int *part_processStart,*part_processEnd,*part_parentid,*part_trackid,*part_pid;
  Int_t *capt_num,*capt_nucleus,*capt_pid,*capt_nphot,*capt_ngamma;
  Double_t *capt_x,*capt_y,*capt_z,*capt_t0,*capt_E;
  Int_t eventcount, nphot, npart, ncapturecount, neutroncount;
  
  //gen file
  TFile* tcf;
  //gen tree
  TTree* tcardtree;
  //gen tree variables
  int mevt,ntrks,nneutrons,mmode,mbeam_id,nMRDlayers;
  double mvtxx,mvtxy,mvtxz,mE,mbeam_px,mbeam_py,mbeam_pz,mMRDtotEdep,mMRDtotEchdep;
  int *mpid,*mMRDhitlayer,*mMRDhitorientation;
  double *mpx,*mpy,*mpz,*mKE,*mMRDhitEdep,*mMRDhitEchdep,*mMRDhitx,*mMRDhity,*mMRDhitz;

  //photon hit variables
  int nhits;
  int *phot_PMTid;

  //Christoph tree variables
  double x_hit,y_hit,z_hit,true_time,true_time_corrected,photon_wavelength;
  double process;
  double evtID;
  int pidd;
  int mprocessStart,mtrackID,mparentID,mprocessStep;
  double pos_x,pos_y,pos_z,dir_x,dir_y,dir_z,tt,Et,Ekin;
  
  //output file and tree
  TFile* outfile;
  TTree* outtree;
  TTree* truetree;
  
  //output text file
  fstream* textout;

  ClassDef(WCLTreeWriter,0)
    
};

#endif
