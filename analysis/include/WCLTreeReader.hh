#ifndef WCLTREEREADER_HH
#define WCLTREEREADER_HH

#include "TObject.h"
#include "TChain.h"

#include "WCSimTrueLight.hh"
#include "WCSimTruePart.hh"
#include "WCSimTrueCapture.hh"

#include <vector>
//#include "TRandom.h"
//#include "TRandom3.h"

class WCLTreeReader : public TObject {

 public: 

  Bool_t TouchData();
  Bool_t TouchEvent();

  void LoadData(const char* file);
  void LoadData(const char* dfile, const char* gfile);
  void LoadEvent(Int_t ievent);
  void SetEventBranchAddys(int iS);
  void SetGenBranchAddys();

  Int_t get_nphot() {return nphot;}
  Int_t get_nhits() {return nhits;}
  Int_t get_npart() {return npart;}
  
  Double_t* get_phot_xStart() {return phot_xStart;}
  Double_t* get_phot_yStart() {return phot_yStart;}
  Double_t* get_phot_zStart() {return phot_zStart;}
  Double_t* get_phot_tStart() {return phot_tStart;}
  Double_t* get_phot_xEnd() {return phot_xEnd;}
  Double_t* get_phot_yEnd() {return phot_yEnd;}
  Double_t* get_phot_zEnd() {return phot_zEnd;}
  Double_t* get_phot_tEnd() {return phot_tEnd;}
  Double_t* get_phot_wavelength() {return phot_wavelength;}
  Int_t* get_phot_processStart() {return phot_processStart;}
  Int_t* get_phot_parentid() {return phot_parentid;}
  Int_t* get_phot_trackid() {return phot_trackid;}
  Int_t* get_phot_isScat() {return phot_isScat;}
  Int_t* get_phot_hit() {return phot_hit;}
  Int_t* get_phot_capnum() {return phot_capnum;}
  Int_t* get_phot_PMTid() {return phot_PMTid;}

  Double_t* get_part_xStart() {return part_xStart;}
  Double_t* get_part_yStart() {return part_yStart;}
  Double_t* get_part_zStart() {return part_zStart;}
  Double_t* get_part_tStart() {return part_tStart;}
  Double_t* get_part_xEnd() {return part_xEnd;}
  Double_t* get_part_yEnd() {return part_yEnd;}
  Double_t* get_part_zEnd() {return part_zEnd;}
  Double_t* get_part_tEnd() {return part_tEnd;}
  Double_t* get_part_pxStart() {return part_pxStart;}
  Double_t* get_part_pyStart() {return part_pyStart;}
  Double_t* get_part_pzStart() {return part_pzStart;}
  Double_t* get_part_KEstart() {return part_KEstart;}
  Double_t* get_part_pxEnd() {return part_pxEnd;}
  Double_t* get_part_pyEnd() {return part_pyEnd;}
  Double_t* get_part_pzEnd() {return part_pzEnd;}
  Double_t* get_part_KEend() {return part_KEend;}
  Int_t* get_part_processStart() {return part_processStart;}
  Int_t* get_part_processEnd() {return part_processEnd;}
  Int_t* get_part_parentid() {return part_parentid;}
  Int_t* get_part_trackid() {return part_trackid;}
  Int_t* get_part_pid() {return part_pid;}

  // capture info

  Int_t get_neutroncount() { return neutroncount; }
  Int_t get_ncapturecount() { return ncapturecount; }
  Double_t* get_capt_x() { return capt_x; }
  Double_t* get_capt_y() { return capt_y; }
  Double_t* get_capt_z() { return capt_z; }
  Double_t* get_capt_t0() { return capt_t0; }
  Double_t* get_capt_E() { return capt_E; }
  Int_t* get_capt_num() { return capt_num; }
  Int_t* get_capt_nucleus() { return capt_nucleus; }
  Int_t* get_capt_pid() { return capt_pid; }
  Int_t* get_capt_nphot() { return capt_nphot; }
  Int_t* get_capt_ngamma()  { return capt_ngamma; }

  // generator info

  Int_t get_genntrks() {return ntrks;}
  Int_t get_genmode() {return mmode;}
  Int_t get_genbeam_id() {return mbeam_id;}
  Int_t get_gennneutrons() {return nneutrons;}
  Int_t get_genevt() {return mevt;}
  Double_t get_genvtxx() {return mvtxx;}
  Double_t get_genvtxy() {return mvtxy;}
  Double_t get_genvtxz() {return mvtxz;}
  Double_t get_genE() {return mE;}
  Double_t get_genbeam_px() {return mbeam_px;}
  Double_t get_genbeam_py() {return mbeam_py;}
  Double_t get_genbeam_pz() {return mbeam_pz;}
  Int_t* get_genpid() {return mpid;}
  Double_t* get_genpx() {return mpx;}
  Double_t* get_genpy() {return mpy;}
  Double_t* get_genpz() {return mpz;}
  Double_t* get_genKE() {return mKE;}
  
  // MRD info

  Int_t get_nMRDlayers() {return nMRDlayers;}
  Double_t get_MRDtotEdep() {return mMRDtotEdep;}
  Double_t get_MRDtotEcgdep() {return mMRDtotEchdep;}
  Int_t* get_MRDhitlayer() {return mMRDhitlayer;}
  Int_t* get_MRDhitorientation() {return mMRDhitorientation;}
  Double_t* get_MRDhitEdep() {return mMRDhitEdep;}
  Double_t* get_MRDhitEchdep() {return mMRDhitEchdep;}
  Double_t* get_MRDhitx() {return mMRDhitx;}
  Double_t* get_MRDhity() {return mMRDhity;}
  Double_t* get_MRDhitz() {return mMRDhitz;}

  WCSimTrueLight* getPhot(int ip);
  WCSimTruePart* getPart(int ip);
  WCSimTrueCapture* getCapt(int ic);

  Int_t GetEntries();

  WCLTreeReader(int filetype, int includegen);
  ~WCLTreeReader();

  

 private:

  void ResetForNewSample();

  TChain* fChain;
  TChain* fChainT;

  int fincgen;

  //main file
  TFile* tf;
  //main event tree
  TTree* evttree;
  //main tree varialbles
  double *phot_xStart,*phot_yStart,*phot_zStart,*phot_tStart,*phot_xEnd,*phot_yEnd,*phot_zEnd,*phot_tEnd,*phot_wavelength;
  int *phot_processStart,*phot_isScat,*phot_parentid,*phot_trackid,*phot_hit,*phot_capnum;
  double *part_xStart,*part_yStart,*part_zStart,*part_tStart,*part_xEnd,*part_yEnd,*part_zEnd,*part_tEnd,*part_KEstart,*part_KEend,*part_pxStart,*part_pyStart,*part_pzStart,*part_pxEnd,*part_pyEnd,*part_pzEnd;
  int *part_processStart,*part_processEnd,*part_parentid,*part_trackid,*part_pid;
  int *capt_num,*capt_nucleus,*capt_pid,*capt_nphot,*capt_ngamma;
  double *capt_x,*capt_y,*capt_z,*capt_t0,*capt_E;
  int eventcount, nphot, npart, ncapturecount, neutroncount;


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

  ClassDef(WCLTreeReader,0)
    
};

#endif
