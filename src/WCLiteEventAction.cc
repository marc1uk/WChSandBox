//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: WCLiteEventAction.cc,v 1.3 2006/06/29 17:54:20 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 

#include "WCLiteEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "WCLiteTrajectory.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "TFile.h"
#include "TString.h"
#include "WCLiteDetectorConstruction.hh"
#include <stdexcept>      // std::out_of_range
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "MRDSD.hh"
#include "mrdPMTSD.hh"
#include "G4HCofThisEvent.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ConvertProcessNameToCode(TString processName);
G4int ConvertParticleNameToCode(TString particleName);
 
WCLiteEventAction::WCLiteEventAction()
{

  eventcount=1;

  textout = new std::fstream("TextOut.txt",std::ios::out);

  no = new TFile("FullEvent.root","RECREATE");
  evttree = new TTree("EventTree","EventTree");

  
  nR = new TRandom3();
  
  phot_xStart = new G4double[knphotmax];
  phot_yStart = new G4double[knphotmax];
  phot_zStart = new G4double[knphotmax];
  phot_tStart = new G4double[knphotmax];
  phot_xEnd = new G4double[knphotmax];
  phot_yEnd = new G4double[knphotmax];
  phot_zEnd = new G4double[knphotmax];
  phot_tEnd = new G4double[knphotmax];
  phot_wavelength = new G4double[knphotmax];
  phot_processStart = new G4int[knphotmax];
  phot_isScat = new G4int[knphotmax];
  phot_parentid = new G4int[knphotmax];
  phot_trackid = new G4int[knphotmax];
  phot_hit = new G4int[knphotmax];
  phot_capnum = new G4int[knphotmax];
  phot_mrdhit = new G4int[knphotmax];											// new!
  phot_mrdprimary = new G4int[knphotmax]; 	
  phot_mrdnumsecs = new G4int[knphotmax];	
  phot_mrdpritrackid = new G4int[knphotmax];
  phot_mrdedep = new G4double[knphotmax];
  phot_mrdstartx = new G4double[knphotmax];	
  phot_mrdstarty = new G4double[knphotmax];
  phot_mrdstartz = new G4double[knphotmax];
  phot_mrddetected = new G4int[knphotmax];
  
  /*
  phot_pxStart = new G4double[knphotmax];
  phot_pyStart = new G4double[knphotmax];
  phot_pyStart = new G4double[knphotmax];
  phot_pxEnd = new G4double[knphotmax];
  phot_pyEnd = new G4double[knphotmax];
  phot_pyEnd = new G4double[knphotmax];
  */
  
  part_xStart = new G4double[knpartmax];
  part_yStart = new G4double[knpartmax];
  part_zStart = new G4double[knpartmax];
  part_tStart = new G4double[knpartmax];
  part_xEnd = new G4double[knpartmax];
  part_yEnd = new G4double[knpartmax];
  part_zEnd = new G4double[knpartmax];
  part_tEnd = new G4double[knpartmax];
  part_KEstart = new G4double[knpartmax];
  part_KEend = new G4double[knpartmax];
  part_pxStart = new G4double[knpartmax];
  part_pyStart = new G4double[knpartmax];
  part_pzStart = new G4double[knpartmax];
  part_pxEnd = new G4double[knpartmax];
  part_pyEnd = new G4double[knpartmax];
  part_pzEnd = new G4double[knpartmax];
  part_processStart = new G4int[knpartmax];
  part_processEnd = new G4int[knpartmax];
  part_parentid = new G4int[knpartmax];
  part_trackid = new G4int[knpartmax];
  part_pid = new G4int[knpartmax];
  part_mrdhit = new G4int[knpartmax];	
  part_mrdprimary = new G4int[knpartmax];	
  part_mrdnumsecs = new G4int[knpartmax];	
  part_mrdpritrackid = new G4int[knpartmax];
  part_mrdedep = new G4double[knpartmax];	
  part_mrdstartx = new G4double[knpartmax];	
  part_mrdstarty = new G4double[knpartmax];
  part_mrdstartz = new G4double[knpartmax];
  part_mrddetected = new G4int[knpartmax];

  capt_num = new G4int[kcapmax];
  capt_nucleus = new G4int[kcapmax];
  capt_pid = new G4int[kcapmax];
  capt_nphot = new G4int[kcapmax];
  capt_ngamma = new G4int[kcapmax];
  capt_x = new G4double[kcapmax];
  capt_y = new G4double[kcapmax];
  capt_z = new G4double[kcapmax];
  capt_t0 = new G4double[kcapmax];
  capt_E = new G4double[kcapmax];

  evttree->Branch("evt",&eventcount);
  evttree->Branch("nphot",&nphot);
  evttree->Branch("npart",&npart);
  evttree->Branch("ncapturecount",&ncapturecount);
  evttree->Branch("neutroncount",&neutroncount);
  
  evttree->Branch("phot_xStart",phot_xStart,"phot_xStart[nphot]/D");
  evttree->Branch("phot_yStart",phot_yStart,"phot_yStart[nphot]/D");
  evttree->Branch("phot_zStart",phot_zStart,"phot_zStart[nphot]/D");
  evttree->Branch("phot_tStart",phot_tStart,"phot_tStart[nphot]/D");
  evttree->Branch("phot_xEnd",phot_xEnd,"phot_xEnd[nphot]/D");
  evttree->Branch("phot_yEnd",phot_yEnd,"phot_yEnd[nphot]/D");
  evttree->Branch("phot_zEnd",phot_zEnd,"phot_zEnd[nphot]/D");
  evttree->Branch("phot_tEnd",phot_tEnd,"phot_tEnd[nphot]/D");
  evttree->Branch("phot_wavelength",phot_wavelength,"phot_wavelength[nphot]/D");
  evttree->Branch("phot_processStart",phot_processStart,"phot_processStart[nphot]/I");
  evttree->Branch("phot_isScat",phot_isScat,"phot_isScat[nphot]/I");
  evttree->Branch("phot_parentid",phot_parentid,"phot_parentid[nphot]/I");
  evttree->Branch("phot_trackid",phot_trackid,"phot_trackid[nphot]/I");
  evttree->Branch("phot_hit",phot_hit,"phot_hit[nphot]/I");
  evttree->Branch("phot_capnum",phot_capnum,"phot_capnum[nphot]/I");
  evttree->Branch("phot_mrdhit",phot_mrdhit,"phot_mrdhit[nphot]/I");
  evttree->Branch("phot_mrdprimary",phot_mrdprimary,"phot_mrdprimary[nphot]/I");
  evttree->Branch("phot_mrdnumsecs",phot_mrdnumsecs,"phot_mrdnumsecs[nphot]/I");
  evttree->Branch("phot_mrdpritrackid",phot_mrdpritrackid,"phot_mrdpritrackid[nphot]/I");
  evttree->Branch("phot_mrdedep",phot_mrdedep,"phot_mrdedep[nphot]/D");
  evttree->Branch("phot_mrdstartx",phot_mrdstartx,"phot_mrdstartx[nphot]/D");
  evttree->Branch("phot_mrdstarty",phot_mrdstarty,"phot_mrdstarty[nphot]/D");
  evttree->Branch("phot_mrdstartz",phot_mrdstartz,"phot_mrdstartz[nphot]/D");
  evttree->Branch("phot_mrddetected",phot_mrddetected,"phot_mrddetected[nphot]/I");

  evttree->Branch("part_xStart",part_xStart,"part_xStart[npart]/D");
  evttree->Branch("part_yStart",part_yStart,"part_yStart[npart]/D");
  evttree->Branch("part_zStart",part_zStart,"part_zStart[npart]/D");
  evttree->Branch("part_tStart",part_tStart,"part_tStart[npart]/D");
  evttree->Branch("part_xEnd",part_xEnd,"part_xEnd[npart]/D");
  evttree->Branch("part_yEnd",part_yEnd,"part_yEnd[npart]/D");
  evttree->Branch("part_zEnd",part_zEnd,"part_zEnd[npart]/D");
  evttree->Branch("part_tEnd",part_tEnd,"part_tEnd[npart]/D");
  evttree->Branch("part_pxStart",part_pxStart,"part_pxStart[npart]/D");
  evttree->Branch("part_pyStart",part_pyStart,"part_pyStart[npart]/D");
  evttree->Branch("part_pzStart",part_pzStart,"part_pzStart[npart]/D");
  evttree->Branch("part_pxEnd",part_pxEnd,"part_pxEnd[npart]/D");
  evttree->Branch("part_pyEnd",part_pyEnd,"part_pyEnd[npart]/D");
  evttree->Branch("part_pzEnd",part_pzEnd,"part_pzEnd[npart]/D");
  evttree->Branch("part_KEstart",part_KEstart,"part_KEstart[npart]/D");
  evttree->Branch("part_KEend",part_KEend,"part_KEend[npart]/D");
  evttree->Branch("part_processStart",part_processStart,"part_processStart[npart]/I");
  evttree->Branch("part_processEnd",part_processEnd,"part_processEnd[npart]/I");
  evttree->Branch("part_parentid",part_parentid,"part_parentid[npart]/I");
  evttree->Branch("part_trackid",part_trackid,"part_trackid[npart]/I");
  evttree->Branch("part_pid",part_pid,"part_pid[npart]/I");
  evttree->Branch("part_mrdhit",part_mrdhit,"part_mrdhit[npart]/I");
  evttree->Branch("part_mrdprimary",part_mrdprimary,"part_mrdprimary[npart]/I");
  evttree->Branch("part_mrdnumsecs",part_mrdnumsecs,"part_mrdnumsecs[npart]/I");
  evttree->Branch("part_mrdpritrackid",part_mrdpritrackid,"part_mrdpritrackid[npart]/I");
  evttree->Branch("part_mrdedep",part_mrdedep,"part_mrdedep[npart]/D");
  evttree->Branch("part_mrdstartx",part_mrdstartx,"part_mrdstartx[npart]/D");
  evttree->Branch("part_mrdstarty",part_mrdstarty,"part_mrdstarty[npart]/D");
  evttree->Branch("part_mrdstartz",part_mrdstartz,"part_mrdstartz[npart]/D");
  evttree->Branch("part_mrddetected",part_mrddetected,"part_mrddetected[npart]/I");
  

  evttree->Branch("capt_x",capt_x,"capt_x[ncapturecount]/D");
  evttree->Branch("capt_y",capt_y,"capt_y[ncapturecount]/D");
  evttree->Branch("capt_z",capt_z,"capt_z[ncapturecount]/D");
  evttree->Branch("capt_t0",capt_t0,"capt_t0[ncapturecount]/D");
  evttree->Branch("capt_E",capt_E,"capt_E[ncapturecount]/D");
  evttree->Branch("capt_num",capt_num,"capt_num[ncapturecount]/I");
  evttree->Branch("capt_pid",capt_num,"capt_pid[ncapturecount]/I");
  evttree->Branch("capt_nucleus",capt_nucleus,"capt_nucleus[ncapturecount]/I");
  evttree->Branch("capt_nphot",capt_nphot,"capt_nphot[ncapturecount]/I");
  evttree->Branch("capt_ngamma",capt_ngamma,"capt_ngamma[ncapturecount]/I");
  
  // ***********************
  // MRD INFORMATION STORAGE
  // ***********************
  // 
  // INFO FROM HITS.   
  
  mrdfile = new TFile("MRDEvents.root","RECREATE");
  mrdtree = new TTree("MRDTree","MRDTree"); 

  mrdhit_x = new G4double[kmrdhitnmax];		// where/when was the hit?
  mrdhit_y = new G4double[kmrdhitnmax];
  mrdhit_z = new G4double[kmrdhitnmax];
  mrdhit_t = new G4double[kmrdhitnmax];	
  mrdhit_process = new G4int[kmrdhitnmax];	// what was the interaction process?
  mrdhit_particleID = new G4int[kmrdhitnmax];	// what was the particle type interacting?
  mrdhit_trackID = new G4int[kmrdhitnmax];	// what was the track ID
  mrdhit_edep = new G4double[kmrdhitnmax]; 	// how much energy was deposited?
  mrdhit_objnum = new G4int[kmrdhitnmax];	// which geometry object was hit?
  mrdhit_copynum = new G4int[kmrdhitnmax];	// which copy was hit?
  
  mrdtree->Branch("evt",&eventcount);
  mrdtree->Branch("hitnum",&hitnum);
  
  mrdtree->Branch("mrdhit_x",mrdhit_x,"mrdhit_x[hitnum]/D");
  mrdtree->Branch("mrdhit_y",mrdhit_y,"mrdhit_y[hitnum]/D");
  mrdtree->Branch("mrdhit_z",mrdhit_z,"mrdhit_z[hitnum]/D");
  mrdtree->Branch("mrdhit_t",mrdhit_t,"mrdhit_t[hitnum]/D");
  mrdtree->Branch("mrdhit_process",mrdhit_process,"mrdhit_process[hitnum]/I");
  mrdtree->Branch("mrdhit_particleID",mrdhit_particleID,"mrdhit_particleID[hitnum]/I");
  mrdtree->Branch("mrdhit_trackID",mrdhit_trackID,"mrdhit_trackID[hitnum]/I");
  mrdtree->Branch("mrdhit_edep",mrdhit_edep,"mrdhit_edep[hitnum]/D");
  mrdtree->Branch("mrdhit_objnum",mrdhit_objnum,"mrdhit_objnum[hitnum]/I");
  mrdtree->Branch("mrdhit_copynum",mrdhit_copynum,"mrdhit_copynum[hitnum]/I");
   
  // **********************************
  // MRD LG/PMT HIT INFORMATION STORAGE
  // **********************************
  // 
  // INFO FROM HITS. 
  
  pmttree = new TTree("PMTTree","PMTTree");
 
  pmthit_x = new G4double[kpmthitnmax];			// where/when was the hit?
  pmthit_y = new G4double[kpmthitnmax];
  pmthit_z = new G4double[kpmthitnmax];
  pmthit_t = new G4double[kpmthitnmax];	
  pmthit_process = new G4int[kpmthitnmax];		// what was the creation process?
  pmthit_trackID = new G4int[kpmthitnmax];		// what was the track ID
  pmthit_parentID = new G4int[kpmthitnmax];		// track ID of parent particle
  pmthit_wavelength = new G4double[kpmthitnmax]; 	// what was the hit wavelength?  
  pmthit_copynum = new G4int[kpmthitnmax];		// which PMT/LG was hit?
  
  pmttree->Branch("evt",&eventcount);
  pmttree->Branch("pmt_hitnum",&pmt_hitnum);
  
  pmttree->Branch("pmthit_x",pmthit_x,"pmthit_x[pmt_hitnum]/D");
  pmttree->Branch("pmthit_y",pmthit_y,"pmthit_y[pmt_hitnum]/D");
  pmttree->Branch("pmthit_z",pmthit_z,"pmthit_z[pmt_hitnum]/D");
  pmttree->Branch("pmthit_t",pmthit_t,"pmthit_t[pmt_hitnum]/D");
  pmttree->Branch("pmthit_process",pmthit_process,"pmthit_process[pmt_hitnum]/I");
  pmttree->Branch("pmthit_trackID",pmthit_trackID,"pmthit_trackID[pmt_hitnum]/I");
  pmttree->Branch("pmthit_parentID",pmthit_parentID,"pmthit_parentID[pmt_hitnum]/I");
  pmttree->Branch("pmthit_wavelength",pmthit_wavelength,"pmthit_wavelength[pmt_hitnum]/D");
  pmttree->Branch("pmthit_copynum",pmthit_copynum,"pmthit_copynum[pmt_hitnum]/I");
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
WCLiteEventAction::~WCLiteEventAction()
{
   no->cd();
   evttree->Write();
   no->Close();
   //delete no;
   
   mrdfile->cd();
   mrdtree->Write();
   pmttree->Write();
   mrdfile->Close();
   delete mrdfile;
   
   textout->close();
   //delete textout;

   delete[] phot_xStart;
   delete[] phot_yStart;
   delete[] phot_zStart;
   delete[] phot_tStart;
   delete[] phot_xEnd;
   delete[] phot_yEnd;
   delete[] phot_zEnd;
   delete[] phot_tEnd;
   delete[] phot_wavelength;
   delete[] phot_processStart;
   delete[] phot_isScat;
   delete[] phot_parentid;
   delete[] phot_trackid;
   delete[] phot_hit;
   delete[] phot_capnum;
   delete[] phot_mrdhit;
   delete[] phot_mrdprimary;
   delete[] phot_mrdnumsecs;
   delete[] phot_mrdpritrackid;
   delete[] phot_mrdedep;
   delete[] phot_mrdstartx;
   delete[] phot_mrdstarty;
   delete[] phot_mrdstartz;
   delete[] phot_mrddetected;
   delete[] part_xStart;
   delete[] part_yStart;
   delete[] part_zStart;
   delete[] part_tStart;
   delete[] part_xEnd;
   delete[] part_yEnd;
   delete[] part_zEnd;
   delete[] part_tEnd;
   delete[] part_KEstart;
   delete[] part_KEend;
   delete[] part_pxStart;
   delete[] part_pyStart;
   delete[] part_pzStart;
   delete[] part_pxEnd;
   delete[] part_pyEnd;
   delete[] part_pzEnd;
   delete[] part_processStart;
   delete[] part_processEnd;
   delete[] part_parentid;
   delete[] part_trackid;
   delete[] part_pid;
   delete[] part_mrdhit;
   delete[] part_mrdprimary;
   delete[] part_mrdnumsecs;
   delete[] part_mrdpritrackid;
   delete[] part_mrdedep;
   delete[] part_mrdstartx;
   delete[] part_mrdstarty;
   delete[] part_mrdstartz;
   delete[] part_mrddetected;

   delete[] capt_num;
   delete[] capt_nucleus;
   delete[] capt_pid;
   delete[] capt_nphot;
   delete[] capt_ngamma;
   delete[] capt_x;
   delete[] capt_y;
   delete[] capt_z;
   delete[] capt_t0;
   delete[] capt_E;
   
   delete[] mrdhit_x;
   delete[] mrdhit_y;
   delete[] mrdhit_z;
   delete[] mrdhit_t;
   delete[] mrdhit_process;
   delete[] mrdhit_particleID;
   delete[] mrdhit_trackID;
   delete[] mrdhit_edep;
   delete[] mrdhit_objnum;
   delete[] mrdhit_copynum;
   

   delete[] pmthit_x;
   delete[] pmthit_y;
   delete[] pmthit_z;
   delete[] pmthit_t;
   delete[] pmthit_process;
   delete[] pmthit_trackID;
   delete[] pmthit_parentID;
   delete[] pmthit_wavelength;
   delete[] pmthit_copynum;


   G4cout<<"End of this run. Number of events: "<< eventcount <<" no evts with muons: "<< totalmucount<< " frac of evts with muons: "<<(totalmucount/eventcount)<<G4endl;
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void WCLiteEventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void WCLiteEventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4int myneutroncount=0;
  //  G4cout <<"**********************"<<G4endl;
  G4cout <<" "<<G4endl <<" "<<G4endl; 
  G4cout <<"START EndOfEventAction"<<G4endl;

  std::vector<double> ngammaspercaptcha;
  std::vector<double> nfinalphotspercaptcha;
  std::vector<double> energypercaptcha;
  std::vector<double> t0ofcaptcha,xofcaptcha,yofcaptcha,zofcaptcha,nucofcaptcha;
  std::vector<double> t0ofcaptcha_nuc,xofcaptcha_nuc,yofcaptcha_nuc,zofcaptcha_nuc,pidofcaptcha_nuc;
  std::vector<double> neutron_xStart,neutron_yStart,neutron_zStart,neutron_tStart,neutron_xEnd,neutron_yEnd,neutron_zEnd,neutron_tEnd;
  std::vector<int> neutron_processStart,neutron_processEnd;
  std::vector<int> parentidofcaptcha, parentidofcaptcha_nuc;
  std::vector<std::vector<double> > photons;
  std::vector<double> vpidtrack;
  std::vector<double> vtrackid;
  
  ncapturecount=0;
  neutroncount=0;   
  nphot=0;
  npart=0;	
  G4int n_trajectories=0, n_earlyEMparts=0, n_earlyMuons=0, photcount=0, mpcount=0, neithercount=0, numloops=0;
    
  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
  if (trajectoryContainer==0){G4cout << "trajectoryContainer null!" <<G4endl;}
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  /*
  for (G4int i=0; i<n_trajectories; i++){ 
    vpidtrack.push_back(0);
  }
  */
  G4double eventstarttime=99999999999;
  G4double eventduration=0;
  for (G4int i=0; i<n_trajectories; i++){ 

    TEStart = -5555; // Initialisation
    //G4cout << "analysing trajectory " << i << " of " << n_trajectories << G4endl;
    //G4Trajectory* trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
    WCLiteTrajectory* trj=0;
    try{
    	trj = (WCLiteTrajectory*)(anEvent->GetTrajectoryContainer()->GetVector()->at(i));
    } 
    catch (const std::out_of_range& oor) {
    	std::cerr << "Out of Range error: " << oor.what() << '\n';
    }
    if (trj==0){G4cout << "trj null!" <<G4endl; break;}
    //G4cout << "retrieved trajectory: " << trj << G4endl;
        
    trackid = trj->GetTrackID();
    //G4cout << "parent trj " << trackid << G4endl;
    parentid = trj->GetParentID();
    
    TString partname = (trj->GetParticleName());
    partcode = ConvertParticleNameToCode(partname);
    
    G4int numPoints = trj->GetPointEntries();
    //G4cout << "trj has " << numPoints << " points." << G4endl;
    WCLiteTrajectoryPoint* mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(numPoints-1));
    //G4cout << "got point " << (numPoints-1) << G4endl;
    if (mypnt==0){G4cout << "mypnt-1 null!" <<G4endl;}
    //G4cout << "mypnt = " << mypnt << G4endl;
    //G4cout << "getting position " << (mypnt->GetPosition()) << G4endl;
    
    /*if(eventcount==3 && i==68557){
	    G4cout << "particle: " << partname << " process " << (TString)(mypnt->GetProcessName());
	    G4cout << " track ID: " << trj->GetTrackID() << " passesthroughMRD: " << trj->GetPassThruMRD();
	    G4cout << " sameasparenttrackid: " << trj->GetSameAsParentTrackID() << " mrdOriginalTrackID: ";
	    G4cout << trj->GetmrdOriginalTrackID() << " total E: " <<  mypnt->GetTotalEnergy();
	    WCLiteTrajectoryPoint* mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(0));
	    G4cout << " zero point position: " << (mypnt->GetPosition()) << G4endl;
	    G4cout << "trajectory retrieved: " << trj << " has " << numPoints << " points. " << G4endl;
	    mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(numPoints-1));
	    G4cout << "point " << (numPoints-1) << " is: " << (mypnt->GetPosition()) << G4endl;
    }*/
  
    xEnd = mypnt->GetPosition().x();
    yEnd = mypnt->GetPosition().y();
    zEnd = mypnt->GetPosition().z();
    xEndDir = mypnt->GetDirection().x();
    yEndDir = mypnt->GetDirection().y();
    zEndDir = mypnt->GetDirection().x();
    tEnd = mypnt->GetGlobalTime();
    TEEnd = mypnt->GetTotalEnergy();
    KEEnd = mypnt->GetKineticEnergy();  
    //G4cout << "got end stats" <<G4endl;
    if(tEnd>eventduration){eventduration=tEnd;}
     
    TString theProcess = mypnt->GetProcessName();
    processEnd = ConvertProcessNameToCode(theProcess);
    //if(processEnd==0){G4cout<<"Particle left via transportation at: " << zEnd<<G4endl;}
    
    //G4cout << "end process is: " << theProcess << G4endl;
     phhit=0;
     if(partcode==100 && mypnt->GetDetectorLocation()=="DetectorHit"){
       phhit=1;
     }

     //   tEnd=5.;

     mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(0));
     if (mypnt==0){G4cout << "mypnt0 null!" <<G4endl;}
     xStart = mypnt->GetPosition().x();
     yStart = mypnt->GetPosition().y();
     zStart = mypnt->GetPosition().z();
     //if(zStart<-2000){G4cout<<"particle with zStart at: " <<zStart<<" ceated with start process " << theProcess<<G4endl;}
     xStartDir = mypnt->GetDirection().x();
     yStartDir = mypnt->GetDirection().y();
     zStartDir = mypnt->GetDirection().z();
     tStart = mypnt->GetGlobalTime();
     TEStart = mypnt->GetTotalEnergy();
     KEStart = mypnt->GetKineticEnergy();
    //G4cout << "KE start is: " << KEStart << G4endl;
    if(tStart<eventstarttime){eventstarttime=tStart;}

     //G4cout << "checking mrd hit stats" <<G4endl;
     mrdhit = trj->GetPassThruMRD();	// actually says does it pass thru MRD, doesn't say if it 'hits' the MRD?
     if(mrdhit!=0){
	     ismrdprimary =trj->GetIsMRDprimary();
	     mrdnumsecs= trj->GetNumSecs();
	     mrdpritrackid=trj->GetmrdOriginalTrackID();
	     mrdedep= 0;
	     G4ThreeVector mrdStartPos=trj->GetMRDstartPos();
	     mrdstartx=mrdStartPos.x();
	     mrdstarty=mrdStartPos.y();
	     mrdstartz=mrdStartPos.z();
     }
     mrddetected=trj->GetMRDdetected();
     //G4cout << "mrd hit stats ok" << G4endl;

//     mrdhit=0;	// actually says does it pass thru MRD, doesn't say if it 'hits' the MRD?
//     ismrdprimary=0;
//     mrdnumsecs=0;
//     mrdedep=0;
//     mrdpritrackid=0;
	
     //     if( (partcode==100) && (tStart<100) ) G4cout<<"hup "<<tStart<<G4endl
     /*;
     if(partcode==2112&&processEnd==11){
       G4cout<<"NEUTRONSTOPPED AND DID CAPTURE!!!"<<G4endl;
       G4cout<<"npoints: "<<numPoints<<G4endl;
       G4cout<<KEStart<<" "<<(tEnd-tStart)<<" "<<(zEnd-zStart)<<G4endl;
     }

     if(partcode==2112&&processEnd!=11){
       G4cout<<"NEUTRONSTOPPED BUT DIDN'T CAPTURE!!!"<<G4endl;
       G4cout<<"npoints: "<<numPoints<<G4endl;
       G4cout<<KEStart<<" "<<(tEnd-tStart)<<" "<<(zEnd-zStart)<<G4endl;
     }
     */

     TString sProcess = mypnt->GetProcessName();
     processStart=ConvertProcessNameToCode(sProcess);
     
     if(processStart==-1){
     //  std::cout << "NEW WARNING: PROCESS START **NOT** DEFINED!!! (" << sProcess << ")" << std::endl;
     }

     if(KEStart>0) {wavelengthStart = (4.13566733*2.99792458 *0.0001 )/KEStart;}

          
     if(processStart==11){
     //std::cout<<"VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"<<std::endl;
     //  std::cout<<"Particle Name = "<<partname<<std::endl;
	 for(int qq=0; qq<numPoints; qq++){
	   WCLiteTrajectoryPoint* tpnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(qq));
	//   std::cout<<qq<<" "<<tpnt->GetProcessName()<<std::endl;
	 }
     //std::cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl;
     }

     //G4cout << "checking partcode 100" << G4endl;
     if(partcode==100){

       double QE50prescale = nR->Rndm();
       if(QE50prescale<0.5){


       isScatPhot=0;
       if( ((xStartDir-xEndDir)!=0) || ((yStartDir-yEndDir)!=0) || ((yStartDir-yEndDir)!=0) ){
	 //G4cout<<"scattered? "<<(xStartDir-xEndDir)<<" "<<(yStartDir-yEndDir)<<" "<<(yStartDir-yEndDir)<<G4endl;
	 isScatPhot=1;  
       }

       if(nphot>knphotmax) G4cout<<"max number of photons in phot array (1000000) has been exeeded"<<G4endl;

       phot_xStart[nphot]=xStart;
       phot_yStart[nphot]=yStart;
       phot_zStart[nphot]=zStart;
       phot_tStart[nphot]=tStart;
       phot_xEnd[nphot]=xEnd;
       phot_yEnd[nphot]=yEnd;
       phot_zEnd[nphot]=zEnd;
       //if(zEnd>2000){G4cout<<"photon with zEnd > 2000 has zStart "<<zStart<<G4endl;}
       phot_tEnd[nphot]=tEnd;
       phot_wavelength[nphot]=wavelengthStart;
       phot_processStart[nphot]= (G4int)processStart;
       phot_isScat[nphot]=  (G4int)isScatPhot;
       phot_parentid[nphot]= (G4int)parentid;
       phot_trackid[nphot]= (G4int)trackid;
       phot_hit[nphot]= (G4int)phhit;
       phot_mrdhit[nphot]=(G4int)mrdhit;
       phot_mrdprimary[nphot] = (G4int)ismrdprimary;
       phot_mrdnumsecs[nphot] = (G4int)mrdnumsecs;
       phot_mrdpritrackid[nphot]=(G4int)mrdpritrackid;
       phot_mrdedep[nphot] = (G4double)mrdedep;
       phot_mrdstartx[nphot] = (G4double)mrdstartx;
       phot_mrdstarty[nphot] = (G4double)mrdstarty;
       phot_mrdstartz[nphot] = (G4double)mrdstartz;
       phot_mrddetected[nphot] = (G4int)mrddetected;
       //mrdphots.push_back(nphot);

       nphot++;	// photon counter: initialised before the loop over trajectories and incremented here, to track position in photon arrays.
       //       if(nphot%100==0) G4cout<<nphot<<" PPPP "<<G4endl;
       }
     }else{

       if(npart>knpartmax) G4cout<<"max number of particles in part array (10000) has been exeeded"<<G4endl;

       part_xStart[npart]=xStart;
       part_yStart[npart]=yStart;
       part_zStart[npart]=zStart;
       part_tStart[npart]=tStart;
       part_xEnd[npart]=xEnd;
       part_yEnd[npart]=yEnd;
       part_zEnd[npart]=zEnd;
       part_tEnd[npart]=tEnd;
       part_pxStart[npart]=xStartDir;
       part_pyStart[npart]=yStartDir;
       part_pzStart[npart]=zStartDir;
       part_pxEnd[npart]=xEndDir;
       part_pyEnd[npart]=yEndDir;
       part_pzEnd[npart]=zEndDir;
       part_KEstart[npart]=KEStart;
       part_KEend[npart]=KEEnd;
       part_processStart[npart]= (G4int)processStart;
       part_processEnd[npart]= (G4int)processEnd;
       part_parentid[npart]= (G4int)parentid;
       part_trackid[npart]= (G4int)trackid;
       part_pid[npart]= (G4int)partcode;
       part_mrdhit[npart]=(G4int)mrdhit;
       part_mrdprimary[npart] = (G4int)ismrdprimary;
       part_mrdnumsecs[npart] = (G4int)mrdnumsecs;
       part_mrdpritrackid[npart]=(G4int)mrdpritrackid;
       part_mrdedep[npart] = (G4double)mrdedep;
       part_mrdstartx[npart] = (G4double)mrdstartx;
       part_mrdstarty[npart] = (G4double)mrdstarty;
       part_mrddetected[npart] = (G4int)mrddetected;
       //mrdparts.push_back(npart);

       //if(npart%5==0) G4cout<<"npart "<<npart<<G4endl;
       
       npart++;
     }

     // Stats on neutrons and neutron capture
     // if(partcode==2112) neutroncount++;
     //G4cout << "checking partcode 2112" << G4endl;
     if(partcode==2112){	// push all neutron start and endpoints into a vector (used for secondary neutron counting later)
				// but completely separate from nCapture counting below. 
       neutron_xStart.push_back(part_xStart[npart]);
       neutron_yStart.push_back(part_yStart[npart]);
       neutron_zStart.push_back(part_zStart[npart]);
       neutron_tStart.push_back(part_tStart[npart]);
       neutron_xEnd.push_back(part_xEnd[npart]);
       neutron_yEnd.push_back(part_yEnd[npart]);
       neutron_zEnd.push_back(part_zEnd[npart]);
       neutron_tEnd.push_back(part_tEnd[npart]);
       neutron_processStart.push_back(part_processStart[npart]);
       neutron_processEnd.push_back(part_processEnd[npart]);
       if((trj->GetSameAsParentTrackID())==0){myneutroncount++;}
     }

     if( (processStart==11) ){	// in summary: for all particles from an nCapture event, store start times and parent IDs
				// unless they already exist. Otherwise, increment energy & gamma count stats for that capture. 
				// num captures is found from capture events that haven't already been identified. 
				// BUT this is modified later.
       //G4cout<<"a captcha"<<G4endl;
       //std::cout << "Neutron capture! Paticle code = " << partcode << std::endl;
	
       if(partcode!=22){				// if particle from nCapture is NOT a gamma
	 t0ofcaptcha_nuc.push_back(tStart);		// store emission point of non-gamma daughter
	 xofcaptcha_nuc.push_back(xStart);
	 yofcaptcha_nuc.push_back(yStart);
	 zofcaptcha_nuc.push_back(zStart);
	 pidofcaptcha_nuc.push_back(partcode);
	 parentidofcaptcha_nuc.push_back(parentid);	// and parent neutron of that track
       }

       if(partcode==22){							// for gammas

	 std::cout<<"parent id: "<<parentid<<std::endl;
	 
	 bool isGd=false;

	 int whichcaptcha=-1;
	 for(int m=0; m<t0ofcaptcha.size(); m++){				// scan non-gamma capture events 

	   if(((int)parentid)==parentidofcaptcha.at(m)) whichcaptcha=m;		// for the same parent neutron.

	   //   if( fabs(xStart-xofcaptcha.at(m))<1.0 && fabs(yStart-yofcaptcha.at(m))<1.0 && fabs(zStart-zofcaptcha.at(m))<1.0 ) whichcaptcha=m;
	   //   if( fabs(tStart-t0ofcaptcha.at(m)) <5.0) whichcaptcha=m;
	 }
	 
	 if(whichcaptcha==-1){							// if not found

	   ncapturecount++;							// increment capture count (gamma-only capture)
	   t0ofcaptcha.push_back(tStart);					// add capture stats to the (same) list
	   xofcaptcha.push_back(xStart);
	   yofcaptcha.push_back(yStart);
	   zofcaptcha.push_back(zStart);
	   energypercaptcha.push_back(TEStart);
	   ngammaspercaptcha.push_back(1.0);					// say it had one gamma from this capture
	   parentidofcaptcha.push_back((int)parentid);				// and it's parent neutron track
	   
	   std::cout <<"Capture energy = " << TEStart << std::endl;
	   std::cout << std::endl;

	 } else{								// if found
	   double theE = energypercaptcha.at(whichcaptcha);			// increment energy from this capture
	   double theNgam = ngammaspercaptcha.at(whichcaptcha);
	   
	   theE+=TEStart;
	   theNgam+=1.0;
	   
	   energypercaptcha.at(whichcaptcha) = theE;
	   ngammaspercaptcha.at(whichcaptcha) = theNgam;

	 }
       }
       else{ //G4cout<<"Non Gamma Capture Particle: "<<partcode<<" "<<partname<<G4endl; 
       }
     }
     //Done with neutroncapture stats
     //G4cout << "done with neutroncapture stats" << G4endl;
     WCLiteTrajectoryPoint* mypntP;

     //Full track info on all non-photon particles
     //G4int numloops=0;
     if(partcode!=100){
       for(int i=0; i<=numPoints-1; i++){
       	 numloops++;
       	 //G4cout << "loop " << i << " of " << numPoints << G4endl;

	 mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(i));
	 if (mypnt==0){G4cout << "mypnt null!" <<G4endl;}

	 xStep = mypnt->GetPosition().x();
	 yStep = mypnt->GetPosition().y();
	 zStep = mypnt->GetPosition().z();
	 tStep = mypnt->GetGlobalTime();
	 TEStep = mypnt->GetTotalEnergy();
	 KEStep = mypnt->GetKineticEnergy();
	 
	 iStep = i;

	 TString theProcess = mypnt->GetProcessName();
	 processStep=ConvertProcessNameToCode(theProcess);
	 
	 if(i<(numPoints-1)){
	   mypntP = (WCLiteTrajectoryPoint*)(trj->GetPoint(i+1));
	   if (mypntP==0){G4cout << "mypntP null!" <<G4endl;}	   
	   double xStepP = mypntP->GetPosition().x();
	   double yStepP = mypntP->GetPosition().y();
	   double zStepP = mypntP->GetPosition().z();
	   double tStepP = mypntP->GetGlobalTime();
	   double TEStepP = mypntP->GetTotalEnergy();
	   double KEStepP = mypntP->GetKineticEnergy();

	   double StepL,TrackL,ProcessStepSub,TrackStatus,TEdep;
	   StepL=0.;
	   TrackL=0.;
	   TEdep=TEStepP-TEStep;
	   ProcessStepSub=0.;
	   TrackStatus=0.;


	   //	   (*textout)<<eventcount<<" 0 0 0 0 "<<iStep<<" "<<partname<<" "<<partcode<<" "<<trackid<<" "<<parentid<<" "<<xStep<<" "<<yStep<<" "<<zStep<<" "<<xStepP<<" "<<yStepP<<" "<<zStepP<<" "<<tStep<<" "<<tStepP<<" "<<KEStep<<" "<<KEStepP<<" "<<TEdep<<" "<<StepL<<" "<<TrackL<<" "<<theProcess<<" "<<processStep<<" "<<ProcessStepSub<<" "<<TrackStatus<<" 0 "<<std::endl;

	 }



	 //	 if(iStep==0) trk->Fill();
       }
     }
     // Done with Full track info on all non-photon particles
     //G4cout << "done looping all non-photons" << G4endl;
  } // Done looping over all trajectories
  //G4cout << "A total of " << numloops << " pointless loops were executed." << G4endl;	
  	
   // G4cout<<ncapturecount<<" captures "<<xofcaptcha.size()<<G4endl;

   if(ncapturecount>kcapmax) G4cout<<"max number of captures in capt array (100) has been exeeded"<<G4endl;

   // now we begin final neutron count
   for(int nunu=0; nunu<neutron_xStart.size(); nunu++){

     int ismatched=0;

     for(int npnp=0; npnp<neutron_xStart.size(); npnp++){

       //std::cout << "Now checking neutron tracks..." << std::endl;

       if( (neutron_xStart.at(nunu)==neutron_xEnd.at(npnp)) && 
	   (neutron_yStart.at(nunu)==neutron_yEnd.at(npnp)) && //Conditional jump or move depends on uninitialised value(s) ??
	   (neutron_zStart.at(nunu)==neutron_zEnd.at(npnp)) ){
	   ismatched=1;
	   
	   //std::cout << "NEUTRON TRACK MATCHED!" << std::endl;

	 }
       }

     //std::cout << "Incrementing neutron count now..." << std::endl;

     if(ismatched==0) neutroncount++;
   }

   for(int cc=0; cc<ncapturecount; cc++){	// i think this just copies from vectors into the arrays used for filling ROOT tree
     //G4cout<<cc<<G4endl;
     capt_num[cc] = cc;
     capt_x[cc] = xofcaptcha.at(cc);
     capt_y[cc] = yofcaptcha.at(cc);
     capt_z[cc] = zofcaptcha.at(cc);
     capt_t0[cc] = t0ofcaptcha.at(cc);
     capt_ngamma[cc] = (G4int)ngammaspercaptcha.at(cc);
     capt_E[cc] = energypercaptcha.at(cc);
 
     capt_nucleus[cc]=0;
  

     //     if( (capt_E[cc]>7.5) && (capt_E[cc]<9.0) ) capt_nucleus[cc] = 1;
     //     if( (capt_E[cc]>2.0) && (capt_E[cc]<2.5) ) capt_nucleus[cc] = 2;

     //G4cout<<"OK..THIS IS CAPTURE "<<cc<<":"<<G4endl;
    
     int nnucs=0;
     for(int lll=0; lll<t0ofcaptcha_nuc.size(); lll++){	// this loop checks that we only got one daughter non-gamma from a nCapture.
       

       //       if((capt_t0[cc]-t0ofcaptcha_nuc.at(lll))<5.0){
	 //	 G4cout<<"non-gamma daughter particle: "<<pidofcaptcha_nuc.at(lll)<<G4endl;
       if(parentidofcaptcha_nuc.at(lll)==parentidofcaptcha.at(cc)){
	 if( pidofcaptcha_nuc.at(lll)!=11 ){
	   nnucs++;
	   capt_nucleus[cc] = pidofcaptcha_nuc.at(lll);
	 }
       }
     }
    
     if(nnucs>1) std::cout<<"NOOOOOOOOO!!! NNUCS IS GT 1!!!!"<<std::endl;


     capt_nphot[cc]=0;
   }
   
   // G4cout<<"done "<<G4endl;

   for(int pc=0; pc<nphot; pc++){

     phot_capnum[pc]=-1;

     //     if(pc%1000==0) G4cout<<"pc "<<pc<<G4endl;

     if(pc>1000000) G4cout<<"MAX number of photons in phot array (1,000,000) has been exeeded!!!"<<G4endl;

     for(int cc=0; cc<ncapturecount; cc++){
       if( (phot_tStart[pc]-capt_t0[cc])<10 ){
	 phot_capnum[pc]=cc;
	 capt_nphot[cc]++;
       } 
     }
   }
   evttree->Fill();
   G4cout <<"*********************************************"<<G4endl;

   /*
   for(int jjj=0; jjj<ncapturecount; jjj++){
     G4cout<<"Capture #"<<jjj+1<<" t0: "<<t0ofcaptcha.at(jjj)<<" ngammas: "<<ngammaspercaptcha.at(jjj)<<" energy: "<<energypercaptcha.at(jjj)<<G4endl;

     //capturecount=jjj;
     totalcapturecount++;
     captenergy = energypercaptcha.at(jjj);
     ncaptgammas = ngammaspercaptcha.at(jjj);
     captT = t0ofcaptcha.at(jjj);
     captx = xofcaptcha.at(jjj);  
     capty = yofcaptcha.at(jjj);  
     captz = zofcaptcha.at(jjj);  

     G4cout<< captT <<"N phots: "<<photons.size()<<G4endl;

 
     nphotspercapt=0;
     for(int iii=0; iii<photons.size(); iii++){

       std::vector<double> ahit;
       ahit = photons.at(iii);
       if( fabs(ahit.at(3) - captT)<35. ){
	 cxStart = ahit.at(0);
	 cyStart = ahit.at(1);
	 czStart = ahit.at(2);
	 ctStart = ahit.at(3);
	 cxEnd = ahit.at(4);
	 cyEnd = ahit.at(5);
	 czEnd = ahit.at(6);
	 ctEnd = ahit.at(7);
	 cwavelengthStart = ahit.at(8);
	 cisScatPhot = ahit.at(9);
	 cparentid = ahit.at(10);
	 cprocessStart = ahit.at(11);
	 //	 capt->Fill();
	 nphotspercapt++;
       }
     }

     //     gcapt->Fill();

     }*/

   if(n_earlyMuons==1){
 
     //    G4cout<<"here2 "<<mxStart<<" "<<myStart<<G4endl;
     totalmucount++;

     G4int nearlycount=0;

     for(int iii=0; iii<photons.size(); iii++){

       std::vector<double> ahit;
       ahit = photons.at(iii);
       if( ahit.at(7) < 18. ){
	 //	 mhxStart = ahit.at(0);
	 //	 mhyStart = ahit.at(1);
	 //	 mhzStart = ahit.at(2);
	 //	 mhtStart = ahit.at(3);
	 //	 mhxEnd = ahit.at(4);
	 //	 mhyEnd = ahit.at(5);
	 //	 mhzEnd = ahit.at(6);
	 //	 mhtEnd = ahit.at(7);
	 //	 mhwavelength = ahit.at(8);
	 //	 mhisScatPhot = ahit.at(9);
	 //	 mhparentid = ahit.at(10);
	 //	 mhprocessStart = ahit.at(11);
	 //	 mul->Fill();
	 nearlycount++;
       }
     }

     //     if(nearlycount>0) gtrk->Fill();

   }
   
  //********* HANDLING MRD HITS *************
  //=========================================
  // Get the ID of the HitCollection we want, based on it's name (if we don't already know the ID(?)).  
  G4SDManager* pSDMan = G4SDManager::GetSDMpointer();
  G4int collectionID = pSDMan->GetCollectionID("mrdHitsCollection");	
  
  // Get current event from event manager, then retrieve desired hit collection for that event
  G4RunManager* pRunMan = G4RunManager::GetRunManager();
  const G4Event* currentEvent = pRunMan->GetCurrentEvent();
  G4HCofThisEvent* HCofEvent = currentEvent->GetHCofThisEvent();
  MRDHitsCollection* myCollection = (MRDHitsCollection*)(HCofEvent->GetHC(collectionID));
  
  G4double totalEnergy = 0.;
  G4int numberOfHits = myCollection->GetSize();
  G4cout << "A total of " << numberOfHits << " hits on MRD were recorded!" << G4endl;
  if (myCollection!=0) {
    for (hitnum=0; hitnum<numberOfHits; hitnum++) {
       MRDHit* aHit = (*myCollection)[hitnum];
       totalEnergy += aHit->GetHitEdeposit();
       G4ThreeVector hitPos = aHit->GetHitPos();
       hitPosx=hitPos.x();
       hitPosy=hitPos.y();
       hitPosz=hitPos.z();
       hitTime = aHit->GetHitTime();
       hitParticleName = aHit->GetHitParticleName();
       hitTrackID = aHit->GetHitTrackID(); 
       hitPartCode = ConvertParticleNameToCode(hitParticleName); // or use hitParticleID - from PDGEncoding?
       //hitPartCode = aHit->GetHitParticleID();
       hitProcessName = aHit->GetHitProcessName();
       hitProcessCode = ConvertProcessNameToCode(hitProcessName);
       hitEdep = aHit->GetHitEdeposit();
       hitCopyNum = aHit->GetHitCopyNum();
       hitPhysical = (std::string)aHit->GetHitPhysical();
       // aHit->Print();
       
       mrdhit_x[hitnum] = hitPosx;
       mrdhit_y[hitnum] = hitPosy;
       mrdhit_z[hitnum] = hitPosz;
       mrdhit_t[hitnum] = hitTime;
       mrdhit_process[hitnum] = hitProcessCode;
       mrdhit_particleID[hitnum] = hitPartCode;
       mrdhit_trackID[hitnum] = hitTrackID;
       mrdhit_edep[hitnum] = hitEdep;
       mrdhit_copynum[hitnum] = hitCopyNum;
       if(hitPhysical=="mrdPaddles"){objnum=0;}	// all hits have this object num.. why...
       else if(hitPhysical=="mrdPlates"){objnum=1;}
       else if(hitPhysical.find("alu")!=0){objnum=2;} 	//found!=std::string::npos
       else{objnum=-1;}
       mrdhit_objnum[hitnum] = objnum;
    }
    hitnum=numberOfHits;
  }
  mrdtree->Fill();
  G4cout<<"Energy deposited in MRD: "<<totalEnergy*MeV<<" MeV"<<G4endl;
     

  //******* DONE HANDLING MRD HITS **********
  //=========================================
  
  //************ HANDLING MRD PMT DETECTIONS ***********
  //====================================================
  // Get the ID of the HitCollection we want, based on it's name (if we don't already know the ID(?)).  
  collectionID = pSDMan->GetCollectionID("mrdpmtHitCollection");	
  mrdPMThitsCollection* myPMTCollection = (mrdPMThitsCollection*)(HCofEvent->GetHC(collectionID));
  
  G4int numberOfPMThits = myPMTCollection->GetSize();
  G4cout << "A total of " << numberOfPMThits << " hits on MRD LGs were recorded!" << G4endl;
  if (myPMTCollection!=0) {
    for (pmt_hitnum=0; pmt_hitnum<numberOfPMThits; pmt_hitnum++) {
       //G4cout<<"adding PMT hit"<<G4endl;
       if(pmt_hitnum>kpmthitnmax){G4cout<<"max num PMT hits exceeded!"<<G4endl;break;}
       mrdPMThit* aHit = (*myPMTCollection)[pmt_hitnum];
       hitTrackID = aHit->GetHitTrackID();
       hitParentID = aHit->GetHitParentID();
       hitPMTnumber = aHit->GetPMTNumber();
       G4ThreeVector hitPos = aHit->GetHitPos();
       hitPosx=hitPos.x();
       hitPosy=hitPos.y();
       hitPosz=hitPos.z();
       hitTime = aHit->GetHitTime();
       hitProcessName = aHit->GetCreationProcess();
       hitProcessCode = ConvertProcessNameToCode(hitProcessName); // need to generate this function      
       hitWavelength = aHit->GetHitWavelength();

       // aHit->Print();  
  
       pmthit_x[pmt_hitnum] = hitPosx;
       pmthit_y[pmt_hitnum] = hitPosy;
       pmthit_z[pmt_hitnum] = hitPosz;
       pmthit_t[pmt_hitnum] = hitTime;
       pmthit_process[pmt_hitnum] = hitProcessCode;
       pmthit_trackID[pmt_hitnum] = hitTrackID;
       pmthit_parentID[pmt_hitnum] = hitParentID;
       pmthit_wavelength[pmt_hitnum] = hitWavelength;
       pmthit_copynum[pmt_hitnum] = hitPMTnumber;

    }
    pmt_hitnum=numberOfPMThits;
  }
  pmttree->Fill();
     

  //********* DONE HANDLING MRD PMT DETECTIONS *********
  //====================================================
  
   G4cout <<"event: "<<eventcount<<G4endl;
   G4cout <<"Number of particles is: "<<n_trajectories<<G4endl;
   G4cout <<"Number of neutrons is: "<<neutroncount<<G4endl;
   G4cout <<"My neutron count is: "<<myneutroncount<<G4endl;
   G4cout<<ncapturecount<<" captures"<<G4endl;
   G4cout <<nphot<<" Photons"<<G4endl;
   G4cout <<"******************************************"<<G4endl;
   G4cout <<"******************************************"<<G4endl;
   G4cout <<"END EndOfEventAction"<<G4endl;
   G4cout <<" "<<G4endl;
   G4cout <<" "<<G4endl;   
   G4cout <<"EVENT STARTED AT: "<<eventstarttime<<" AND ENDED AT "<<eventduration<<G4endl;
   eventcount++;
   
  // Trajectory drawing now done by vis mananger under vis comamnds.
  // See vis.mac.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int ConvertParticleNameToCode(TString particleName){
    G4int particleCode = -1;
    if(particleName=="opticalphoton") particleCode=100;
    if(particleName=="e+") particleCode=-11;
    if(particleName=="e-") particleCode=11;
    if(particleName=="gamma") particleCode=22;
    if(particleName=="mu+") particleCode=-13;
    if(particleName=="mu-") particleCode=13;
    if(particleName=="pi0") particleCode=111;
    if(particleName=="pi+") particleCode=211;
    if(particleName=="pi-") particleCode=-211;
    if(particleName=="neutron") particleCode = 2112 ;
    if(particleName=="proton") particleCode = 2212 ;
    if(particleName=="nu_mu") particleCode=14;
    if(particleName=="nu_e") particleCode=12;
    if(particleName=="anti_nu_e") particleCode=-12;
    if(particleName=="anti_nu_mu") particleCode=-14;
    if(particleName=="alpha") particleCode=3328;
    if(particleName=="deuteron") particleCode=3329;
    if(particleName=="triton") particleCode=3330;
    if(particleName=="Li7[0.0]") particleCode=3351;
    if(particleName=="C10[0.0]") particleCode=3331;
    if(particleName=="B11[0.0]") particleCode=3345;
    if(particleName=="C12[0.0]") particleCode=3332;
    if(particleName=="C13[0.0]") particleCode=3350;
    if(particleName=="N13[0.0]") particleCode=3349;
    if(particleName=="N14[0.0]") particleCode=3340;
    if(particleName=="N15[0.0]") particleCode=3333;
    if(particleName=="N16[0.0]") particleCode=3334;
    if(particleName=="O16[0.0]") particleCode=3335;
    if(particleName=="Al27[0.0]") particleCode=3346;
    if(particleName=="Fe54[0.0]") particleCode=3341;
    if(particleName=="Mn54[0.0]") particleCode=3348;    
    if(particleName=="Mn55[0.0]") particleCode=3342;
    if(particleName=="Mn56[0.0]") particleCode=3352;
    if(particleName=="Fe56[0.0]") particleCode=3343;
    if(particleName=="Fe57[0.0]") particleCode=3344;
    if(particleName=="Fe58[0.0]") particleCode=3347; 
    if(particleName=="Eu154[0.0]") particleCode=3353;
    if(particleName=="Gd158[0.0]") particleCode=3336;
    if(particleName=="Gd156[0.0]") particleCode=3337;
    if(particleName=="Gd157[0.0]") particleCode=3338;
    if(particleName=="Gd155[0.0]") particleCode=3339;
    
    if(particleCode==-1) G4cout<<"unaccounted for particle: "<<particleName<<G4endl;
    return particleCode;
}

G4int ConvertProcessNameToCode(TString processName){
    G4int processCode=-1;
    if(processName=="Transportation") processCode=0;
    if(processName=="Scintillation") processCode=1;
    if(processName=="Cerenkov") processCode=2;
    if(processName=="phot") processCode=3;
    if(processName=="OpAbsorption") processCode=4;
    if(processName=="eIoni") processCode=5;
    if(processName=="muIoni") processCode=6;
    if(processName=="hIoni") processCode=7;
    if(processName=="eBrem") processCode=8;
    if(processName=="muBrem") processCode=9;
    if(processName=="HadronElastic") processCode=10;
    if(processName=="nCapture") processCode=11;
    if(processName=="compt") processCode=12;
    if(processName=="Decay") processCode=13;
    if(processName=="muMinusCaptureAtRest") processCode=14;
    if(processName=="NeutronInelastic") processCode=15;
    if(processName=="ProtonInelastic") processCode=16;
    if(processName=="conv") processCode=17;
    if(processName=="annihil") processCode=18;
    if(processName=="PiMinusAbsorptionAtRest") processCode=19;
    if(processName=="PionMinusInelastic") processCode=20;
    if(processName=="PionPlusInelastic") processCode=21;
    if(processName=="Electromagnetic") processCode=22;
    if(processName=="msc") processCode=23;	//msc = multiple scattering
    
    if(processCode<0) G4cout<<"Process unaccounted for "<<processName<<G4endl;
    return processCode;
}

