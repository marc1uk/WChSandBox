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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
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
  phot_capnum = new G4int[knphotmax];  /*
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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
WCLiteEventAction::~WCLiteEventAction()
{
   no->cd();
   evttree->Write();
   no->Close();

   textout->close();

   delete phot_xStart;
   delete phot_yStart;
   delete phot_zStart;
   delete phot_tStart;
   delete phot_xEnd;
   delete phot_yEnd;
   delete phot_zEnd;
   delete phot_tEnd;
   delete phot_wavelength;
   delete phot_processStart;
   delete phot_isScat;
   delete phot_parentid;
   delete phot_trackid;
   delete phot_hit;
   delete phot_capnum;
   delete part_xStart;
   delete part_yStart;
   delete part_zStart;
   delete part_tStart;
   delete part_xEnd;
   delete part_yEnd;
   delete part_zEnd;
   delete part_tEnd;
   delete part_KEstart;
   delete part_KEend;
   delete part_pxStart;
   delete part_pyStart;
   delete part_pzStart;
   delete part_pxEnd;
   delete part_pyEnd;
   delete part_pzEnd;
   delete part_processStart;
   delete part_processEnd;
   delete part_parentid;
   delete part_trackid;
   delete part_pid;

   delete capt_num;
   delete capt_nucleus;
   delete capt_pid;
   delete capt_nphot;
   delete capt_ngamma;
   delete capt_x;
   delete capt_y;
   delete capt_z;
   delete capt_t0;
   delete capt_E;

  

   G4cout<<"End of this run. Number of events: "<< eventcount <<" no evts with muons: "<< totalmucount<< " frac of evts with muons: "<<(totalmucount/eventcount)<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void WCLiteEventAction::BeginOfEventAction(const G4Event*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void WCLiteEventAction::EndOfEventAction(const G4Event* anEvent)
{
 
  //  G4cout <<"**********************"<<G4endl;
  //  G4cout <<"**********************"<<G4endl;
  G4cout <<" "<<G4endl;
  G4cout <<" "<<G4endl; 
  G4cout <<"START EndOfEventAction"<<G4endl;

  ncapturecount=0;
  neutroncount=0;   
  nphot=0;
  npart=0;

  std::vector<double> ngammaspercaptcha;
  std::vector<double> nfinalphotspercaptcha;
  std::vector<double> energypercaptcha;
  std::vector<double> t0ofcaptcha,xofcaptcha,yofcaptcha,zofcaptcha,nucofcaptcha;
  std::vector<double> t0ofcaptcha_nuc,xofcaptcha_nuc,yofcaptcha_nuc,zofcaptcha_nuc,pidofcaptcha_nuc;
  std::vector<double> neutron_xStart,neutron_yStart,neutron_zStart,neutron_tStart,neutron_xEnd,neutron_yEnd,neutron_zEnd,neutron_tEnd;
  std::vector<int> neutron_processStart,neutron_processEnd;


  std::vector<int> parentidofcaptcha, parentidofcaptcha_nuc;
  std::vector<std::vector<double> > photons;

  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  
  G4int n_earlyEMparts = 0;
  G4int n_earlyMuons = 0;

  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  G4int photcount=0, mpcount=0, neithercount=0;

 
  std::vector<double> vpidtrack;
  std::vector<double> vtrackid;

  /*
  for (G4int i=0; i<n_trajectories; i++){ 
    vpidtrack.push_back(0);
  }
  */

  for (G4int i=0; i<n_trajectories; i++){ 

    TEStart = -5555; // Initialisation

    G4Trajectory* trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
    
    trackid = trj->GetTrackID();
    parentid = trj->GetParentID();
    
    TString partname = (trj->GetParticleName());
    partcode = -1;
    if(partname=="opticalphoton") partcode=100;
    if(partname=="e+") partcode=-11;
    if(partname=="e-") partcode=11;
    if(partname=="gamma") partcode=22;
    if(partname=="mu+") partcode=-13;
    if(partname=="mu-") partcode=13;
    if(partname=="pi0") partcode=111;
    if(partname=="pi+") partcode=211;
    if(partname=="pi-") partcode=-211;
    if(partname=="neutron") partcode = 2112 ;
    if(partname=="proton") partcode = 2212 ;
    if(partname=="nu_mu") partcode=14;
    if(partname=="nu_e") partcode=12;
    if(partname=="anti_nu_e") partcode=-12;
    if(partname=="anti_nu_mu") partcode=-14;
    if(partname=="alpha") partcode=3328;
    if(partname=="deuteron") partcode=3329;
    if(partname=="triton") partcode=3330;
    if(partname=="C10[0.0]") partcode=3331;
    if(partname=="C12[0.0]") partcode=3332;
    if(partname=="N15[0.0]") partcode=3333;
    if(partname=="N16[0.0]") partcode=3334;
    if(partname=="O16[0.0]") partcode=3335;
    if(partname=="Gd158[0.0]") partcode=3336;
    if(partname=="Gd156[0.0]") partcode=3337;
    if(partname=="Gd157[0.0]") partcode=3338;
    if(partname=="Gd155[0.0]") partcode=3339;
    
    if(partcode==-1) G4cout<<"unaccounted for particle: "<<partname<<G4endl;

    G4int numPoints = trj->GetPointEntries();
    WCLiteTrajectoryPoint* mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(numPoints-1));
    xEnd = mypnt->GetPosition().x();
    yEnd = mypnt->GetPosition().y();
    zEnd = mypnt->GetPosition().z();
    xEndDir = mypnt->GetDirection().x();
    yEndDir = mypnt->GetDirection().y();
    zEndDir = mypnt->GetDirection().x();
    tEnd = mypnt->GetGlobalTime();
    TEEnd = mypnt->GetTotalEnergy();
    KEEnd = mypnt->GetKineticEnergy();
     
    TString theProcess = mypnt->GetProcessName();
    processEnd=-1;
    if(theProcess=="Transportation") processEnd=0;
    if(theProcess=="Scintillation") processEnd=1;
    if(theProcess=="Cerenkov") processEnd=2;
    if(theProcess=="phot") processEnd=3;
    if(theProcess=="OpAbsorption") processEnd=4;
    if(theProcess=="eIoni") processEnd=5;
    if(theProcess=="muIoni") processEnd=6;
    if(theProcess=="hIoni") processEnd=7;
    if(theProcess=="eBrem") processEnd=8;
    if(theProcess=="muBrem") processEnd=9;
    if(theProcess=="HadronElastic") processEnd=10;
    if(theProcess=="nCapture") processEnd=11;
    if(theProcess=="compt") processEnd=12;
    if(theProcess=="Decay") processEnd=13;
    if(theProcess=="muMinusCaptureAtRest") processEnd=14;
    if(theProcess=="NeutronInelastic") processEnd=15;
    if(theProcess=="ProtonInelastic") processEnd=16;
    if(theProcess=="conv") processEnd=17;
    if(theProcess=="annihil") processEnd=18;
    if(theProcess=="PiMinusAbsorptionAtRest") processEnd=19;
    if(theProcess=="PionMinusInelastic") processEnd=20;
    if(theProcess=="PionPlusInelastic") processEnd=21;
    
    if(processEnd<0) G4cout<<"Process unaccounted for "<<theProcess<<G4endl;

     phhit=0;
     if(partcode==100 && mypnt->GetDetectorLocation()=="DetectorHit"){
       // G4cout<<mypnt->GetDetectorLocation()<<G4endl;
       phhit=1;
     }

     //   tEnd=5.;

     mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(0));
     xStart = mypnt->GetPosition().x();
     yStart = mypnt->GetPosition().y();
     zStart = mypnt->GetPosition().z();
     xStartDir = mypnt->GetDirection().x();
     yStartDir = mypnt->GetDirection().y();
     zStartDir = mypnt->GetDirection().z();
     tStart = mypnt->GetGlobalTime();
     TEStart = mypnt->GetTotalEnergy();
     KEStart = mypnt->GetKineticEnergy();

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
     processStart=-1;
     if(sProcess=="Transportation") processStart=0;
     if(sProcess=="Scintillation") processStart=1;
     if(sProcess=="Cerenkov") processStart=2;
     if(sProcess=="phot") processStart=3;
     if(sProcess=="OpAbsorption") processStart=4;
     if(sProcess=="eIoni") processStart=5;
     if(sProcess=="muIoni") processStart=6;
     if(sProcess=="hIoni") processStart=7;
     if(sProcess=="eBrem") processStart=8;
     if(sProcess=="muBrem") processStart=9;
     if(sProcess=="HadronElastic") processStart=10;
     if(sProcess=="nCapture") processStart=11;
     if(sProcess=="compt") processStart=12;
     if(sProcess=="Decay") processStart=13;
     if(sProcess=="muMinusCaptureAtRest") processStart=14;
     if(sProcess=="NeutronInelastic") processStart=15;
     if(sProcess=="ProtonInelastic") processStart=16;
     if(sProcess=="conv") processStart=17;
     if(sProcess=="annihil") processStart=18;
     if(sProcess=="PiMinusAbsorptionAtRest") processStart=19;
     if(sProcess=="PionMinusInelastic") processStart=20;
     if(sProcess=="PionPlusInelastic") processStart=21;

     if(processStart==-1){
       std::cout << "WARNING: PROCESS START **NOT** DEFINED!!! (" << sProcess << ")" << std::endl;
     }

     if(KEStart>0) {wavelengthStart = (4.13566733*2.99792458 *0.0001 )/KEStart;}

          
     if(processStart==11){
     std::cout<<"VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV"<<std::endl;
       std::cout<<"Particle Name = "<<partname<<std::endl;
	 for(int qq=0; qq<numPoints; qq++){
	   WCLiteTrajectoryPoint* tpnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(qq));
	   std::cout<<qq<<" "<<tpnt->GetProcessName()<<std::endl;
	 }
     std::cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl;
     }


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
       phot_tEnd[nphot]=tEnd;
       phot_wavelength[nphot]=wavelengthStart;
       phot_processStart[nphot]= (G4int)processStart;
       phot_isScat[nphot]=  (G4int)isScatPhot;
       phot_parentid[nphot]= (G4int)parentid;
       phot_trackid[nphot]= (G4int)trackid;
       phot_hit[nphot]= (G4int)phhit;

       nphot++;
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

       //if(npart%5==0) G4cout<<"npart "<<npart<<G4endl;
       
       npart++;
     }

     // Stats on neutrons and neutron capture
     // if(partcode==2112) neutroncount++;

     if(partcode==2112){

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
     }

     if( (processStart==11) ){

       //G4cout<<"a captcha"<<G4endl;
       std::cout << "Neutron capture! Paticle code = " << partcode << std::endl;

       if(partcode!=22){
	 t0ofcaptcha_nuc.push_back(tStart);
	 xofcaptcha_nuc.push_back(xStart);
	 yofcaptcha_nuc.push_back(yStart);
	 zofcaptcha_nuc.push_back(zStart);
	 pidofcaptcha_nuc.push_back(partcode);
	 parentidofcaptcha_nuc.push_back(parentid);
       }

       if(partcode==22){

	 std::cout<<"parent id: "<<parentid<<std::endl;
	 
	 bool isGd=false;

	 int whichcaptcha=-1;
	 for(int m=0; m<t0ofcaptcha.size(); m++){

	   if(((int)parentid)==parentidofcaptcha.at(m)) whichcaptcha=m;

	   //   if( fabs(xStart-xofcaptcha.at(m))<1.0 && fabs(yStart-yofcaptcha.at(m))<1.0 && fabs(zStart-zofcaptcha.at(m))<1.0 ) whichcaptcha=m;
	   //   if( fabs(tStart-t0ofcaptcha.at(m)) <5.0) whichcaptcha=m;
	 }
	 
	 if(whichcaptcha==-1){

	   ncapturecount++;
	   t0ofcaptcha.push_back(tStart);
	   xofcaptcha.push_back(xStart);
	   yofcaptcha.push_back(yStart);
	   zofcaptcha.push_back(zStart);
	   energypercaptcha.push_back(TEStart);
	   ngammaspercaptcha.push_back(1.0);
	   parentidofcaptcha.push_back((int)parentid);
	   
	   std::cout <<"Capture energy = " << TEStart << std::endl;
	   std::cout << std::endl;

	 } else{
	   double theE = energypercaptcha.at(whichcaptcha);
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

     WCLiteTrajectoryPoint* mypntP;

     //Full track info on all non-photon particles
     if(partcode!=100){
       for(int i=0; i<=numPoints-1; i++){

	 mypnt = (WCLiteTrajectoryPoint*)(trj->GetPoint(i));

	 xStep = mypnt->GetPosition().x();
	 yStep = mypnt->GetPosition().y();
	 zStep = mypnt->GetPosition().z();
	 tStep = mypnt->GetGlobalTime();
	 TEStep = mypnt->GetTotalEnergy();
	 KEStep = mypnt->GetKineticEnergy();
	 
	 iStep = i;

	 TString theProcess = mypnt->GetProcessName();
	 processStep=-1;
	 if(theProcess=="Transportation") processStep=0;
	 if(theProcess=="Scintillation") processStep=1;
	 if(theProcess=="Cerenkov") processStep=2;
	 if(theProcess=="phot") processStep=3;
	 if(theProcess=="OpAbsorption") processStep=4;
	 if(theProcess=="eIoni") processStep=5;
	 if(theProcess=="muIoni") processStep=6;
	 if(theProcess=="hIoni") processStep=7;
	 if(theProcess=="eBrem") processStep=8;
	 if(theProcess=="muBrem") processStep=9;
	 if(theProcess=="HadronElastic") processStep=10;
	 if(theProcess=="nCapture") processStep=11;
	 if(theProcess=="compt") processStep=12;
	 if(theProcess=="Decay") processStep=13;
	 if(theProcess=="muMinusCaptureAtRest") processStep=14;
	 if(theProcess=="NeutronInelastic") processStep=15;
	 if(theProcess=="ProtonInelastic") processStep=16;
	 
	 if(i<(numPoints-1)){
	   mypntP = (WCLiteTrajectoryPoint*)(trj->GetPoint(i+1));	   
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
  } // Done looping over all trajectories

     //   G4cout<<ncapturecount<<" captures "<<xofcaptcha.size()<<G4endl;

   if(ncapturecount>kcapmax) G4cout<<"max number of captures in capt array (100) has been exeeded"<<G4endl;

   // now we begin final neutron count
   for(int nunu=0; nunu<neutron_xStart.size(); nunu++){

     int ismatched=0;

     for(int npnp=0; npnp<neutron_xStart.size(); npnp++){

       if( (neutron_xStart.at(nunu)==neutron_xEnd.at(npnp)) && 
	   (neutron_yStart.at(nunu)==neutron_yEnd.at(npnp)) &&
	   (neutron_zStart.at(nunu)==neutron_zEnd.at(npnp)) ){
	   ismatched=1;
	 }
       }
     if(ismatched==0) neutroncount++;
   }


   for(int cc=0; cc<ncapturecount; cc++){
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
     for(int lll=0; lll<t0ofcaptcha_nuc.size(); lll++){
       

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

   G4cout <<"event: "<<eventcount<<G4endl;
   G4cout <<"Number of particles is: "<<n_trajectories<<G4endl;
   G4cout<<ncapturecount<<" captures"<<G4endl;
   G4cout <<nphot<<" Photons"<<G4endl;
   G4cout <<"******************************************"<<G4endl;
   G4cout <<"******************************************"<<G4endl;
   G4cout <<"END EndOfEventAction"<<G4endl;
   G4cout <<" "<<G4endl;
   G4cout <<" "<<G4endl;   
   
   
   eventcount++;
   
  // Trajectory drawing now done by vis mananger under vis comamnds.
  // See vis.mac.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
