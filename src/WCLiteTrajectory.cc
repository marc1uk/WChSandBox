#include "WCLiteTrajectory.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "WCLiteTrackInformation.hh"
#include "MRDSD.hh"

#include <sstream>
#include <string>
#include <iostream>

G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;

std::string trackStatuses[] = {"fAlive", "fStopButAlive", "fStopAndKill", "fKillTrackAndSecondaries", "fSuspend", "fPostponeToNextEvent"};
std::string TrackStatusToString(G4TrackStatus status){
	return trackStatuses[(int)status];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WCLiteTrajectory::WCLiteTrajectory()
  :  positionRecord(0), fParentID(0), fTrackID(0), PDGEncoding( 0 ), PDGCharge(0.0), 
     ParticleName(""), initialMomentum( G4ThreeVector() ), SaveIt(false), creatorProcess(""),
     globalTime(0.0), stoppingPoint( G4ThreeVector() ), stoppingVolume(0), 
     sameasparenttrackid(0), numSecondaries(0), GdOriginalTrackID(-1), passesthroughNCV(false), 
     ncvEntryPos( G4ThreeVector() ), ncvEntryTime(0.), ncvExitPos( G4ThreeVector() ), ncvExitTime(0.)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WCLiteTrajectory::WCLiteTrajectory(const G4Track* aTrack)
{
 // ************************************************
 // ================================================
 // The following is for the first trajectory point
 // ================================================
 
  // Fixed trajectory properties set during stepping: set to default for now
  // ========================================================================
  sameasparenttrackid=0;
  passesthroughNCV=0;
  numSecondaries=0;
  GdOriginalTrackID=-1;
  ncvEntryPos = G4ThreeVector();
  ncvEntryTime = 0.;
  ncvExitPos = G4ThreeVector();
  ncvExitTime = 0.;
  stoppingPoint = G4ThreeVector();
  stoppingVolume = 0;
  
  // Fixed trajectory properties from G4Track
  // =========================================
  G4ParticleDefinition * fpParticleDefinition = aTrack->GetDefinition();
  ParticleName = fpParticleDefinition->GetParticleName();
  PDGCharge = fpParticleDefinition->GetPDGCharge();
  PDGEncoding = fpParticleDefinition->GetPDGEncoding();
  fTrackID = aTrack->GetTrackID();
  fParentID = aTrack->GetParentID();
  initialMomentum = aTrack->GetMomentum();
  positionRecord = new TrajectoryPointContainer();
  
  // Properties for the first AppendStep() call
  // ==========================================
  globalTime = aTrack->GetGlobalTime();
  stoppingPoint  = aTrack->GetPosition();
  stoppingVolume = aTrack->GetVolume();
  
  G4String theLocation = "Default";
  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="ExitWorld";
  }
  
  if(theLocation=="NCVliquid_phys"){
    ncvEntryPos = aTrack->GetPosition();
    ncvEntryTime=aTrack->GetGlobalTime();
  }
  
  if ( aTrack->GetUserInformation() != 0 ){
    SaveIt = true;
  } else {
    SaveIt = false;
  }
  
  G4String theProcess;
  if (aTrack->GetCreatorProcess() != 0 ) {
      const G4VProcess* tempproc = aTrack->GetCreatorProcess();
      creatorProcess = tempproc->GetProcessName();
      if(creatorProcess=="nCapture"){
        GdOriginalTrackID=fParentID;
      }
  } else {
    creatorProcess = "Primary";
    //G4cout<<"creation of track "<<fTrackID<<" for "<<ParticleName<<" with null creator process and partentid "<<fParentID<<", energy "<<aTrack->GetKineticEnergy()<<", momentum "<<initialMomentum<<G4endl;
  }
//  if(fParentID==0){
//  	G4cout<<"creation of primary track "<<fTrackID<<" for "<<ParticleName<<G4endl;
//  }
  
  positionRecord->push_back(
     new WCLiteTrajectoryPoint(
        aTrack->GetPosition(),
        aTrack->GetMomentumDirection(),
        aTrack->GetGlobalTime(),
        aTrack->GetTotalEnergy(),
        aTrack->GetKineticEnergy(),
        creatorProcess,
        theLocation
     )
  );
  
 // ================================================
 // End setting values for first trajectory point
 // ================================================
  // ************************************************
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WCLiteTrajectory::WCLiteTrajectory(WCLiteTrajectory & right):G4VTrajectory()
{
// ================
// Copy Constructor - probably never gets called
// ================
  ParticleName = right.ParticleName;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  initialMomentum = right.initialMomentum;
  positionRecord = new TrajectoryPointContainer();
  
  stoppingPoint  = right.stoppingPoint;
  stoppingVolume = right.stoppingVolume;
  SaveIt = right.SaveIt;
  creatorProcess = right.creatorProcess;
  
  sameasparenttrackid=right.sameasparenttrackid;
  numSecondaries=right.numSecondaries;
  GdOriginalTrackID=right.GdOriginalTrackID;
  passesthroughNCV=right.passesthroughNCV;  
  ncvEntryPos=right.ncvEntryPos;
  ncvEntryTime=right.ncvEntryTime;
  ncvExitPos=right.ncvExitPos;
  ncvExitTime=right.ncvExitTime;
  
  for(size_t i=0;i<right.positionRecord->size();i++)
  {
    WCLiteTrajectoryPoint* rightPoint = (WCLiteTrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->push_back(new WCLiteTrajectoryPoint(*rightPoint));
  }
  globalTime = right.globalTime;
  
  // =====================
  // End Copy Constructor
  // =====================
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WCLiteTrajectory::~WCLiteTrajectory()
{
  //  positionRecord->clearAndDestroy();
  size_t i;
  for(i=0;i<positionRecord->size();i++){
    delete  (*positionRecord)[i];
  }
  positionRecord->clear();
  
  delete positionRecord;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::ShowTrajectory(std::ostream& os) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::ShowTrajectory(os);
  // ... or override with your own code here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::DrawTrajectory(G4int i_mode) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::DrawTrajectory();
  // ... or override with your own code here.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::map<G4String,G4AttDef>* WCLiteTrajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("WCLiteTrajectory",isNew);
  if (isNew) {
    
    G4String ID("ID");
    (*store)[ID] = G4AttDef(ID,"Track ID","Bookkeeping","","G4int");
    
    G4String PID("PID");
    (*store)[PID] = G4AttDef(PID,"Parent ID","Bookkeeping","","G4int");
    
    G4String PN("PN");
    (*store)[PN] = G4AttDef(PN,"Particle Name","Physics","","G4String");
    
    G4String Ch("Ch");
    (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","","G4double");
    
    G4String PDG("PDG");
    (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Physics","","G4int");
    
    G4String IMom("IMom");
    (*store)[IMom] = G4AttDef(IMom, "Momentum of track at start of trajectory",
			      "Physics","","G4ThreeVector");   
    G4String NTP("NTP");
    (*store)[NTP] = G4AttDef(NTP,"No. of points","Physics","","G4int");
    
  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4AttValue>* WCLiteTrajectory::CreateAttValues() const
{
  char c[100];
  //std::ostrstream s(c,100);
  std::ostringstream s(c);
  
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
  
  s.seekp(std::ios::beg);
  s << fTrackID << std::ends;
  values->push_back(G4AttValue("ID",c,""));
  
  s.seekp(std::ios::beg);
  s << fParentID << std::ends;
  values->push_back(G4AttValue("PID",c,""));
  
  values->push_back(G4AttValue("PN",ParticleName,""));
  
  s.seekp(std::ios::beg);
  s << PDGCharge << std::ends;
  values->push_back(G4AttValue("Ch",c,""));

  s.seekp(std::ios::beg);
  s << PDGEncoding << std::ends;
  values->push_back(G4AttValue("PDG",c,""));

  s.seekp(std::ios::beg);
  s << G4BestUnit(initialMomentum,"Energy") << std::ends;
  values->push_back(G4AttValue("IMom",c,""));

  s.seekp(std::ios::beg);
  s << GetPointEntries() << std::ends;
  values->push_back(G4AttValue("NTP",c,""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::AppendStep(const G4Step* aStep)
{
  G4Track* aTrack = aStep->GetTrack();

  // get the location of the next step
  G4String theLocation;
  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
    //G4cout<<"AppendStep: DetectorHit"<<G4endl;
  }
  
  // get the step process
  G4String theProcess;
  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
    if(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!=""){
    	 theProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    } else {
    	theProcess = "MysteryProcess";
    }
  } else {
    theProcess+="UserLimit";
  }
  
  G4TrackStatus currentstatus = aTrack->GetTrackStatus();		// fStopButAlive if it has?
  G4int numsecondaries = aStep->GetNumberOfSecondariesInCurrentStep();
  G4double depositedenergythisstep = aStep->GetTotalEnergyDeposit();
  
  // check if this step puts the track in the NCV
  G4VSensitiveDetector* sensdet = aStep->GetPostStepPoint()->GetSensitiveDetector();
  G4String sdetectorname;
  if(sensdet!=0){
    sdetectorname=sensdet->GetName(); 
  } else {
    sdetectorname = "";
  }

		// if we have a need for user information...
		if(sdetectorname=="NeutronCaptureVolume"){
			// check if the track has user information already
			WCLiteTrackInformation* userinfo = (WCLiteTrackInformation*)aTrack->GetUserInformation();
			if(userinfo==0){
				// no user info: make it first
				userinfo = new WCLiteTrackInformation();
				userinfo->SetIsPrimaryN(false);
				aTrack->SetUserInformation(userinfo);
			}
			// Then record those things that needed it: ensure passage through NCV is noted if needed
//			if(sdetectorname=="MuonRangeDetector"){
//				//G4cout<<"particle step in MRD: track zpos: " << aTrack->GetPosition().z();
//				//G4cout<<" particle type: "<< aTrack->GetParticleDefinition()->GetParticleName()<<G4endl;
//				userinfo->SetPassThruMRD(1);
//			}
			if(sdetectorname=="NeutronCaptureVolume"){
				userinfo->SetPassThruNCV(1);
				userinfo->SetNCVentryPos(aStep->GetPostStepPoint()->GetPosition());
				userinfo->SetNCVentryTime(aStep->GetPostStepPoint()->GetGlobalTime());
				//userinfo->SetNCVentryEnergy(aStep->GetPostStepPoint()->GetTotalEnergy()); ?? worth implmenting?
			}
		}
  
  
  if(sdetectorname=="MRDPMTSD"){
  	G4cout<<"Step in the LG! "<<aTrack->GetParticleDefinition()->GetParticleName()<<G4endl;
  }
  
  
  // append the step to the trajectory
  /* JUST PRINTING
  if(aTrack->GetDefinition()->GetParticleName()=="neutron"){
  
 	G4double pntx = aStep->GetPostStepPoint()->GetPosition().x();
 	G4double pnty = aStep->GetPostStepPoint()->GetPosition().y();
 	G4double pntz = aStep->GetPostStepPoint()->GetPosition().z();
 	G4int trackiddd = aTrack->GetTrackID();
 	G4int numpoints=this->GetPointEntries();
  	G4cout << "Appending point (" << pntx << ", " <<  pnty << ", " << pntz << ") from trackID " << trackiddd;
  	G4cout << " to trajectory " << this << " which currently has " << numpoints << " points." << G4endl;
  } 
  */
  
  /*
  G4String lastprocess="noexistingpoints";
  G4String secondlastprocess="noexistingpoints";
  G4int existingpoints = positionRecord->size();
  if(existingpoints>0){
  	WCLiteTrajectoryPoint* lastpoint = (WCLiteTrajectoryPoint*)positionRecord->at(existingpoints-1);
  	lastprocess = lastpoint->GetProcessName();
  }
  if(existingpoints>1){
    WCLiteTrajectoryPoint* seclastpoint = (WCLiteTrajectoryPoint*)positionRecord->at(existingpoints-2);
  	secondlastprocess = seclastpoint->GetProcessName();
  }
  G4String aParticleType = aTrack->GetDefinition()->GetParticleName();
  
  	if(theProcess=="Scintillation"&&!(theLocation=="waterTank"||theLocation=="Hall"||theLocation=="DetectorHit")){
		G4cout<<"scintillation in "<<theLocation
		<<" by "<<aParticleType
		<<" with edep "<<(depositedenergythisstep/eV)
		<<" generating "<<numsecondaries<<" secondaries,"
		<<" with status "<<TrackStatusToString(currentstatus)
		<<" at traj point "<<existingpoints
		<<" with last process "<<lastprocess
		<<" and before that "<<secondlastprocess
		<<" with KE "<<aStep->GetPostStepPoint()->GetKineticEnergy()
		<<" and position ("<<aStep->GetPostStepPoint()->GetPosition()[0]<<","<<aStep->GetPostStepPoint()->GetPosition()[1]<<","<<aStep->GetPostStepPoint()->GetPosition()[2]<<")"
		<<G4endl;
	}
	*/
  
  // push this point back into the trajectory - but ignore steps with process Scintillation, and no secondaries or energy deposit
  if(!(theProcess=="Scintillation"&&numsecondaries==0&&depositedenergythisstep==0)){
	  positionRecord->push_back( 
		   new WCLiteTrajectoryPoint(
		   aStep->GetPostStepPoint()->GetPosition(), 
		   aStep->GetPostStepPoint()->GetMomentumDirection(),
		   aStep->GetPostStepPoint()->GetGlobalTime(), 
		   aStep->GetPostStepPoint()->GetTotalEnergy(), 
		   aStep->GetPostStepPoint()->GetKineticEnergy(),
		   theProcess,
		   theLocation
		   )
	  );
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ParticleDefinition* WCLiteTrajectory::GetParticleDefinition()
{
  return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
// ================================================================================================
// Merge sometimes gets called by Geant4 when an inelastic scattering gives a particle two 'tracks'
// ================================================================================================
  if(!secondTrajectory) return;

  WCLiteTrajectory* seco         = (WCLiteTrajectory*)secondTrajectory;

  //fTrackID			// all acquired from Track which is maintained.
  //fParentID
  //PDGEncoding
  //PDGCharge
  //ParticleName
  //creatorProcess
  //globalTime
  G4ThreeVector emptypos = G4ThreeVector();
  
  if(SaveIt||seco->GetSaveFlag()){SaveIt=true;} 			// is this still needed?
  initialMomentum = seco->GetInitialMomentum();
  sameasparenttrackid   =  seco->GetSameAsParentTrackID();
  numSecondaries += (seco->GetNumSecs());
  GdOriginalTrackID = seco->GetGdOriginalTrackID();
  if(this->GetPassThruNCV() ||  seco->GetPassThruNCV()){passesthroughNCV=1;}
  if(this->GetNCVentryPos()==emptypos){ncvEntryPos=seco->GetNCVentryPos(); ncvEntryTime = seco->GetNCVentryTime();}
  if(this->GetNCVexitPos()==emptypos){ncvExitPos=seco->GetNCVexitPos(); ncvExitTime = seco->GetNCVexitTime();}
  stoppingPoint  	= 	seco->stoppingPoint;
  stoppingVolume 	= 	seco->stoppingVolume;
  
  /*if(fTrackID==44953){
	  G4cout << "merging : " << this << " with " << seco << G4endl;
	  G4cout << this << " currently has " << GetPointEntries() << " points, " << seco << " has " << seco->GetPointEntries();
	  G4cout << " additional points." << G4endl;
  }*/
  G4int ent = seco->GetPointEntries();
  for(G4int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    positionRecord->push_back((*(seco->positionRecord))[i]);
    //    positionRecord->push_back(seco->positionRecord->removeAt(1));
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();
  
  /*if(fTrackID==44953){
	G4cout << this << " now has " << GetPointEntries() << " points. " << G4endl;
	WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(GetPoint(0));
	G4double pntx = firstpnt->GetPosition().x();
	G4double pnty = firstpnt->GetPosition().y();
	G4double pntz = firstpnt->GetPosition().z();
	G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
	WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(GetPoint(GetPointEntries()-1));
	pntx = lastpnt->GetPosition().x();
	pnty = lastpnt->GetPosition().y();
	pntz = lastpnt->GetPosition().z();
	G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
  }*/
  
  // ============================
  // End Merge Trajectories
  // ============================
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::CopyInfo(WCLiteTrackInformation* infoin){
// =================================================================================
// Copy Information from UserTrackInformation into Trajectory when Finished Tracking
// =================================================================================
  //G4cout << "copying from: " << infoin << G4endl;
  //infoin->Print();
  //this->Print();
  
  sameasparenttrackid=infoin->GetSameAsParentTrackID();
  numSecondaries=infoin->GetNumSecs();
  GdOriginalTrackID=infoin->GetGdOriginalTrackID();
  passesthroughNCV=infoin->GetPassThruNCV();
  ncvEntryPos=infoin->GetNCVentryPos();
  ncvEntryTime=infoin->GetNCVentryTime();
  ncvExitPos=infoin->GetNCVexitPos();
  ncvExitTime=infoin->GetNCVexitTime();

// =================================================================================
// End Copying from UserTrackInformation
// =================================================================================
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WCLiteTrajectory::Print() const
{
    G4cout<< "fTrackID = " << fTrackID << G4endl;	
    G4cout<< "fParentID = " << fParentID << G4endl;
    G4cout<< "PDGEncoding = " << PDGEncoding << G4endl;
    G4cout<< "PDGCharge = " << PDGCharge << G4endl;
    G4cout<< "ParticleName = " << ParticleName << G4endl;
    G4cout<< "initialMomentum = " << initialMomentum << G4endl;
    G4cout<< "stoppingPoint = " << stoppingPoint << G4endl;
    //G4cout<< "stoppingVolume = " << stoppingVolume->GetName() << G4endl;
    G4cout<< "SaveIt = " << SaveIt << G4endl;
    G4cout<< "creatorProcess = " << creatorProcess << G4endl;
    G4cout<< "globalTime = " << globalTime << G4endl; 
    G4cout<< "sameasparenttrackid = " << sameasparenttrackid << G4endl;
    G4cout<< "numSecondaries = " << numSecondaries << G4endl;
    G4cout<< "GdOriginalTrackID = " << GdOriginalTrackID << G4endl;
    G4cout<< "passesthroughNCV = " << passesthroughNCV << G4endl;
    G4cout<< "ncvEntryPos = " << ncvEntryPos << G4endl;
    G4cout<< "ncvEntryTime = " << ncvEntryTime << G4endl;
    G4cout<< "ncvExitPos = " << ncvExitPos << G4endl;
    G4cout<< "ncvExitTime = " << ncvExitTime << G4endl;
    G4cout<< "PositionRecord has " << this->GetPointEntries() << " points." << G4endl;
}
 
