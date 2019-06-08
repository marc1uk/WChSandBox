#include "WCLiteTrajectory.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "WCLiteTrackInformation.hh"
#include "MRDSD.hh"

#include <sstream>

//G4Allocator<WCLiteTrajectory> aTrajectoryAllocator;
G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;

WCLiteTrajectory::WCLiteTrajectory()
  :  positionRecord(0), fParentID(0), fTrackID(0), PDGEncoding( 0 ), PDGCharge(0.0), 
     ParticleName(""), initialMomentum( G4ThreeVector() ), SaveIt(false), creatorProcess(""),
     globalTime(0.0), stoppingPoint( G4ThreeVector() ), stoppingVolume(0), 
     sameasparenttrackid(0), passesthroughMRD(0), isMRDprimary(0), 
     numSecondaries(0), mrdOriginalTrackID(0), mrdStartPos( G4ThreeVector() ), mrdDetected(0)
{;}

WCLiteTrajectory::WCLiteTrajectory(const G4Track* aTrack)
{
  G4ParticleDefinition * fpParticleDefinition = aTrack->GetDefinition();
  ParticleName = fpParticleDefinition->GetParticleName();
  PDGCharge = fpParticleDefinition->GetPDGCharge();
  PDGEncoding = fpParticleDefinition->GetPDGEncoding();
  fTrackID = aTrack->GetTrackID();
  fParentID = aTrack->GetParentID();
  initialMomentum = aTrack->GetMomentum();
  positionRecord = new TrajectoryPointContainer();
  sameasparenttrackid=0;
  passesthroughMRD=0;
  isMRDprimary=0;
  numSecondaries=0;
  mrdOriginalTrackID=0;
  mrdStartPos = G4ThreeVector();
  mrdDetected=0;
  stoppingPoint = G4ThreeVector();
  stoppingVolume = 0;
   
 //Following is for the first trajectory point
  G4String theProcess;
  G4String theLocation;

  theLocation = "Default";
 
  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
  stoppingPoint  = aTrack->GetPosition();
  stoppingVolume = aTrack->GetVolume();
  if ( aTrack->GetUserInformation() != 0 ) 
    SaveIt = true;
  else SaveIt = false;
  globalTime = aTrack->GetGlobalTime();
  if (aTrack->GetCreatorProcess() != 0 )
    {
      const G4VProcess* tempproc = aTrack->GetCreatorProcess();
      creatorProcess = tempproc->GetProcessName();
    }
  else 
    creatorProcess = "";
  positionRecord->push_back(new WCLiteTrajectoryPoint(aTrack->GetPosition(),aTrack->GetMomentumDirection(),aTrack->GetGlobalTime(),aTrack->GetTotalEnergy(),aTrack->GetKineticEnergy(),creatorProcess,theLocation));

}

WCLiteTrajectory::WCLiteTrajectory(WCLiteTrajectory & right):G4VTrajectory()
{
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
  passesthroughMRD=right.passesthroughMRD;
  isMRDprimary=right.isMRDprimary;
  numSecondaries=right.numSecondaries;
  mrdOriginalTrackID=right.mrdOriginalTrackID;
  mrdStartPos=right.mrdStartPos;
  mrdDetected=right.mrdDetected;
  
  for(size_t i=0;i<right.positionRecord->size();i++)
  {
    WCLiteTrajectoryPoint* rightPoint = (WCLiteTrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->push_back(new WCLiteTrajectoryPoint(*rightPoint));
  }
  globalTime = right.globalTime;


}

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

void WCLiteTrajectory::ShowTrajectory(std::ostream& os) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::ShowTrajectory(os);
  // ... or override with your own code here.
}

void WCLiteTrajectory::DrawTrajectory(G4int i_mode) const
{
  // Invoke the default implementation in G4VTrajectory...
  G4VTrajectory::DrawTrajectory(i_mode);
  // ... or override with your own code here.
}

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

void WCLiteTrajectory::AppendStep(const G4Step* aStep)
{
  G4String theProcess;
  G4String theLocation;

  G4Track* aTrack = aStep->GetTrack();

  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
    //    G4cout<<"AppendStep: DetectorHit"<<G4endl;
  }

  
  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
    theProcess = ((G4VProcess*) (aStep->GetPostStepPoint()->GetProcessDefinedStep()))->GetProcessName();
  } else {
    theProcess+="UserLimit";
  }  
  
  // if it hasn't already, check if this step puts the track in the MRD and flag its user info if it does. 
  WCLiteTrackInformation* userinfo = (WCLiteTrackInformation*)aTrack->GetUserInformation();
  if(userinfo!=0){
  	if(userinfo->GetPassThruMRD()==0){	// has this track passed through the MRD so far?
  	  G4VSensitiveDetector* sensdet= aStep->GetPostStepPoint()->GetSensitiveDetector();
	    if(sensdet!=0){
	  	  G4String sdetectorname=sensdet->GetName(); 
		    if(sdetectorname=="MuonRangeDetector"){		// does this step put the track in the MRD?
			    //G4cout<<"particle step in MRD: track zpos: " << aTrack->GetPosition().z();
			    //G4cout<<" particle type: "<< aTrack->GetParticleDefinition()->GetParticleName()<<G4endl;
			    userinfo->SetPassThruMRD(1);
			    G4ThreeVector emptypos = G4ThreeVector();
			    if((userinfo->GetMRDstartPos())==emptypos){userinfo->SetMRDstartPos(aTrack->GetPosition());}
			    if((userinfo->GetmrdOriginalTrackID())==0){
				    userinfo->SetIsMRDprimary(1);
				    userinfo->SetmrdOriginalTrackID(aTrack->GetTrackID());
			    }
		    }
	    }
	  } // else: it's already recorded as passing through the MRD, it's fine.
  } else {	// no user info
  	  MRDSD* mrdsds=(MRDSD*)aStep->GetPostStepPoint()->GetSensitiveDetector();
	    if(mrdsds!=0){
	  	  G4String sdetector = mrdsds->GetName();
		    if(sdetector=="MuonRangeDetector"){		// does this step put the track in the MRD?
		  	  WCLiteTrackInformation* newinfo = new WCLiteTrackInformation();
		  	  newinfo->SetIsPrimaryN(false);
		  	  newinfo->SetPassThruMRD(1);
		  	  //G4cout<<"particle step in MRD: track zpos: " << aTrack->GetPosition().z();
			    //G4cout<<" particle type: "<< aTrack->GetParticleDefinition()->GetParticleName()<<G4endl;
		  	  newinfo->SetIsMRDprimary(1);
		  	  newinfo->SetmrdOriginalTrackID(aTrack->GetTrackID());
		  	  newinfo->SetMRDstartPos(aTrack->GetPosition());
		  	  aTrack->SetUserInformation(newinfo);
		    }
	   }  
  }
  
  G4VSensitiveDetector* sensdet2= aStep->GetPostStepPoint()->GetSensitiveDetector();
  if(sensdet2!=0){
  	G4String sdetectorname2=sensdet2->GetName();
  	if(sdetectorname2=="MRDPMTSD"){G4cout<<"Step in the LG! "<<aTrack->GetParticleDefinition()->GetParticleName()<<G4endl;}
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
  positionRecord->push_back( new WCLiteTrajectoryPoint(aStep->GetPostStepPoint()->GetPosition(), 
						       aStep->GetPostStepPoint()->GetMomentumDirection(),
						       aStep->GetPostStepPoint()->GetGlobalTime(), 
						       aStep->GetPostStepPoint()->GetTotalEnergy(), 
						       aStep->GetPostStepPoint()->GetKineticEnergy(),
						       theProcess,
						       theLocation
						       ));
}

G4ParticleDefinition* WCLiteTrajectory::GetParticleDefinition()
{
  return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void WCLiteTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  WCLiteTrajectory* seco         = (WCLiteTrajectory*)secondTrajectory;

  //fTrackID			// all acquired from Track which is maintained.
  //fParentID
  //PDGEncoding
  //PDGCharge
  //ParticleName
  //creatorProcess
  //globalTime
  if(SaveIt||seco->GetSaveFlag()){SaveIt=true;} 			// is this still needed?
  initialMomentum	=	seco->GetInitialMomentum();
  sameasparenttrackid   = 	seco->GetSameAsParentTrackID();
  if(this->GetPassThruMRD() ||  seco->GetPassThruMRD()){passesthroughMRD=1;}
  if(this->GetIsMRDprimary() || seco->GetIsMRDprimary()){isMRDprimary=1;}
  numSecondaries	+=	(seco->GetNumSecs());
  mrdOriginalTrackID	=	seco->GetmrdOriginalTrackID();
  mrdDetected		= 	seco->GetMRDdetected();	// not really necessary since only applicable for photons, which get killed.
  G4ThreeVector emptypos = G4ThreeVector();
  if(this->GetMRDstartPos()==emptypos){mrdStartPos=seco->GetMRDstartPos();}
  
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
}

void WCLiteTrajectory::CopyInfo(WCLiteTrackInformation* infoin){
  //G4cout << "copying from: " << infoin << G4endl;
  //infoin->Print();
  //this->Print();
  
  sameasparenttrackid=infoin->GetSameAsParentTrackID();
  passesthroughMRD=infoin->GetPassThruMRD();
  isMRDprimary=infoin->GetIsMRDprimary();
  numSecondaries=infoin->GetNumSecs();
  mrdOriginalTrackID=infoin->GetmrdOriginalTrackID();
  mrdStartPos=infoin->GetMRDstartPos();
  mrdDetected=infoin->GetMRDdetected();
  //exitProcess=infoin->GetExitProcess();
}

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
    G4cout<< "passesthroughMRD = " << passesthroughMRD << G4endl;
    G4cout<< "isMRDprimary = " << isMRDprimary << G4endl;
    G4cout<< "numSecondaries = " << numSecondaries << G4endl;
    G4cout<< "mrdOriginalTrackID = " << mrdOriginalTrackID << G4endl;
    G4cout<< "mrdStartPos = " << mrdStartPos << G4endl;
    G4cout<< "mrdDetected = " << mrdDetected << G4endl;
    G4cout<< "PositionRecord has " << this->GetPointEntries() << " points." << G4endl;
}
 
