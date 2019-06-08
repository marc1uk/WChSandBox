#include "WCLiteTrajectory.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

#include <sstream>

//G4Allocator<WCLiteTrajectory> aTrajectoryAllocator;
G4Allocator<WCLiteTrajectory> myTrajectoryAllocator;

WCLiteTrajectory::WCLiteTrajectory()
  :  positionRecord(0), fTrackID(0), fParentID(0),
     PDGEncoding( 0 ), PDGCharge(0.0), ParticleName(""),
     initialMomentum( G4ThreeVector() ),SaveIt(false),creatorProcess(""),
     globalTime(0.0)
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

  stoppingPoint  = seco->stoppingPoint;
  stoppingVolume = seco->stoppingVolume;
  
  G4int ent = seco->GetPointEntries();
  for(G4int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  { 
    positionRecord->push_back((*(seco->positionRecord))[i]);
    //    positionRecord->push_back(seco->positionRecord->removeAt(1));
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();
}


