// ====================================================================
//   SBsimMRDResponse.cc
//   
//   2006/09/04 K. Hiraide
// ====================================================================
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "SBsimMRDResponse.hh"
#include "SBsimMRDDB.hh"

/////////////////////////////////////////////////////////////////////////////
SBsimMRDResponse::SBsimMRDResponse(WCLiteDetectorConstruction* detector)	// marcus: replaced SBsimDetectorConstruction
/////////////////////////////////////////////////////////////////////////////
{
  //sbcard  = detector->GetInputCard()->GetSciBarCard(); marcus: hopefully can work around this
  //sbdb    = detector->GetSciBarDB(); marcus: try to work around this
  mrddb   = detector->GetMRDDB();
}

///////////////////////////////////////////
SBsimMRDResponse::~SBsimMRDResponse()
///////////////////////////////////////////
{

}

//////////////////////////////////////////////////////////////////////////////////
void SBsimMRDResponse::ApplyScintiResponse(G4double* edeposit, G4Track* aTrack
					   ,G4StepPoint* preStepPoint)
//////////////////////////////////////////////////////////////////////////////////
{
  // quenching
  BirksSaturation(edeposit,aTrack);

  // Attenuation length
  AttenuationLength(edeposit,aTrack,preStepPoint);

  return;
}

//////////////////////////////////////////////////////////////////////////////
void SBsimMRDResponse::BirksSaturation(G4double* edeposit, G4Track* aTrack)
//////////////////////////////////////////////////////////////////////////////
{
  //const G4double CBIRKS = sbcard->Birks; // tentatively, same as SciBar scinti. - marcus: replaced by const so we don't need sbcard
  const G4double CBIRKS = 0.0208;

  G4double              kineticE = aTrack->GetKineticEnergy();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4Material*           material = aTrack->GetMaterial();

  if(particle->GetPDGCharge()==0) return;

  G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/cm);
  *edeposit /=  (1.+CBIRKS*dedx);

  return;
}

//////////////////////////////////////////////////////////////////////////////
void SBsimMRDResponse::AttenuationLength(G4double* edeposit, G4Track* aTrack, 
					 G4StepPoint* preStepPoint)
//////////////////////////////////////////////////////////////////////////////
{

  G4ParticleDefinition* particle = aTrack->GetDefinition();
  if(particle->GetPDGCharge()==0) return;

  
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector worldPosition = preStepPoint->GetPosition();
  G4ThreeVector localPosition = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint(worldPosition);

  // Hit positions in the view of each counter
  G4double hitposx = localPosition[0]*mm;
  G4double hitposy = localPosition[1]*mm;
  G4double hitposz = localPosition[2]*mm;

  MRDModule*   mrdmod = mrddb->GetMRDModuleInfo();
  G4double Horizontallength = mrdmod->HScintiSize[0]*cm;
  G4double   Verticallength = mrdmod->VScintiSize[1]*cm;
  G4double      Taperlength = mrdmod->TScintiSize[4]*cm;

  G4double AttLength_VScinti = mrdmod->Attlength[0]*cm;
  G4double AttLength_HScinti = mrdmod->Attlength[1]*cm;

  G4double Ldir, Lref;
  Ldir = 0; // Path length of photon directoly entering PMT 
  Lref = 0; // Path length of photon reflectively entering PMT

  if (theTouchable->GetVolume()->GetName() == "MRDHScinti_PV") {
    if (worldPosition[0] < 0) { // Horizontal Scintillator is placed without rotation
      Ldir = Horizontallength + hitposx + 2*Taperlength;
      Lref = Horizontallength - hitposx + 2*Horizontallength + 2*Taperlength;
    } else { 
      Ldir = Horizontallength - hitposx + 2*Taperlength;
      Lref = Horizontallength + hitposx + 2*Horizontallength + 2*Taperlength;
    }     
  } else if (theTouchable->GetVolume()->GetName() == "MRDVScinti_PV") {
    if (worldPosition[1] < 0) { // Vertical Scintillator is placed without rotation
      Ldir = Verticallength + hitposy + 2*Taperlength;
      Lref = Verticallength - hitposy + 2*Verticallength + 2*Taperlength;
    } else { 
      Ldir = Verticallength - hitposy + 2*Taperlength;
      Lref = Verticallength + hitposy + 2*Verticallength + 2*Taperlength;
    }     
  } else if (theTouchable->GetVolume()->GetName() == "MRDTHScinti_PV") {
    Ldir = Taperlength - hitposz;
    Lref = Taperlength + hitposz + 4*Horizontallength + 2*Taperlength;
  } else if (theTouchable->GetVolume()->GetName() == "MRDTVScinti_PV") {
    Ldir = Taperlength - hitposz;
    Lref = Taperlength + hitposz + 4*Verticallength + 2*Taperlength;
  }
  

  // Different attenuation length is applied for horizontal and vertical counter
  if (theTouchable->GetVolume()->GetName() == "MRDHScinti_PV" || 
      theTouchable->GetVolume()->GetName() == "MRDTHScinti_PV") {
    *edeposit *= (exp(-1*Ldir/AttLength_HScinti) + exp(-1*Lref/AttLength_HScinti));
  } else {
    *edeposit *= (exp(-1*Ldir/AttLength_VScinti) + exp(-1*Lref/AttLength_VScinti));
  }

  return;
}
