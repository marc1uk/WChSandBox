#include "WCLitePhysicsMessenger.hh"
#include "WCLitePhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"

WCLitePhysicsMessenger::WCLitePhysicsMessenger(WCLitePhysicsList* WCLitePhys)
  :WCLitePhysics(WCLitePhys)
{
  
  WCLiteDir = new G4UIdirectory("/WCLite/physics/secondaries/");
  WCLiteDir->SetGuidance("Commands to change secondary interaction model for protons");

  hadmodelCmd = new G4UIcmdWithAString("/WCLite/physics/secondaries/model",this);
  hadmodelCmd->SetGuidance("Available options: GHEISHA BERTINI BINARY");
  hadmodelCmd->SetGuidance("Description:");
  hadmodelCmd->SetGuidance("GHEISHA = standard, fast G4 hadronic interaction model");
  hadmodelCmd->SetGuidance("BERTINI = Bertini cascade model");
  hadmodelCmd->SetGuidance("BINARY  = Binary cascade model (2KM default)");
  hadmodelCmd->SetParameterName("secondaries", true, false);
  hadmodelCmd->SetDefaultValue("BINARY");
  hadmodelCmd->SetCandidates("GHEISHA BERTINI BINARY");

}

WCLitePhysicsMessenger::~WCLitePhysicsMessenger()
{
  delete hadmodelCmd;
  delete WCLiteDir;
}

void WCLitePhysicsMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == hadmodelCmd)
    WCLitePhysics->SetSecondaryHad(newValue);

}
