#ifndef WCLitePhysicsMessenger_h
#define WCLitePhysicsMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class WCLitePhysicsList;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;

class WCLitePhysicsMessenger: public G4UImessenger
{
public:
  WCLitePhysicsMessenger(WCLitePhysicsList*);
  ~WCLitePhysicsMessenger();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:
  WCLitePhysicsList* WCLitePhysics;

  G4UIdirectory*      WCLiteDir;
  G4UIcmdWithAString* hadmodelCmd;

};

#endif
