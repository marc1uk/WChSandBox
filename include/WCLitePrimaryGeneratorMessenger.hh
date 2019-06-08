#ifndef WCLitePrimaryGeneratorMessenger_h
#define WCLitePrimaryGeneratorMessenger_h 1

class WCLitePrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"

class WCLitePrimaryGeneratorMessenger: public G4UImessenger
{
 public:
  WCLitePrimaryGeneratorMessenger(WCLitePrimaryGeneratorAction* mpga);
  ~WCLitePrimaryGeneratorMessenger();
  
 public:
  void     SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);
  
 private:
  WCLitePrimaryGeneratorAction* myAction;
  
 private: //commands
  G4UIdirectory*      mydetDirectory;
  G4UIcmdWithAString* genCmd;
  G4UIcmdWithAString* fileNameCmd;
  
};

#endif


