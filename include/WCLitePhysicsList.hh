
#ifndef WCLitePhysicsList_h
#define WCLitePhysicsList_h 1
//#include "G4VUserPhysicsList.hh"
//#include "WCLitePhysicsMessenger.hh"
#include "globals.hh"

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VModularPhysicsList.hh"
//#include "CompileTimeConstraints.hh"

class WCLitePhysicsList: public G4VModularPhysicsList
{
  public:
    WCLitePhysicsList();
    virtual ~WCLitePhysicsList();

  private:
//    WCLitePhysicsMessenger* PhysicsMessenger;
//  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

  protected:
    // Construct particle and physics process
//    void ConstructParticle();
//    void ConstructProcess();
    virtual void SetCuts();

  protected:
    // these methods Construct particles 

  protected:
    // these methods Construct physics processes and register them

};

#endif







