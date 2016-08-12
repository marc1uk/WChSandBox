
#ifndef WCLitePhysicsList_h
#define WCLitePhysicsList_h 1
#include "G4VUserPhysicsList.hh"
#include "WCLitePhysicsMessenger.hh"
#include "globals.hh"

class WCLitePhysicsList: public G4VUserPhysicsList
{
  public:
    WCLitePhysicsList();
    ~WCLitePhysicsList();

  //=================================
  // Added by JLR 2005-07-05
  //=================================
  // This method sets the model for 
  // hadronic secondary interactions
  void SetSecondaryHad(G4String hadval);

  private:
    WCLitePhysicsMessenger* PhysicsMessenger;

    G4String SecondaryHadModel;

    G4bool gheishahad;
    G4bool bertinihad;
    G4bool binaryhad;

  protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
    // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();
    void ConstructOp();
    void ConstructHad();
    void ConstructlArStepLimiter();


};

#endif







