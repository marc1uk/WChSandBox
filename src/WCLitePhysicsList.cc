
#include "WCLitePhysicsList.hh"
//#include "WCLitePhysicsMessenger.hh"

#include <iomanip>   

#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"

#include <G4UnitsTable.hh>
#include "G4StepLimiterPhysics.hh"

//WCLitePhysicsList::WCLitePhysicsList():  G4VUserPhysicsList(), PhysicsMessenger(0)
WCLitePhysicsList::WCLitePhysicsList(): G4VModularPhysicsList() //, PhysicsMessenger(0)
{
  G4DataQuestionaire it(photon, neutron);
  G4cout << "<<< Geant4 Physics List simulation engine: FTFP_BERT_HP 2.0"<<G4endl;
  G4cout <<G4endl;
  defaultCutValue = 0.7*mm;
  SetVerboseLevel(1);
  G4int ver=1;

//PhysicsMessenger = new WCLitePhysicsMessenger(this);
 
  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics(ver) );	//this->

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays 
  RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );

   // Hadron Physics
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(ver) );

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver) );
  
  // Step limiter physics?? Not part of FTFP_BERT_HP but often included externally in main.cc when this physics list is used
  RegisterPhysics( new G4StepLimiterPhysics(ver) );
  
  // Optical Physics
  G4OpticalPhysics* opticalPhysics =  new G4OpticalPhysics(0);
  RegisterPhysics(opticalPhysics);	// don't be verbose - it spouts stuff all over the screen
  
  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
}

WCLitePhysicsList::~WCLitePhysicsList()
{
//  delete PhysicsMessenger;
//  PhysicsMessenger = 0;
}


//----set cut values----

void WCLitePhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "WCLitePhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  /*
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(0.0, "proton");   
  */
  
  // sets the default cut value for all particle types
  // (method inherited from G4VUserPhysicsList)
  //this->SetCutsWithDefault();
  G4VUserPhysicsList::SetCuts();

  if (verboseLevel>0) DumpCutValuesTable();
}

