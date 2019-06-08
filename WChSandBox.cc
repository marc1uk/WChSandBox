//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: exampleN06.cc,v 1.14 2006/06/29 17:53:52 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// mail:        gum@triumf.ca
//     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4ios.hh"

#include "WCLiteDetectorConstruction.hh"
#include "WCLitePhysicsList.hh"
#include "WCLitePrimaryGeneratorAction.hh"
#include "WCLiteRunAction.hh"
#include "WCLiteEventAction.hh"
#include "WCLiteStackingAction.hh"
#include "WCLiteSteppingVerbose.hh"
#include "WCLiteTrackingAction.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Seed the random number generator manually
  //
  G4long myseed = 345354;
  CLHEP::HepRandom::setTheSeed(myseed);
  
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new WCLiteSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  //
  WCLiteDetectorConstruction* detector = new WCLiteDetectorConstruction;
  runManager-> SetUserInitialization(detector);

  //
  G4VUserPhysicsList* physics = new WCLitePhysicsList;
  runManager-> SetUserInitialization(physics);
  
#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // UserAction classes
  //
  G4UserRunAction* run_action = new WCLiteRunAction;
  runManager->SetUserAction(run_action);

  // Set user action classes
  WCLitePrimaryGeneratorAction* myGeneratorAction = new WCLitePrimaryGeneratorAction(detector);
  runManager->SetUserAction(myGeneratorAction);
  
  //
  G4UserEventAction* event_action = new WCLiteEventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserTrackingAction* tracking_action = new WCLiteTrackingAction;
  runManager->SetUserAction(tracking_action);
  //
  G4UserStackingAction* stacking_action = new WCLiteStackingAction;
  runManager->SetUserAction(stacking_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
    
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
   
  if (argc==1)   // Define UI session for interactive mode
    {
      G4UIsession* session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    
      UI->ApplyCommand("/control/execute vis.mac"); 
      UI->ApplyCommand("/mygen/vecfile fluxesandtables/numu_center.txt");
      session->SessionStart();
      delete session;
   }
   
  else         // Batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
   }
   
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
