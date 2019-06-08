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
// $Id: WCLiteSteppingAction.hh,v0, 26/12/2015 $
// GEANT4 tag $Name: geant4-09-04-patch-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WCLiteSteppingAction_H
#define WCLiteSteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4OpBoundaryProcess.hh"

class WCLiteEventAction;
class WCLiteTrackingAction;

enum PMThitStatus { active=1, hitPMT=2, absorbed=4, boundaryAbsorbed=8, hitSphere=16, inactive=14};
/* PMThitStatus:
   active: still being tracked
   hitPMT: stopped by being detected in a PMT
   absorbed: stopped by being absorbed with G4OpAbsorbtion
   boundaryAbsorbed: stopped by being aborbed with G4OpAbsorbtion
   hitSphere: track hit the sphere at some point
   inactive: track is stopped for some reason
    -This is the sum of all stopped flags so can be used to remove stopped flags  
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WCLiteSteppingAction : public G4UserSteppingAction
{
  public:
    WCLiteSteppingAction();//WCLiteEventAction* eventAction);
    virtual ~WCLiteSteppingAction();
    virtual void UserSteppingAction(const G4Step*);

  private:   
   G4OpBoundaryProcessStatus fExpectedNextStatus;
};
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ToName(G4OpBoundaryProcessStatus boundaryStatus);
G4String ToName2(G4StepStatus stepStatus);
G4bool LocateParticle(G4ThreeVector inpos);
#endif

