
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
#include "WCLiteTrajectory.hh"
#include "WCLiteTrackingAction.hh"
#include "WCLiteTrackInformation.hh"
#include "WCLiteDetectorConstruction.hh"
//#include "RecorderBase.hh" 
#include <iostream>
#include <fstream>
#include "G4TrackingManager.hh"
#include "G4Trajectory.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "G4TrajectoryContainer.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4RunManager.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::WCLiteTrackingAction()
{

  callcount==0;

  // Tree Declarations

  HitTreeEnd = new TTree("hittreeend","hittreeend");

  HitTreeEnd->Branch("trackID0",&TrackID_S);
  HitTreeEnd->Branch("pdgid",&partcode);
  HitTreeEnd->Branch("parentID0",&parentID_S);
  HitTreeEnd->Branch("Et0",&TEStart);
  HitTreeEnd->Branch("E0",&KEStart);
  HitTreeEnd->Branch("x0",&xStart);
  HitTreeEnd->Branch("y0",&yStart);
  HitTreeEnd->Branch("z0",&zStart);
  HitTreeEnd->Branch("t0",&tStart);
  HitTreeEnd->Branch("processStart",&processStart);
  HitTreeEnd->Branch("trackIDf",&TrackID_E);
  HitTreeEnd->Branch("parentIDf",&parentID_E);
  HitTreeEnd->Branch("Etf",&TEEnd);
  HitTreeEnd->Branch("Ef",&KEEnd);
  HitTreeEnd->Branch("xf",&xEnd);
  HitTreeEnd->Branch("yf",&yEnd);
  HitTreeEnd->Branch("zf",&zEnd);
  HitTreeEnd->Branch("tf",&tEnd);
  HitTreeEnd->Branch("wavelength",&wavelength);
  HitTreeEnd->Branch("processEnd",&processEnd);
  HitTreeEnd->Branch("phhit",&phhit);
  HitTreeEnd->Branch("isScatPhot",&isScatPhot);


  // Variable Decalarations
  line=0;
  numhits=0;
  numphots=0;
  nummuplus=0;
  nummuminus=0;
  numeplus=0;
  numeminus=0;
  numgamma=0;
  back_count=0;
  Etot=0;
  iter=1;
  iter2=1;
  numtrajects=0;	//marcus

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void WCLiteTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);
  //Use custom trajectory class
  //  fpTrackingManager->SetTrajectory(new WCLiteTrajectory(aTrack));
  
  G4String theProcess;
  G4String theLocation;

  isScatPhot=0.0;

  if( aTrack->GetNextVolume() != 0 ) {
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
  G4int stepno = (aTrack->GetCurrentStepNumber());
  G4Step* aStep = (G4Step*) (aTrack->GetStep());
 
  theProcess+="UserLimit";

  processStart = -1;
 
  detectorStart = theLocation;
  
  G4String thisparticlename = aTrack->GetDefinition()->GetParticleName();
  if (thisparticlename == "neutronDISABLED"){
		G4cout << "Opening neutron track. Position is: " << (aTrack->GetPosition()) << ". ";
		
		//  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() ){}
		WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());
		if (info == 0){
			G4cout << "No user info. Marking as a primary neutron. TrackID is: " << (aTrack->GetTrackID()) << G4endl;
			WCLiteTrackInformation* infoNew = new WCLiteTrackInformation();
			infoNew->SetTrackID(aTrack->GetTrackID());
			WCLiteTrajectory* newTrajectory = new WCLiteTrajectory(aTrack);		// create a new trajectory
			numtrajects++;
			fpTrackingManager->SetTrajectory(newTrajectory);			// and use that
			infoNew->SetParentTrajectory(newTrajectory);
			G4Track* theTrack = (G4Track*)aTrack;					// the received track pointer is a const,
												// we can't modify it. Need to work around that.
			theTrack->SetUserInformation(infoNew);
		} else { 									// does our track have user information?
			G4cout << " Has user info. This TrackID is: " << (aTrack->GetTrackID()) << ", ";
			G4cout << "while stored trackID is: " << info->GetTrackID() << ". ";
			if (info->GetTrackID()==0){info->SetTrackID(aTrack->GetTrackID());}
			if (info->GetIsScatteredN() && info->GetIsPrimaryN()){G4cout << "oh dear. Says it's both primary and scattered." << G4endl;}
			if (info->GetIsScatteredN()){						// does that information say it's a scattered neutron?
				G4cout << "Info says it's a scattered neutron. Parent tack is: " << (info->GetOriginalTrackID());
				WCLiteTrajectory* parenttraj = info->GetParentTrajectory();	// retrieve the parent trajectory
				G4cout << ", parent trajectory is: " << parenttraj;									
				G4int numPoints = parenttraj->GetPointEntries();									
				G4cout << " with " << numPoints << " points in it so far. ";// << G4endl;
				WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(0));
			 	G4double pntx = firstpnt->GetPosition().x();
 				G4double pnty = firstpnt->GetPosition().y();
 				G4double pntz = firstpnt->GetPosition().z();
     				G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
     				WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(numPoints-1));
			 	pntx = lastpnt->GetPosition().x();
 				pnty = lastpnt->GetPosition().y();
 				pntz = lastpnt->GetPosition().z();
     				G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
				fpTrackingManager->SetTrajectory(parenttraj);				// and tell the tracking manager to use that
			} else if (info->GetIsPrimaryN()){
				WCLiteTrajectory* parenttraj = info->GetParentTrajectory();
				G4cout << "Info says it's a primary neutron. Using existing trajectory, " << parenttraj << " which has ";
				G4int numPoints = parenttraj->GetPointEntries();									
				G4cout << numPoints << " points in it so far. " << G4endl;
				WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(0));
			 	G4double pntx = firstpnt->GetPosition().x();
 				G4double pnty = firstpnt->GetPosition().y();
 				G4double pntz = firstpnt->GetPosition().z();
     				G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
     				WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(numPoints-1));
			 	pntx = lastpnt->GetPosition().x();
 				pnty = lastpnt->GetPosition().y();
 				pntz = lastpnt->GetPosition().z();
     				G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
				fpTrackingManager->SetTrajectory(parenttraj);		// and use that.
			}
		}
	} else {									// not a neutron
		WCLiteTrajectory* thisTrajectory = new WCLiteTrajectory(aTrack);	// create a new trajectory
		numtrajects++;
		fpTrackingManager->SetTrajectory(thisTrajectory);			// and use that.
	}

		fpTrackingManager->SetStoreTrajectory(true);
		//    }
		//  else 
		//    fpTrackingManager->SetStoreTrajectory(false);
    //if (numtrajects > 135500){G4cout << "number of trajectories so far: " << numtrajects << G4endl;}
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void WCLiteTrackingAction::PostUserTrackingAction(const G4Track* aTrack){ 

  callcount++;

  phhit=0;

  G4String theProcess;
  G4String theLocation;
  G4String sProcess;
  
  // G4Step* aStep = aStep->GetTrack();


  if( aTrack->GetNextVolume() != 0 ) { 
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
  G4int stepno = (aTrack->GetCurrentStepNumber());
  G4Step* aStep = (G4Step*) (aTrack->GetStep());

  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
    theProcess = ((G4VProcess*) (aStep->GetPostStepPoint()->GetProcessDefinedStep()))->GetProcessName();
  } else {
    theProcess+="UserLimit";
  }

  if(line%100000==0) G4cout<<"  Post tracking call number: "<<line<<" "<<detectorEnd<<G4endl;
  
  if(partname=="opticalphoton") numphots++;
  line++;  

  //  G4cout<<"post-end"<<G4endl;
  // Marcus addition:
  G4String newPartType = aTrack->GetDefinition()->GetParticleName();				// Get this particle type
  if (newPartType == "neutronDISABLED") {								// If it's a neutron...
  	G4cout << G4endl << "finished tracking neutron; final position " << (aTrack->GetPosition()); 
  	G4cout << ", poststeppoint process = " << theProcess << ". ";
	G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();		// Get list of secondaries
	WCLiteTrajectory* parenttraj = (WCLiteTrajectory*)fpTrackingManager->GimmeTrajectory();	// Get parent trajectory
	G4int numPoints = parenttraj->GetPointEntries();									
	G4cout << "Trajectory has " << numPoints << " points in it so far. " << G4endl;
	WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(0));
 	G4double pntx = firstpnt->GetPosition().x();
	G4double pnty = firstpnt->GetPosition().y();
	G4double pntz = firstpnt->GetPosition().z();
	G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
	WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(numPoints-1));
 	pntx = lastpnt->GetPosition().x();
	pnty = lastpnt->GetPosition().y();
	pntz = lastpnt->GetPosition().z();
	G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
	if(secondaries) { 								// If it has any
		size_t nSeco = secondaries->size();
	  if(nSeco>0) {									// If it has any...
	  	G4cout << "it had " << nSeco << " daughters; ";
	  	G4int numDaughterNeutrons = 0;						// count how many daughter neutrons
	  	G4Track* daughterTrack;							// and hold the track if only one
	    for(size_t i=0;i<nSeco;i++) { 						// scan through secondaries
	    	G4String daughterPartType = (*secondaries)[i]->GetDefinition()->GetParticleName();	
											// is it a neutron
	    	if (daughterPartType == "neutron"){					// if so
	    		numDaughterNeutrons++;						// increment count
	    		if(numDaughterNeutrons > 1){					// do we now have > 1?
	    			break;
	    			G4cout << "but more than one neutron.";// << G4endl;
	    		}														
	    		daughterTrack = (*secondaries)[i];				// if not store the track.
	    	}									// end of neutron processing
	    }										// end of loop over daughters
	    if (numDaughterNeutrons==1){						// if we retrieved one and only one daughter
	    	G4cout << "of which one was a neutron. Storing parent trajectory: ";
	        WCLiteTrackInformation* infoNew = new WCLiteTrackInformation(aTrack);	// create a new UserInformation object
	    	//infoNew->SetIsScatteredN(true);					// set the scattered neutron property to true
	    	infoNew->SetParentTrajectory(parenttraj);				// set the parent trajectory pointer
	    	G4cout << parenttraj << ". ";// << G4endl;
	    	daughterTrack->SetUserInformation(infoNew);				// and link it to the daughter.
	    } //else { G4cout << G4endl;}
	 }
	} //else { G4cout << G4endl;}
  }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::~WCLiteTrackingAction() {}
















