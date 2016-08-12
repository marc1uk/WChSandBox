
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
#include "G4TrackStatus.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::WCLiteTrackingAction()
{

//  callcount==0;

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

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void WCLiteTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);
  //Use custom trajectory class
  //  fpTrackingManager->SetTrajectory(new WCLiteTrajectory(aTrack));
  
  // if this particle was a secondary generated via ncapture, record it. this will be propagated to secondaries in post-action
  G4String theProcess="primary";
  if(aTrack->GetCreatorProcess() != 0){
  	theProcess = aTrack->GetCreatorProcess()->GetProcessName();
  }
  if(theProcess=="nCapture"){
  	WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());
  	if(info==0){
  		info = new WCLiteTrackInformation();
  		G4Track* theTrack = (G4Track*)aTrack;	// the received track pointer is a const, we can't modify it. Need to work around that.
			theTrack->SetUserInformation(info);
  	}
  	info->SetGdOriginalTrackID(aTrack->GetParentID());
  }
  
  G4String theLocation;

  isScatPhot=0.0;

  if( aTrack->GetNextVolume() != 0 ) {
    theLocation = aTrack->GetVolume()->GetName();
  } else {
    theLocation+="DetectorHit";
  }
  
//  G4int stepno = (aTrack->GetCurrentStepNumber());
//  G4Step* aStep = (G4Step*) (aTrack->GetStep());

  processStart = -1;
 
  detectorStart = theLocation;
  //  if ( aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() ){
  		
  G4String thisparticlename = aTrack->GetDefinition()->GetParticleName();
  if (thisparticlename == "neutron"){	// change neutron to something random to circumvent this and give new trajectories to everything.
		//G4cout << "Opening neutron track. Position is: " << (aTrack->GetPosition()) << ". ";
		WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());
		if (info == 0){
			//G4cout << "No user info. Marking as a primary neutron. TrackID is: " << (aTrack->GetTrackID()) << G4endl;
			WCLiteTrackInformation* infoNew = new WCLiteTrackInformation();
			infoNew->SetTrackID(aTrack->GetTrackID());
			WCLiteTrajectory* newTrajectory = new WCLiteTrajectory(aTrack);			// create a new trajectory
			fpTrackingManager->SetTrajectory(newTrajectory);										// and use that
			//infoNew->SetParentTrajectory(newTrajectory);
			G4Track* theTrack = (G4Track*)aTrack;	// the received track pointer is a const, we can't modify it. Need to work around that.
			theTrack->SetUserInformation(infoNew);
		}	else { 																												      // does our track have user information?
			//G4cout << " Has user info. This TrackID is: " << (aTrack->GetTrackID()) << ", ";
			//G4cout << "while stored trackID is: " << info->GetTrackID() << ". ";
			if (info->GetTrackID()==0){info->SetTrackID(aTrack->GetTrackID());}
			/* JUST PRINTING (previously would select whether to use parent trajectory; no longer done due to memory issues.)
			if (info->GetIsScatteredN() && info->GetIsPrimaryN()){G4cout << "oh dear. Says it's both primary and scattered." << G4endl;}
			if (info->GetIsScatteredN()){																				// does that information say it's a scattered neutron?
				G4cout << "Info says it's a scattered neutron."; 
			} else if (info->GetIsPrimaryN()){
				G4cout << "Info says it's a primary neutron. ";
			}
			WCLiteTrajectory* parenttraj = info->GetParentTrajectory();				// retrieve the parent trajectory
			G4cout << "Parent tack is: " << (info->GetOriginalTrackID()) << ". Parent trajectory is: " << parenttraj;						
			G4int numPoints = parenttraj->GetPointEntries();									
			G4cout << " with " << numPoints << " points in it so far. ";// << G4endl;
			WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(0));
			if(firstpnt!=0){
			 	G4double pntx = firstpnt->GetPosition().x();
				G4double pnty = firstpnt->GetPosition().y();
				G4double pntz = firstpnt->GetPosition().z();
		 		G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
		 	}
   		WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(parenttraj->GetPoint(numPoints-1));
   		if(lastpnt!=0){
			 	pntx = lastpnt->GetPosition().x();
				pnty = lastpnt->GetPosition().y();
				pntz = lastpnt->GetPosition().z();
		 		G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;	
		 	}
			// previously neutrons would have retrieved parent trajectory to pass to tracking manager.     		
			*/
			WCLiteTrajectory* newTrajectory = new WCLiteTrajectory(aTrack);
		  fpTrackingManager->SetTrajectory(newTrajectory);	
		}
	} else {		// not a neutron
		WCLiteTrajectory* newTrajectory = new WCLiteTrajectory(aTrack);	    // create a new trajectory
		fpTrackingManager->SetTrajectory(newTrajectory);									    // and use that.
	}
			
	/*if(aTrack->GetTrackID()==44953){
		// problematic track
		WCLiteTrackInformation* infoNew = new WCLiteTrackInformation();
		infoNew->SetTrackID(aTrack->GetTrackID());
		G4Track* theTrack = (G4Track*)aTrack;	// the received track pointer is a const, we can't modify it. Need to work around that.
		theTrack->SetUserInformation(infoNew);
	}*/

		fpTrackingManager->SetStoreTrajectory(true);
		//    }
		//  else 
		//    fpTrackingManager->SetStoreTrajectory(false);
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
  
//  G4int stepno = (aTrack->GetCurrentStepNumber());
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
  // all this now does is, for neutrons with one (or more) daughter neutrons, set the first daughter neutron as the 
  // same particle, by recording the parent track ID in user info under sameasparenttrackid. No trajectory use will be changed.
  G4String newPartType = aTrack->GetDefinition()->GetParticleName();			// Get this particle type
  if (newPartType == "neutron") {							// If it's a neutron... CHANGE TO CIRCUMVENT
  	/* JUST PRINTING
  	WCLiteTrajectory* thistraj = (WCLiteTrajectory*)fpTrackingManager->GimmeTrajectory();	// Get used trajectory
  	G4cout << G4endl << "finished tracking neutron; final position " << (aTrack->GetPosition()); 
  	G4cout << ", poststeppoint process = " << theProcess << ". ";
		G4int numPoints = thistraj->GetPointEntries();									
		G4cout << "Trajectory Manager trajectory has " << numPoints << " points in it so far. " << G4endl;
		WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(thistraj->GetPoint(0));
	 	G4double pntx = firstpnt->GetPosition().x();
		G4double pnty = firstpnt->GetPosition().y();
		G4double pntz = firstpnt->GetPosition().z();
		G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
		WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(thistraj->GetPoint(numPoints-1));
	 	pntx = lastpnt->GetPosition().x();
		pnty = lastpnt->GetPosition().y();
		pntz = lastpnt->GetPosition().z();
		G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
		*/
		WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());	// Get track info
		/* JUST PRINTING
		WCLiteTrajectory* storedtraj = info->GetParentTrajectory();	// Get stored trajectory
		numPoints = storedtraj->GetPointEntries();									
		G4cout << "Stored trajectory has " << numPoints << " points in it so far. " << G4endl;
		firstpnt = (WCLiteTrajectoryPoint*)(storedtraj->GetPoint(0));
		pntx = firstpnt->GetPosition().x();
		pnty = firstpnt->GetPosition().y();
		pntz = firstpnt->GetPosition().z();
		G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
		lastpnt = (WCLiteTrajectoryPoint*)(storedtraj->GetPoint(numPoints-1));
		pntx = lastpnt->GetPosition().x();
		pnty = lastpnt->GetPosition().y();
		pntz = lastpnt->GetPosition().z();
		G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
		*/ 
		G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();		// Get list of secondaries
		if(secondaries) { 																										// If it has any
			size_t nSeco = secondaries->size();
			if(info){info->AddNumSecs(nSeco);}																	// if existing info, add these secondaries
		  if(nSeco>0) {																												// If it has any...
		  	//G4cout << "it had " << nSeco << " daughters; ";
		  	G4int numDaughterNeutrons = 0;																	// count how many daughter neutrons
		  	G4Track* daughterTrack=0;																					// and hold the track if only one
		    for(size_t i=0;i<nSeco;i++) { 																		// scan through secondaries
		    	G4String daughterPartType = (*secondaries)[i]->GetDefinition()->GetParticleName();	// is it a neutron
		    	if (daughterPartType == "neutron"){															// if so
		    		numDaughterNeutrons++;																			// increment count
		    		//if(numDaughterNeutrons > 1){																// do we now have > 1?
		    		//	G4cout << "and more than one neutron.";// << G4endl;
		    		//}														
		    		daughterTrack = (*secondaries)[i];														// if not store the track.
		    		break;
		    	}																																// end of neutron processing
		    }																																	// end of loop over daughters
		    if (daughterTrack!=0 && (G4int)(aTrack->GetTrackStatus())==2){	// if we retrieved a daughter neutron
		    	//G4cout << "of which at least one was a neutron. Storing parent trajectory: ";
		      WCLiteTrackInformation* infoNew = new WCLiteTrackInformation(aTrack);	// create a new UserInformation object
		    	//infoNew->SetParentTrajectory(thistraj);												// set the parent trajectory pointer
		    	infoNew->SetSameAsParentTrackID(aTrack->GetTrackID());					// declare this daughter the same as it's parent.
		    	//G4cout << thistraj << ". ";// << G4endl;
		    	daughterTrack->SetUserInformation(infoNew);											// and link it to the daughter.
		  	} //else { G4cout << G4endl;}
		 	}
		} //else { G4cout << G4endl;}
	} //else {fpTrackingManager->SetStoreTrajectory(false);}		//marcus: try to stem memory leak
	
	
// for all particles with GdOriginalTrackID!=0, copy this information to all secondaries
// -------------------------------------------------------------------------------------
	WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());
	if (info!=0){		// primary particles that have not passed through the NCV may not have user info
		G4int gdorig = info->GetGdOriginalTrackID();															// did this or an ancestor result from n-capture?
		if(gdorig>0) {														// this track or one of it's ancestors was a secondary of n-capture
			G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();			// Get list of secondaries
			if(secondaries){	 																											// If it has any
				size_t nSeco = secondaries->size();																		// how many?
				if(nSeco>0) {																													// if really has any
					info->AddNumSecs(nSeco);																						// add number of secondaries to accumulator for this particle (track)
					for(size_t i=0;i<nSeco;i++){																				// for each secondary
						G4Track* secondarytrack = (*secondaries)[i];											// get it's track
						WCLiteTrackInformation* infoSeco;
						if(secondarytrack->GetUserInformation()){													// check for user info (created above for neutrons)
							infoSeco = (WCLiteTrackInformation*)(secondarytrack->GetUserInformation());	// use if exists
						} else {
							infoSeco = new WCLiteTrackInformation();												// create a new UserInformation object
							secondarytrack->SetUserInformation(infoSeco);										// and link it to the daughter.
							infoSeco->SetIsPrimaryN(false);
						}
						infoSeco->SetGdOriginalTrackID(gdorig);													// set originator
					}
				}
			}
		} // else: neither this nor any ancestor tracks is a product of n-capture
	}
	// we're at the end of processing a Track that had some useful information.
	// If we don't want to lose it, we need to copy to the trajectory.
	if(info!=0){
		WCLiteTrajectory* thistraj = (WCLiteTrajectory*)fpTrackingManager->GimmeTrajectory();
		if((G4int)(aTrack->GetTrackStatus())==2){	// status 2 = fStopAndKill
			//G4cout << " end tracking with track status: " << aTrack->GetTrackStatus() << G4endl;
			//G4cout << "copying ";
			thistraj->CopyInfo(info);
			//thistraj->Print();
		}
	}
 //if(aTrack->GetTrackStatus()==2 && aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="Transportation"){
 // G4cout<<"particle stopped and killed with transportation process at: "<<aTrack->GetPosition().z()<<G4endl;}
 /*if(aTrack->GetTrackID()==44953){
 	//problematic track
 	WCLiteTrajectory* thistraj = (WCLiteTrajectory*)fpTrackingManager->GimmeTrajectory();
 	G4int numPoints = thistraj->GetPointEntries();
 	G4cout << "track 44953 has: " << numPoints << " points." << G4endl;
 	WCLiteTrajectoryPoint* firstpnt = (WCLiteTrajectoryPoint*)(thistraj->GetPoint(0));
 	G4double pntx = firstpnt->GetPosition().x();
	G4double pnty = firstpnt->GetPosition().y();
	G4double pntz = firstpnt->GetPosition().z();
	G4cout << "firstpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")";
	WCLiteTrajectoryPoint* lastpnt = (WCLiteTrajectoryPoint*)(thistraj->GetPoint(numPoints-1));
 	pntx = lastpnt->GetPosition().x();
	pnty = lastpnt->GetPosition().y();
	pntz = lastpnt->GetPosition().z();
	G4cout << ", lastpoint: (" << pntx << ", " <<  pnty << ", " << pntz << ")" << G4endl;
	G4cout << "stopped tracking with status: " << (aTrack->GetTrackStatus()) << G4endl;
	G4cout << "tracking used trajectory: " << thistraj << G4endl;
}*/
/*
track statuses:
0: fAlive: Continue the tracking.
1: fStopButAlive: The track has come to zero kinetic energy, but still AtRest process to occur.
2: fStopAndKill: The track has lost its identity because it has decayed, interacted or gone beyond the world boundary. Secondaries will be pushed to the stack.
3: fKillTrackAndSecondaries: Kill the current track and also associated secondaries.
4: fSuspend: Suspend processing of the current track and push it and its secondaries to the stack.
5: fPostponeToNextEvent: Postpone processing of the current track to the next event. Secondaries are still being processed within the current event.
*/
		
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
WCLiteTrackingAction::~WCLiteTrackingAction() {}













