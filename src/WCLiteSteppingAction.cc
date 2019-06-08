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
// $Id: WCLiteSteppingAction.cc,v0, 26/12/2015 $
// GEANT4 tag $Name: geant4-09-04-patch-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "WCLiteSteppingAction.hh"
//#include "MRDSD.hh"
#include "mrdPMTSD.hh"
#include "WCLiteEventAction.hh"
//#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
//#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "WCLiteTrackInformation.hh"
//#include "G4TrackStatus.hh"
//#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ios.hh"
//#include "G4RunManager.hh"
//#include "G4LogicalVolume.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteSteppingAction::WCLiteSteppingAction()//WCLiteEventAction* eventAction)// : G4UserSteppingAction()
{fExpectedNextStatus = Undefined;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteSteppingAction::~WCLiteSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool onscintboundary, onzboundary, onxboundary, onyboundary, frontface;
G4int scintlayer;
G4double xpos, ypos, zpos, xextent, yextent;
extern G4int numpanels;
extern G4double steelfullzlen;
extern G4double scintfullzlen;
extern G4double alufullzlen;	
extern G4double scintfullxlen;
extern G4double scintvfullylen;
extern G4double scinthfullylen;
G4double layerthickness;
std::vector<G4double> scintzedges;

void WCLiteSteppingAction::UserSteppingAction(const G4Step* theStep){
  G4Track* aTrack = theStep->GetTrack();
      
  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume(); 

  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();     
      
  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  //use 'static' variable - will retain it's value between calls ...
  static G4OpBoundaryProcess* boundary=NULL;
  // ... so that once boundary has a value, the following will be skipped:
  // Retrieve pointer to the 'optical boundary process' so that we may retrive it's status
  // (invoked, what type of action etc) at the end of each step
  if(!boundary){
    G4ProcessManager* pm = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	G4cout<<"optical boundary process: "<<boundary<<G4endl;
	break;
      }
    }
  }

  //if(fOneStepPrimaries&&thePrePV->GetName()=="scintillator"){
  //  aTrack->SetTrackStatus(fStopAndKill);
  //}
  
  if(!thePostPV){
    //out of world
    fExpectedNextStatus=Undefined;
    return;
  }

 	static G4bool doonce = true;
  	if(doonce){
		scintzedges.push_back(alufullzlen);				// first panel is offset by depth of one alu struct. 
		scintzedges.push_back(alufullzlen+scintfullzlen);		// rear face fo first panel
		G4cout<<"layer one at z="<<scintzedges.at(0)<<" and "<<scintzedges.at(1)<<G4endl;
		// all others are at even intervals thereafter.
		layerthickness = steelfullzlen + alufullzlen + scintfullzlen;
		for(G4int i=1;i<numpanels;i++){
			scintzedges.push_back(i*layerthickness);		// front face
			scintzedges.push_back(i*layerthickness+scintfullzlen);	// rear face
			G4cout<<"layer "<<i<<" at z="<<scintzedges.at(i*2)<<" and "<<scintzedges.at((i*2)+1)<<G4endl;
		}
		doonce=false;
	}
	
  G4ParticleDefinition* particleType = aTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){ 	//Optical photon only

      // just for informational purposes******************
/*  	if((thePrePV->GetName()!="waterTank")&&(thePrePV->GetName()!="Hall")){

			G4cout<<G4endl;
			G4cout<<"poststepproc : "<<thePostPoint->GetProcessDefinedStep()->GetProcessName()<<G4endl;
			G4cout<<"poststepstat : "<<ToName2(thePostPoint->GetStepStatus())<<G4endl;
			if(thePostPoint->GetStepStatus()){G4cout<<"boundarystat : "<<ToName(boundaryStatus=boundary->GetStatus())<<G4endl;}
			G4cout<<"pre position : ("<<thePrePoint->GetPosition().x()<<","<<thePrePoint->GetPosition().y()<<","<<thePrePoint->GetPosition().z()<<")"<<G4endl;
			G4cout<<"pre momentum : ("<<thePrePoint->GetMomentumDirection().x()<<","<<thePrePoint->GetMomentumDirection().y()<<","<<thePrePoint->GetMomentumDirection().z()<<")"<<G4endl;
			G4cout<<"post position : ("<<thePostPoint->GetPosition().x()<<","<<thePostPoint->GetPosition().y()<<","<<thePostPoint->GetPosition().z()<<")"<<G4endl;
			G4cout<<"post momentum : ("<<thePostPoint->GetMomentumDirection().x()<<","<<thePostPoint->GetMomentumDirection().y()<<","<<thePostPoint->GetMomentumDirection().z()<<")"<<G4endl;
			G4cout<<"from "<<thePrePoint->GetTouchableHandle()->GetVolume()->GetName()<<" ("<<thePrePoint->GetTouchableHandle()->GetCopyNumber()<<") to "<<thePostPoint->GetTouchableHandle()->GetVolume()->GetName()<<" ("<<thePostPoint->GetTouchableHandle()->GetCopyNumber()<<")"<<G4endl;
			G4cout<<G4endl;

	  	G4ThreeVector inpos = thePostPoint->GetPosition();
	  	G4bool inpaddle = LocateParticle(inpos);
	    if(onscintboundary==true){//&&onzboundary==false
				G4ThreeVector worldPosition = thePostPoint->GetPosition();
				G4ThreeVector localPosition = thePostPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
		  	if(onscintboundary){
		  		if(onxboundary){G4cout<<"on the x boundary ";}
		  		else if(onyboundary){G4cout<<"on the y boundary ";}
		  		else if(onzboundary){G4cout<<"on the z boundary ";}
		  		G4cout<<"of layer "<<scintlayer<<G4endl;
		  	} else {
		  		if(inpaddle){G4cout<<"in layer "<<scintlayer<<G4endl;}
		  		else {G4cout<<"not in the paddles.";}
		  	}
	  	}
  	}
  	// end informational purposes******************
*/  

   //WCLiteTrackInformation* info = (WCLiteTrackInformation*)(aTrack->GetUserInformation());
  
   // Kill photons entering lg_phys.  
   if(thePostPV->GetName()=="lg_phys"){
   		aTrack->SetTrackStatus(fStopAndKill);
   		G4cout<<"killing optical track entering LG "<<thePostPoint->GetTouchableHandle()->GetCopyNumber()<<G4endl;
   }
   if(thePostPoint->GetProcessDefinedStep()->GetProcessName() == "OpAbsorption"){
       //if(info){info->SetMRDdetected(absorbed);}
   }
    
   boundaryStatus=boundary->GetStatus();
   // double check to ensure particle at boundary (may not be valid for vers <Geant4.6.0-p1? omittable >?)
   if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
         	G4cout<<"expected status step too small, but boundary status isn't!"<<G4endl;
          /*G4ExceptionDescription ed; 		// exceptiondescription not in geant4.9.4
          ed << "SteppingAction::UserSteppingAction(): "
             << "No reallocation step after reflection!"
             << G4endl;
          G4Exception("SteppingAction::UserSteppingAction()", "Expl01", FatalException,ed,
          "Something is wrong with the surface normal or geometry"); */
          G4Exception("SteppingAction::UserSteppingAction()","Expl01",FatalException,"No reallocation step after reflection!");
          // G4Exception(const char* issure, const char* errorCode, G4ExceptionSeverity severity, const char* comments);
        }
      } 
      fExpectedNextStatus=Undefined;
            
      // perform actions based on boundary status... 
      switch(boundaryStatus){
	      case Absorption:
	      	  G4cout<<"Absorption of optical photon on boundary of "<<thePostPV->GetName()<<G4endl;
		  //info->SetMRDdetected(boundaryAbsorbed);
		  break;
	      case Detection: {
	          //G4cout<<"Detection of optical photon on boundary of "<<thePostPV->GetName()<<G4endl;
		  G4SDManager* SDman = G4SDManager::GetSDMpointer();
		  G4String sdName="MRDPMTSD";
		  mrdPMTSD* pmtSD = (mrdPMTSD*)SDman->FindSensitiveDetector(sdName);
		  if(pmtSD){
		    pmtSD->ProcessHits_PhotonHit(theStep, NULL);
		    //G4cout<<"called process photon hit!"<<G4endl;
		    //info->SetMRDdetected(hitPMT);
		  }
		  break;
	      }
	      case FresnelReflection:
	      case TotalInternalReflection:
	      case SpikeReflection:
	      	 // trackInformation->IncReflections();
	      	 break;
	      case LambertianReflection:
	      case LobeReflection:
	      case BackScattering: {
		 // trackInformation->IncReflections();
		 fExpectedNextStatus=StepTooSmall;
		 //G4cout<<"boundary status backscattering - ignore the next step (too small)"<<G4endl;
		 break;
	      }
	      default:
		break;
      }
   } // post-step status not on boundary
  } // not an optical photon
  
}

G4String ToName(G4OpBoundaryProcessStatus boundaryStatus){
G4String processnames[14]={"Undefined","FresnelRefraction","FresnelReflection","TotalInternalReflection","LambertianReflection","LobeReflection","SpikeReflection","BackScattering","Absorption","Detection","NotAtBoundary","SameMaterial","StepTooSmall","NoRINDEX"};
return processnames[boundaryStatus];}

G4String ToName2(G4StepStatus stepStatus){
G4String processnames[8]={"fWorldBoundary","fGeomBoundary","fAtRestDoItProc","fAlongStepDoItProc","fPostStepDoItProc","fUserDefinedLimit","fExclusivelyForcedProc","fUndefined"};
return processnames[stepStatus];}


G4bool LocateParticle(G4ThreeVector inpos){
	zpos=inpos.z();
	zpos=zpos-2.0*m;	// remove offset to put MRD z front face at z=0;
	xpos=inpos.x();
	ypos=inpos.y();
	onscintboundary=false;
	onzboundary=false;
	onxboundary=false;
	onyboundary=false;
	scintlayer=-1;
	xextent=-1;
	yextent=-1;
	// scan through front/back faces of scint paddles, see if we're on one of these.
	for(G4int i=0;i<(numpanels*2);i++){
		if(zpos==scintzedges.at(i)){
			onscintboundary=true;
			onzboundary=true;	// exactly equal - we're on a front/back face
			scintlayer=floor(i/2.);	// floor to find which panel we're in
			frontface=(((i/2.)-floor(i/2.))==0);	// mod to find whether we're on front or back face.
			if(((scintlayer/2.)-floor(scintlayer/2.))==0){	// layer mod to find if we're in horiz or vert. oriented panel.
				xextent = scinthfullylen;	// horizontal orientation first.
				yextent = scintfullxlen*15;
			} else {
				xextent = scintfullxlen*15;
				yextent = scintvfullylen;
			}
			return true;	// we've found what boundary we're at and details we're interested in.
		}
	}
	// got this far - aren't on a front/back face. find which panel we're IN, to know top/bottom, then check those boundaries
	for(G4int i=0;i<numpanels;i+2){
		if(zpos>scintzedges.at(i)&&zpos<scintzedges.at(i+1)){
			scintlayer=i/2;		// 
			if(((scintlayer/2.)-floor(scintlayer/2.))==0){	// layer mod to find if we're in horiz or vert. oriented panel.
				xextent = scinthfullylen;	// horizontal orientation first.
				yextent = scintfullxlen*15;
			} else {
				xextent = scintfullxlen*15;
				yextent = scintvfullylen;
			}
			if(abs(xpos)==xextent){
				onscintboundary=true;
				onxboundary=true;
			} else if(abs(ypos)==yextent){
				onscintboundary=true;
				onyboundary=true;
			}
			return true;	// found what panel we're in and set the corresponding stuff.
		}
	}
	// apparently not in a scintillator panel. 
	return false;
}	

