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
// $Id: WCLiteDetectorConstruction.cc,v 1.15 2006/06/29 17:54:17 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "WCLiteDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "MRDDetectorConstruction.hh"   
#include "MRDSD.hh"
#include "G4SDManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpBoundaryProcess.hh" 
#include "mrdPMTSD.hh"

// **********SciBooNE integration
#include "SBsimMRDSD.hh"
#include "SBsimMRDDB.hh"
// **********/SciBooNE integration
void ConstructMRD(G4LogicalVolume* expHall_log, G4double mrdZoffset);
void ConstructNCV(G4LogicalVolume* waterTank_log);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteDetectorConstruction::WCLiteDetectorConstruction() // skipping arg from SciBooNE integration
{
  expHall_x = 50*m;
  expHall_y = expHall_z = 500*m;
  //  expHall_x = 10*m;
  //  expHall_y = expHall_z = 10*m;
  //  tank_x    = tank_y    = tank_z    = 5*m;
  //  bubble_x  = bubble_y  = bubble_z  = 0.5*m;
  
  // **************SciBooNE integration
  // InputCard = aInputCard;		skipping input card from SciBooNE integration
  // *************/SciBooNE integration
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteDetectorConstruction::~WCLiteDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WCLiteDetectorConstruction::Construct()
{

  ConstructMaterials();

//
//	------------- Volumes --------------

// The experimental Hall
//

  //G4Box* expHall_box = new G4Box("Hall",expHall_x,expHall_y,expHall_z);
  G4Box* expHall_box = new G4Box("Hall",90*m,90*m,25*m);		//just..nicer for now.
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,G4Material::GetMaterial("Air"),"Hall",0,0,0);
  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"Hall",0,false,0);

/*
  G4Tubs* expHall_box
    = new G4Tubs("waterTank", 
                 0.*cm, 
                 31.11*m, 
                 39.98*m, 
                 0.*deg, 
                 360.*deg); 
*/

  G4bool isWater = true;
  G4bool WCAddGd = false;
  G4bool isNCV = true; // Construct the Neutron Capture Volume

  //Decide if adding Gd
  G4String Scintillator = "Scintillator";
  G4String water = "Water";

  if (!isWater)
    { water = "Scintillator"; }

  if (WCAddGd && isWater)
    {water = "Doped Water";}
  
  if (WCAddGd && !isWater)
    {water = "Doped Scintillator";}

  
  //G4Sphere* expHall_box 
  //= new G4Sphere("waterTank",0.*m,10.*m,0.*deg,360.*deg,0.*deg,180.*deg);
  	
  tankouterRadius= 1.98*m;
  tankhz = 3.0*m;
  G4Tubs* waterTank_box = new G4Tubs("waterTank",0.*m,tankouterRadius,tankhz,0.*deg,360.*deg);	//prev dims: 5.5m radius, 11m height.
  //G4Tubs* aTube = new G4Tubs("Name", innerRadius, outerRadius, hz, startAngle, spanningAngle);
  
  // G4Box* expHall_box = new G4Box("waterTank",1.5*m,1.5*m,1.5*m); ANNIE
  // G4Box* expHall_box = new G4Box("waterTank",10.0*m,20.0*m,20.0*m);

  G4LogicalVolume* waterTank_log
    = new G4LogicalVolume(waterTank_box,G4Material::GetMaterial(water),"waterTank",0,0,0);
    
  G4RotationMatrix* rotm = new G4RotationMatrix();					// declare a rotation matrix
  rotm->rotateX(90*deg); 							// rotate Y
  
  G4VPhysicalVolume* waterTank_phys
    = new G4PVPlacement(rotm,G4ThreeVector(),waterTank_log,"waterTank",expHall_log,false,0);

 
  /* 
// The Water Tank
//	
  G4Box* waterTank_box = new G4Box("Tank",tank_x,tank_y,tank_z);

  G4LogicalVolume* waterTank_log
    = new G4LogicalVolume(waterTank_box,Water,"Tank",0,0,0);

  G4VPhysicalVolume* waterTank_phys
    = new G4PVPlacement(0,G4ThreeVector(),waterTank_log,"Tank",
                        expHall_log,false,0);

 
// The Air Bubble
//   
  G4Box* bubbleAir_box = new G4Box("Bubble",bubble_x,bubble_y,bubble_z);

  G4LogicalVolume* bubbleAir_log
    = new G4LogicalVolume(bubbleAir_box,Air,"Bubble",0,0,0);

//G4VPhysicalVolume* bubbleAir_phys =
      new G4PVPlacement(0,G4ThreeVector(0,2.5*m,0),bubbleAir_log,"Bubble",
                        waterTank_log,false,0);
  */

//	------------- Surfaces --------------
//
// Water Tank
//
  G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
  OpWaterSurface->SetType(dielectric_dielectric);
  OpWaterSurface->SetFinish(ground);
  OpWaterSurface->SetModel(unified);

  //  G4LogicalBorderSurface* WaterSurface = 
  //                               new G4LogicalBorderSurface("WaterSurface",
  //                               expHall_phys,expHall_phys,OpWaterSurface);

  //  if(WaterSurface->GetVolume1() == waterTank_phys) G4cout << "Equal" << G4endl;
  //  if(WaterSurface->GetVolume2() == expHall_phys  ) G4cout << "Equal" << G4endl;

  /*
// Air Bubble
//
  G4OpticalSurface* OpAirSurface = new G4OpticalSurface("AirSurface");
  OpAirSurface->SetType(dielectric_dielectric);
  OpAirSurface->SetFinish(polished);
  OpAirSurface->SetModel(glisur);
 
  G4LogicalSkinSurface* AirSurface = 
	  new G4LogicalSkinSurface("AirSurface", bubbleAir_log, OpAirSurface);

  if(AirSurface->GetLogicalVolume() == bubbleAir_log) G4cout << "Equal" << G4endl;
  ((G4OpticalSurface*)
  (AirSurface->GetSurface(bubbleAir_log)->GetSurfaceProperty()))->DumpInfo();

//
// Generate & Add Material Properties Table attached to the optical surfaces
//
  const G4int num = 2;
  G4double Ephoton[num] = {2.038*eV, 4.144*eV};

  //OpticalWaterSurface 
  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);

  //OpticalAirSurface
  G4double Reflectivity[num] = {0.3, 0.5};
  G4double Efficiency[num]   = {0.8, 1.0};

  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();

  myST2->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
  myST2->AddProperty("EFFICIENCY",   Ephoton, Efficiency,   num);

  OpAirSurface->SetMaterialPropertiesTable(myST2);
  */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  if (isNCV){
    G4ThreeVector NCVposition;
    NCVposition.setX(0.*m); // Position of the NCV
    NCVposition.setY(0.*m);
    NCVposition.setZ(0.*m);
    ConstructNCV(waterTank_log);
    G4cout << "************ Neutron Capture Volume will be included in the simulation ! ************\n";
  } else{
    G4cout << "************ Neutron Capture Volume will not be included in the simulation ! ************\n";
  }
  
  ConstructMRD(expHall_log, tankouterRadius);			// marcus' MRD construction
  expHall_log->SetVisAttributes (G4VisAttributes::Invisible);	// set hall volume invisible
  
  // ******SciBooNE integration
  G4PVPlacement* expHall= new G4PVPlacement(0, G4ThreeVector(), "EXP_HALL_PV", expHall_log, 0, false, 0);
  //DefineMRD(expHall);	
  // *****/SciBooNE integration
  
  //always return the physical World
  return expHall_phys;

}


void ConstructNCV(G4LogicalVolume* waterTank_log){
  //===============================================================================================================
  //NEUTRON CAPTURE VOLUME DEFINITION
  //===============================================================================================================
  // NCV is defined as two cylinders, one filled with liquid, the other with acrylic. Two 'caps' are added to the top of the NCV as well. 
  // The metal structure will been added as well (from pictures) since it will induce n-Fe captures
  //===============================================================================================================

  // Dimensions 
  G4double NCVliquid_radius = 25.*cm;
  G4double NCVliquid_height = 50.*cm;
  G4double NCVvessel_thickness = 1.*cm;
  G4double NCVvesselcap_thickness = 3.*cm;
  
  // NCV acrylic vessel
  // Cylinder
  G4VSolid* NCVvessel_tub
    = new G4Tubs("NCVvessel_tub",
		    0.,
		    NCVliquid_radius + NCVvessel_thickness, // r1, r2
		    NCVliquid_height/2., //Half height of the cylinder
		    0.0, 2.0*M_PI         // phi0, delta_phi
		    );

  G4LogicalVolume* NCVvessel_log
    = new G4LogicalVolume(NCVvessel_tub,
			  G4Material::GetMaterial("Acrylic"),
			  "NCVvessel_log",
			  0,
			  0,
			  0);
  G4VisAttributes* NCVvessel_vis
    = new G4VisAttributes(G4Color(0.1,0.,1.0,0));
  NCVvessel_log -> SetVisAttributes(NCVvessel_vis);


  G4VPhysicalVolume* NCVvessel_phys
  = new G4PVPlacement(0,  // no rotation
			G4ThreeVector(0.,0.,0.), // shifted in the water tank
			NCVvessel_log,
			"NCVvessel_phys",
			waterTank_log,          // mother
			false,
			0);
  
  // Caps
  G4VSolid* NCVvesselcap_box
    = new G4Box("NCVvesselcap_box",
		    NCVliquid_radius, // x
		    NCVliquid_radius, // y
		    NCVvesselcap_thickness/2. // thickness
		    );
    G4LogicalVolume* NCVvesselcap_log
    = new G4LogicalVolume(NCVvesselcap_box,
			  G4Material::GetMaterial("Acrylic"),
			  "NCVvesselcap_log",
			  0,
			  0,
			  0);
  G4VisAttributes* NCVvesselcap_vis
    = new G4VisAttributes(G4Color(0.1,0.,1.0,0));
  NCVvesselcap_log -> SetVisAttributes(NCVvesselcap_vis);

  G4VPhysicalVolume* NCVvesselcap1_phys
  = new G4PVPlacement(0,  // no rotation
			G4ThreeVector(0.,0.,NCVliquid_height/2. + NCVvesselcap_thickness/2.), // top of the cylinder
			NCVvesselcap_log,
			"NCVvesselcap1_phys",
			waterTank_log,          // mother
			false,
			0);
  
  G4VPhysicalVolume* NCVvesselcap2_phys
  = new G4PVPlacement(0,  // no rotation
			G4ThreeVector(0.,0., - NCVliquid_height/2. - NCVvesselcap_thickness/2.), // bottom of the cylinder
			NCVvesselcap_log,
			"NCVvesselcap2_phys",
			waterTank_log,          // mother
			false,
			0);
    
  
  // NCV liquid volume
     G4VSolid* NCVliquid_tub
      = new G4Tubs("NCVliquid_tub", //name
  		     0., // tube radius
  		     NCVliquid_radius, // r1, r2
		     NCVliquid_height/2., //Half height of the cylinder
  		     0.0, 2.0*M_PI         // phi0, delta_phi
  		     );

    G4LogicalVolume* NCVliquid_log
      = new G4LogicalVolume(NCVliquid_tub,
  			    G4Material::GetMaterial("NCVliquid"),
  			    "NCVliquid_log",
  			    0,
  			    0,
  			    0);
    G4VisAttributes* NCVliquid_vis
      = new G4VisAttributes(G4Color(0.1,1.0,0.,0));
    NCVliquid_log -> SetVisAttributes(NCVliquid_vis);
  
    G4VPhysicalVolume* NCVliquid_phys
     = new G4PVPlacement(0,  // no rotation
  			  G4ThreeVector(), // centered in the vessel
  			  NCVliquid_log,
			  "NCVliquid_phys",
  			  NCVvessel_log,          // mother
  			  false,
  			  0);
}

// ********************************************SciBooNE integration
///////////////////////////////////////////////////////////////
void WCLiteDetectorConstruction::DefineMRD(G4PVPlacement* expHall)
///////////////////////////////////////////////////////////////
{

  //#include "DefineMRD.icc"
}
// *******************************************/SciBooNE integration

void ConstructMRD(G4LogicalVolume* expHall_log, G4double tankRadius){     

//===============================================================================================================
  //MRD DETECTOR DEFINITION
//===============================================================================================================

  // Set up some sort of offset for the MRD from the tank. Leave the tank centred on the origin. 
  // N.B. Hall is 50*500*500m

  // need for the mother volume size	
  extern G4double mrdZlen; 
  extern G4double maxwidth;
  extern G4double maxheight;
  
  G4double mrdZoffset = tankRadius + (mrdZlen/2) + 2*cm;	
  // offset by half of length of both so edges touch + 2cm offset, with 1.98m tank radius puts MRD at z = 2.0*m.
  G4cout<<"MRD z start: "<<mrdZoffset-(mrdZlen/2.)<<" and total length: "<<mrdZlen<<G4endl;

  G4Box* totMRD_box = new G4Box("totMRDbox",maxwidth/2,maxheight/2,mrdZlen/2);
  G4LogicalVolume* totMRD_log = new G4LogicalVolume(totMRD_box, G4Material::GetMaterial("Air"),"totMRDlog",0,0,0);
  G4VPhysicalVolume* totMRD_phys = new G4PVPlacement(0,G4ThreeVector(0,0,mrdZoffset),totMRD_log,"totMRDphys",expHall_log,false,0);

  // ADD SCINTS & STEEL TO MRD
  // =========================
  PlacePaddles(totMRD_log);
  PlaceTapers(totMRD_log);
  PlaceLGs(totMRD_log);
  PlaceSteels(totMRD_log);

  // ADD ALU STRUCTURE
  // =================
  G4AssemblyVolume* aluMRDassembly = new G4AssemblyVolume();
  makeAlu(aluMRDassembly);
  G4RotationMatrix rot (0,0,0);
  G4ThreeVector trans(0,0,0);
  aluMRDassembly->MakeImprint(totMRD_log,trans,&rot);

  // END OF ALU STRUCTURE
  //=====================

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(false);
  totMRD_log->SetVisAttributes(simpleBoxVisAtt);

  // Retrieve Logical volumes and associate with sensitive detectors
  //================================================================
  extern G4LogicalVolume* hpaddle_log;
  extern G4LogicalVolume* vpaddle_log;
  extern G4LogicalVolume* taper_log;
  extern G4LogicalVolume* lg_log;
  
    // Get pointer to detector manager
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  // Create a new instance of MRD sensitive detector
  G4VSensitiveDetector* mrdSD = new MRDSD("MuonRangeDetector"); 
  // Register detector with manager
  SDman->AddNewDetector(mrdSD);
  // Attach detector to volume defining scintillator paddles
  hpaddle_log->SetSensitiveDetector(mrdSD);
  vpaddle_log->SetSensitiveDetector(mrdSD);
  taper_log->SetSensitiveDetector(mrdSD);
  //steel_log->SetSensitiveDetector(mrdSD);
  //aluStructsMRD_log->SetSensitiveDetector(mrdSD);
  
  // Create sensitive detector for light-guides
  G4VSensitiveDetector* mrdpmtSD = new mrdPMTSD("MRDPMTSD"); 
  // Register detector with manager
  SDman->AddNewDetector(mrdpmtSD);
  // Attach to light-guide logical volume
  //lg_log->SetSensitiveDetector(mrdpmtSD);	<- gets called manually and tracks get killed before they enter it....
  
  // Define Light Guide boundaries for optical photon detection
  // ==========================================================

  extern G4int numpaddlesperpanel;
  extern G4int numpanels;
  extern std::vector<G4VPhysicalVolume*> tapers_phys;
  extern std::vector<G4VPhysicalVolume*> lgs_phys;
  G4OpticalSurface* lgSurface_op = new G4OpticalSurface("lgopsurface",glisur, polished, dielectric_metal); 
  // must be dielectric_metal to invoke absorption/detection process - but is this overridden if both volumes have a ref index?
  // for dielectric_metal transmittance isn't possible, so either reflection or absorption with probability from mat. properties. 
  for(G4int i=0;i<(numpaddlesperpanel*numpanels);i++){
  	G4LogicalBorderSurface* lgSurface_log = new G4LogicalBorderSurface("lgborderlog",tapers_phys.at(i),lgs_phys.at(i),lgSurface_op);
  }
  
  const G4int lgmptentries = 2;
  G4double Ephoton[lgmptentries] = {2.038*eV, 4.144*eV};
  G4double lgsurf_EFF[lgmptentries]={1.,1.};
  G4double lgsurf_REFL[lgmptentries]={0.,0.};
  G4MaterialPropertiesTable* lgsurf_MPT = new G4MaterialPropertiesTable();
  lgsurf_MPT->AddProperty("EFFICIENCY",  Ephoton, lgsurf_EFF,  lgmptentries);
  lgsurf_MPT->AddProperty("REFLECTIVITY",Ephoton, lgsurf_REFL, lgmptentries);
  lgSurface_op->SetMaterialPropertiesTable(lgsurf_MPT);


}

