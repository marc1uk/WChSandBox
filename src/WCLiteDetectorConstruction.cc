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

// **********SciBooNE integration
#include "SBsimMRDSD.hh"
#include "SBsimMRDDB.hh"
// **********/SciBooNE integration
void ConstructMRD(G4LogicalVolume* expHall_log, G4double mrdZoffset);

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
  G4Box* expHall_box = new G4Box("Hall",20*m,20*m,25*m);		//just..nicer for now.
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
  G4bool WCAddGd = true;

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
  	

  G4Tubs* waterTank_box = new G4Tubs("waterTank",0.*m,5.5*m,11.*m,0.*deg,360.*deg);
  //  G4Box* expHall_box = new G4Box("waterTank",1.5*m,1.5*m,1.5*m); ANNIE
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

	


  G4double tankRadius = waterTank_box->GetOuterRadius();
  ConstructMRD(expHall_log, tankRadius);			// marcus' MRD construction
  expHall_log->SetVisAttributes (G4VisAttributes::Invisible);	// set hall volume invisible
  
  // ******SciBooNE integration
  G4PVPlacement* expHall= new G4PVPlacement(0, G4ThreeVector(), "EXP_HALL_PV", expHall_log, 0, false, 0);
  //DefineMRD(expHall);	
  // *****/SciBooNE integration
  
  //always return the physical World
  return expHall_phys;

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
extern G4int numpaddlesperpanel;	// paddles per scintillator panel
extern G4int numpanels;			// scintillator panels
extern G4int numrpcs;			// rpc panels
extern G4int numplates;			// steel plates
extern G4int numalustructs;		// number of supporting structs. We may be dropping one as we have fewer scintillators?
// need to import these (from MRDDetectorConstruction.cc translation unit, via extern) so we know how many objects are needed to parameterise the MRD. 
	
extern G4double mrdZlen; 
extern G4double maxwidth;
extern G4double maxheight;
// need for the mother volume (scintandsteelMRD_log) size

        
// Define solids
//==============  
// G4Box* variableName = new G4Box("SolidName", x_length, y_length, z_length);
  	
  	G4Box* generic_box = new G4Box("genericBox",1*cm,1*cm,1*cm); // dims arbitrary - overwritten by parameterisation
  	G4Box* totMRD_box = new G4Box("totMRD",maxwidth/2,maxheight/2,mrdZlen/2);			

// Define logical volumes
//=======================
// G4LogicalVolume* variableName = new G4LogicalVolume(solidVariableName, Material, "logicalVolName");

	G4LogicalVolume* generic_log
	= new G4LogicalVolume(generic_box,G4Material::GetMaterial("Air"),"genericBox",0,0,0);
/*	
	G4LogicalVolume* sciMRDhpaddle_log
	= new G4LogicalVolume(sciMRDhpaddle_box,G4Material::GetMaterial("Scintillator"),"scintHpaddle",0,0,0);
	
	G4LogicalVolume* sciMRDvpaddle_log
	= new G4LogicalVolume(sciMRDvpaddle_box,G4Material::GetMaterial("Scintillator"),"scintVpaddle",0,0,0);
*/	
	G4LogicalVolume* scintAndSteelMRD_log
	= new G4LogicalVolume(totMRD_box, G4Material::GetMaterial("Air"),"totMRD",0,0,0);
	
	G4LogicalVolume* aluStructsMRD_log
	= new G4LogicalVolume(totMRD_box, G4Material::GetMaterial("Air"),"totMRD2",0,0,0);
	// need a second mother volume for alu struct imprint; a mother may contain a repeated (replica/parameterisation) volume, OR placements (including imprints), not both.


// MRD PHYSICAL PLACEMENT
// Paramaterisations are placed within the MRD logical volume; we need to place the logical volume. Can do this before placing parameterisations. 

  	G4double mrdZoffset = tankRadius + totMRD_box->GetZHalfLength()+1*cm;	//offset by half of length of both so edges touch
  	G4VPhysicalVolume* scintAndSteelMRD_phys = new G4PVPlacement(0,G4ThreeVector(0,0,mrdZoffset),scintAndSteelMRD_log,"MRD",expHall_log,false,0);
  	G4VPhysicalVolume* aluStructsMRD_phys = new G4PVPlacement(0,G4ThreeVector(0,0,mrdZoffset),aluStructsMRD_log,"MRD",expHall_log,false,0);

//------------------------------------------
// ADD SCINTS & STEEL TO MRD PLACEMENT
// ===================================
	// 2. Declare an instance of the parameterisation class
	// This will encapsulate functions that retrieves specifics of each member (e.g. position)
	G4VPVParameterisation* mrdObjectPicker = new mrdParamaterisedModel();//mrdParamaterisedModel
	
	// 3. Create physical volume from G4VParameterised class - this will "represent" all physical steel plates
	G4VPhysicalVolume* mrdModel = new G4PVParameterised("mrdModel", generic_log, scintAndSteelMRD_log, kZAxis, 
	(numplates+(numpaddlesperpanel*numpanels)), mrdObjectPicker);
	
// END OF ADDING SCINTS & STEEL
// ===================================

// ADD ALU STRUCTURE
// ====================
	G4AssemblyVolume* aluMRDassembly = new G4AssemblyVolume();
	makeAlu(aluMRDassembly);
	G4RotationMatrix  rot (0,0,0);
	G4ThreeVector  trans(0,0,0);
	aluMRDassembly->MakeImprint(aluStructsMRD_log,trans,&rot);
// END OF ALU STRUCTURE
//=====================

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);		
  scintAndSteelMRD_log->SetVisAttributes(simpleBoxVisAtt);
  aluStructsMRD_log->SetVisAttributes(simpleBoxVisAtt);
  
  
  // Create a new BetamTestEmCalorimeter sensitive detector
  G4VSensitiveDetector* mrdSD = new MRDSD("MuonRangeDetector"); 
 
  // Get pointer to detector manager
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
 
 // Register detector with manager
  SDman->AddNewDetector(mrdSD);
 
  // Attach detector to volume defining scintillator paddles (also steel panels through parameterisation)
  generic_log->SetSensitiveDetector(mrdSD);

}


