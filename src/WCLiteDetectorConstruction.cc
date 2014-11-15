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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteDetectorConstruction::WCLiteDetectorConstruction()
{
  expHall_x = 50*m;
  expHall_y = expHall_z = 500*m;
  //  expHall_x = 10*m;
  //  expHall_y = expHall_z = 10*m;
  //  tank_x    = tank_y    = tank_z    = 5*m;
  //  bubble_x  = bubble_y  = bubble_z  = 0.5*m;
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

 //-------------------------------------------------------------------------
  // Set detector geometry
  //
  // Modified 2014-06-13 by MSM ?!? 
  //  - TITUS geometry added
  
  // ANNIE geometry (3m cube)
  G4Box* expHall_box = new G4Box("waterTank",1.5*m,1.5*m,1.5*m); 

  // TITUS geometry (cylinder w/ 5.5m radius, 22m length)
  // G4Tubs* expHall_box = new G4Tubs("waterTank",0.*m,5.5*m,11.*m,0.*deg,360.*deg);

  //G4Sphere* expHall_box = new G4Sphere("waterTank",0.*m,10.*m,0.*deg,360.*deg,0.*deg,180.*deg);
  //-------------------------------------------------------------------------  
  


  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,G4Material::GetMaterial(water),"waterTank",0,0,0);

  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"waterTank",0,false,0);

 
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

//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
