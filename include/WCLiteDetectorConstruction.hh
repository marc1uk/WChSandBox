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
// $Id: WCLiteDetectorConstruction.hh,v 1.5 2006/06/29 17:53:55 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WCLiteDetectorConstruction_h
#define WCLiteDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


// **************SciBooNE integration
//#include "SBsimInputCard.hh"	--crossed out by marcus, try to simplify by extracting just mrd struct
#include "G4UserLimits.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"

#include "SBsimMRDDB.hh"
static const G4double INCH = 2.54*cm;
// *************/SciBooNE integration

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WCLiteDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    WCLiteDetectorConstruction();
   ~WCLiteDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();
    
    inline SBsimMRDDB*     GetMRDDB()     {return mrddb;};	// marcus: added

  private:
    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;

    G4double tank_x;
    G4double tank_y;
    G4double tank_z;

    G4double bubble_x;
    G4double bubble_y;
    G4double bubble_z;

    void  ConstructMaterials();
    
    // ****************SciBooNE integration
    
    SBsimMRDDB*     mrddb;	// marcus: added
      
      
    void DefineMRD(G4PVPlacement* expHall);
    G4Material *Air, *Al, *Fe, *Pb;
    G4Material *Steel, *Scinti, *Lucite;
    G4Material *ECScinti, *MRDIron;
    G4Material *Plywood;
    
	// for MRD
	G4VSolid *MRDIron_Solid[12];
	G4VSolid *MRDVScinti_Solid, *MRDHScinti_Solid;
	G4VSolid *MRDTScinti_Solid;
	G4VSolid *MRDLG_Solid;

	G4LogicalVolume *MRDIron_LV[12];
	G4LogicalVolume *MRDVScinti_LV, *MRDHScinti_LV;
	G4LogicalVolume *MRDTScinti_LV;
	G4LogicalVolume *MRDLG_LV;

	// for MRD Al support
	G4SubtractionSolid *MRDAlV1_Solid;
	G4SubtractionSolid *MRDAlV2_Solid;
	G4SubtractionSolid *MRDAlV3_Solid;
	G4SubtractionSolid *MRDAlV4_Solid;
	G4SubtractionSolid *MRDAlV5_Solid;
	G4SubtractionSolid *MRDAlH1_Solid;
	G4SubtractionSolid *MRDAlH2_Solid;
	G4SubtractionSolid *MRDAlH3_Solid;

	G4VSolid *MRDAlV1_Outer;
	G4VSolid *MRDAlV2_Outer;
	G4VSolid *MRDAlV3_Outer;
	G4VSolid *MRDAlV4_Outer;
	G4VSolid *MRDAlV5_Outer;
	G4VSolid *MRDAlH1_Outer;
	G4VSolid *MRDAlH2_Outer;
	G4VSolid *MRDAlH3_Outer;

	G4VSolid *MRDAlV1_Inner;
	G4VSolid *MRDAlV2_Inner;
	G4VSolid *MRDAlV3_Inner;
	G4VSolid *MRDAlV4_Inner;
	G4VSolid *MRDAlV5_Inner;
	G4VSolid *MRDAlH1_Inner;
	G4VSolid *MRDAlH2_Inner;
	G4VSolid *MRDAlH3_Inner;

	G4LogicalVolume *MRDAlV1_LV;
	G4LogicalVolume *MRDAlV2_LV;
	G4LogicalVolume *MRDAlV3_LV;
	G4LogicalVolume *MRDAlV4_LV;
	G4LogicalVolume *MRDAlV5_LV;
	G4LogicalVolume *MRDAlH1_LV;
	G4LogicalVolume *MRDAlH2_LV;
	G4LogicalVolume *MRDAlH3_LV;

	G4AssemblyVolume * MRDAlSupportV;
	G4AssemblyVolume * MRDAlSupportH;
    // *********/SciBooNE integration
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*WCLiteDetectorConstruction_h*/



