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
#include "G4SystemOfUnits.hh"

#include "SBsimMRDDB.hh"
static const G4double INCH = 2.54*cm;
// *************/SciBooNE integration

    
  // *************** WCSim PMT integration
  // =====================================
  
#include <map>
#include <vector>
#include "G4OpticalSurface.hh"
//#include <hash_map.h> // warning : hash_map is not part of the standard
#include <ext/hash_map>
#include "WCSimPmtInfo.hh"
#include "WCSimPMTObject.hh"

//class G4LogicalVolume;
//class G4AssemblyVolume;
//class G4VPhysicalVolume;
class WCSimWCSD;

namespace __gnu_cxx  {
  template<> struct hash< std::string >
  {
    size_t operator()( const std::string& x ) const
    {
      return hash< const char* >()( x.c_str() );
    }
  };
}

using __gnu_cxx::hash;
using __gnu_cxx::hashtable;
using __gnu_cxx::hash_map;
using __gnu_cxx::hash_multimap;

  // End WCSim PMT Integration
  // ==========================

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class WCLiteDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    WCLiteDetectorConstruction(G4String gdmlfname = "../WChSandBox_v1/src/annie_v01.gdml");
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
    
    //G4double tankinnerRadius, tankouterRadius, tankhz;

    G4double bubble_x;
    G4double bubble_y;
    G4double bubble_z;

    void  ConstructMaterials();
    
    //================ Rob Hatcher integration
    G4String GDMLFilename;
    G4int doOverlapCheck;
    
    // ****************SciBooNE integration
    
    SBsimMRDDB*     mrddb;	// marcus: added
      
      
    void DefineMRD(G4PVPlacement* expHall);
    void ConstructMRD(G4LogicalVolume* expHall_log, G4VPhysicalVolume* expHall_phys);
    void ConstructVETO(G4LogicalVolume* expHall_log, G4VPhysicalVolume* expHall_phys);
    void ConstructNCV(G4LogicalVolume* waterTank_log);
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

		// *************** WCSim PMT integration
		// =====================================
		public: 
		G4double GetPMTSize()           {return WCPMTRadius;}
		G4String GetPMTName()           {return WCPMTName;}
		G4int    GetTotalNumPmts()      {return totalNumPMTs;}
		G4int    GetPMT_QE_Method()     {return PMT_QE_Method;}
		G4int    UsePMT_Coll_Eff()      {return PMT_Coll_Eff;}
		G4String GetDetectorName()      {return WCDetectorName;}
		
		G4float GetPMTQE(G4String,G4float, G4int, G4float, G4float, G4float);
		G4float GetPMTCollectionEfficiency(G4float theta_angle, G4String CollectionName) {
		  return GetPMTPointer(CollectionName)->GetCollectionEfficiency(theta_angle); 
		};
		
		void SetANNIEGeometry();
		
		WCSimPMTObject *CreatePMTObject(G4String, G4String);
		std::map<G4String, WCSimPMTObject*>  CollectionNameMap; 
		WCSimPMTObject * PMTptr;
		
		void SetPMTPointer(WCSimPMTObject* PMT, G4String CollectionName){
		  CollectionNameMap[CollectionName] = PMT;
		}
		WCSimPMTObject* GetPMTPointer(G4String CollectionName){
		  PMTptr = CollectionNameMap[CollectionName];
		  if (PMTptr == NULL) {G4cout << CollectionName << " is not a recognized hit collection. Exiting WCSim." << G4endl; exit(1);}
		  return PMTptr;
		}
		
		// Related to the WC tube ID
		static G4int GetTubeID(std::string tubeTag){return tubeLocationMap[tubeTag];}
		static G4Transform3D GetTubeTransform(int tubeNo){return tubeIDMap[tubeNo];}
		void   SetPMT_QE_Method(G4int choice){PMT_QE_Method = choice;}
		void   SetPMT_Coll_Eff(G4int choice){PMT_Coll_Eff = choice;}
		
		void   SetVis_Choice(G4String choice){Vis_Choice = choice;}
		G4String GetVis_Choice() {return Vis_Choice;}
		
		std::vector<WCSimPmtInfo*>* Get_Pmts() {return &fpmts;}

		G4String GetIDCollectionName(){return WCIDCollectionName;}
		
		private:
		WCSimWCSD* aWCPMT;
		G4LogicalVolume* ConstructPMT(G4String,G4String);
		
		// Funcs for traversing the geometry
		typedef void (WCLiteDetectorConstruction::*DescriptionFcnPtr) (G4VPhysicalVolume*, int, int, const G4Transform3D&);
		void DumpGeometryTableToFile();
		void PrintGeometryTree(G4VPhysicalVolume*, int, int, const G4Transform3D&);
		void TraverseReplicas(G4VPhysicalVolume*, int, const G4Transform3D&, DescriptionFcnPtr);
		void DescribeAndRegisterPMT(G4VPhysicalVolume*, int, int, const G4Transform3D&);
		void DescribeAndDescendGeometry(G4VPhysicalVolume*, int, int, const G4Transform3D&, DescriptionFcnPtr);
		void GetWCGeom(G4VPhysicalVolume*, int, int /*ReplicaNo*/, const G4Transform3D&);
		G4double GetGeo_Dm(G4int);
		
		// XQ 08/17/10
		//   PMT_QE_Method == 1
		//   Only use it in the stacking function (no WLS)
		//   PMT_QE_Method == 2
		//   Use Part of it in the stacking function (MAX QE)
		//   Put the rest of it in the sensitive detector according to QE/Max_QE
		//   PMT_QE_Method == 3
		//   Put all of it in the sensitive detector according to QE
		//   Good for low energy photons
		G4int PMT_QE_Method;
		
		//XQ 03/31/11
		// 0 to not use collection efficiency
		// 1 to use
		G4int PMT_Coll_Eff;

		//NP 06/17/15
		// "OGLSX" for classic visualization
		// "RayTracer" for RayTracer visualization
		G4String Vis_Choice;
		
		// Hit collection name parameters
		G4String WCDetectorName;
		G4String WCIDCollectionName;
		G4String WCODCollectionName;
		
		// WC PMT parameters
		G4String WCPMTName;
		typedef std::pair<G4String, G4String> PMTKey_t;
		typedef std::map<PMTKey_t, G4LogicalVolume*> PMTMap_t;
		
		// for dumping geometry info to file. Even if we don't use it, the class definitions need these.
		std::ofstream geoFile;   // File for text output
		G4int totalNumPMTs;      // The number of PMTs for this configuration     
		G4double WCCylInfo[3];    // Info for the geometry tree: radius & length or mail box, length, width and depth
		G4double WCPMTSize;       // Info for the geometry tree: pmt size
		G4ThreeVector WCOffset;   // Info for the geometry tree: WC center offset
		G4double innerradius;
		
		
		static PMTMap_t PMTLogicalVolumes;
		G4double WCPMTRadius;
		G4double WCPMTExposeHeight;
		// Tube map information
		static std::map<int, G4Transform3D> tubeIDMap;
		//  static std::map<int, cyl_location> tubeCylLocation;
		static hash_map<std::string, int, hash<std::string> >  tubeLocationMap; 
		std::vector<WCSimPmtInfo*> fpmts;
		  //Water, Blacksheet surface
  G4OpticalSurface * OpWaterBSSurface;
  //Glass, Cathode surface in PMTs
  G4OpticalSurface * OpGlassCathodeSurface;
  //Tyvek surface - jl145
  G4OpticalSurface * OpWaterTySurface;
		
		G4bool isHyperK;
		G4bool checkOverlaps;
		
		// End WCSim PMT Integration
		// ==========================

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*WCLiteDetectorConstruction_h*/



