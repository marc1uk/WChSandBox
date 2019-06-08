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
#include "FACCSD.hh"
#include "faccPMTSD.hh"
#include "NCVSD.hh"
#include "G4SystemOfUnits.hh"
#include "G4GDMLParser.hh"

// **********SciBooNE integration
#include "SBsimMRDSD.hh"
#include "SBsimMRDDB.hh"
// **********/SciBooNE integration

G4double tankouterRadius= 1.524*m;	// 120" exactly (TSW blueprint) = 3.048m diameter
G4double tankhz = 1.98*m;		// 13ft exactly (TSW blueprint) = 3.96m tall
G4double tankzoffset = 15.70*cm;	//15.70*cm
// tank is 7GA (7-gauge steel) = (0.1793"normal sheet steel or) 0.1875" stainless steel = 4.7625mm thick 'A36' steel (density 7,800 kg/m3 )
// ^ gauge represents the thickness of 1 ounce of copper rolled out to an area of 1 square foot
//void ConstructMRD(G4LogicalVolume* expHall_log, G4VPhysicalVolume* expHall_phys);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteDetectorConstruction::WCLiteDetectorConstruction(G4String gdmlfname) // skipping arg from SciBooNE integration
{
  expHall_x = 50*m;
  expHall_y = expHall_z = 500*m;
  GDMLFilename = gdmlfname;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WCLiteDetectorConstruction::~WCLiteDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WCLiteDetectorConstruction::Construct()
{
  
  ConstructMaterials();
  
  //Decide if adding Gd
  G4cout<<"Tank is full of ";
  G4bool isWater = true;
  G4bool WCAddGd = false;
  G4bool isNCV = true; // Construct the Neutron Capture Volume
  G4String Scintillator = "Scintillator";
  G4String water = "Water";
  if (!isWater)
    { water = "Scintillator";
    G4cout<<"***SCINTILLATOR***"<<G4endl; }
  if (WCAddGd && isWater)
    {water = "Doped Water";
    G4cout<<"***Dope Water***"<<G4endl;}
  if (WCAddGd && !isWater)
    {water = "Doped Scintillator";
    G4cout<<"***Dope Scintillator***"<<G4endl;}
  if (water=="Water"){
  G4cout<<"refreshing water"<<G4endl;}
  
  
//  //  ===== Rob Hatcher's integration
//  G4GDMLParser parser;  // Read GDML file
//  parser.SetOverlapCheck(0);
//  G4cout << "Read " << GDMLFilename << " (overlap check = " << (doOverlapCheck?"true":"false") << ")" << G4endl;
//  parser.Read ( GDMLFilename );
//  G4VPhysicalVolume* expHall_phys = parser.GetWorldVolume();
//  // if we wish to set visualisation properties, get logical volume
//  G4LogicalVolume* expHall_log = expHall_phys->GetLogicalVolume();
//  //  ===== Rob Hatcher's integration
  
    // Create Experimental Hall
  //G4Box* expHall_box = new G4Box("Hall",expHall_x,expHall_y,expHall_z);
  G4Box* expHall_box = new G4Box("Hall",90*m,90*m,25*m);		//just..nicer for now.
  G4LogicalVolume* expHall_log = new G4LogicalVolume(expHall_box,G4Material::GetMaterial("Air"),"Hall",0,0,0);
  G4VPhysicalVolume* expHall_phys = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"Hall",0,false,0);
  
  // Create Water Tank
  //G4Tubs* aTube = new G4Tubs("Name", innerRadius, outerRadius, hz, startAngle, spanningAngle);
  G4Tubs* waterTank_tubs = new G4Tubs("waterTank",0.*m,tankouterRadius,tankhz,0.*deg,360.*deg);	//prev dims: 5.5m radius, 11m height.
  G4LogicalVolume* waterTank_log = new G4LogicalVolume(waterTank_tubs,G4Material::GetMaterial(water),"waterTank",0,0,0);
  G4RotationMatrix* rotm = new G4RotationMatrix();				// declare a rotation matrix
  rotm->rotateX(90*deg); 							// rotate Y
  G4VPhysicalVolume* waterTank_phys = new G4PVPlacement(rotm,G4ThreeVector(0,-144.64875*mm,tankouterRadius+tankzoffset),waterTank_log,"waterTank",expHall_log,false,0);

//  // Create MRD
//  ConstructMRD(expHall_log, expHall_phys);			// marcus' MRD construction
  
  // ===== SciBooNE integration
  //G4PVPlacement* expHall= new G4PVPlacement(0, G4ThreeVector(), "EXP_HALL_PV", expHall_log, 0, false, 0);
  DefineMRD((G4PVPlacement*)expHall_phys);	// used to take in expHall above
  //  ===== /SciBooNE integration
  
  // Create FACC
  ConstructVETO(expHall_log, expHall_phys);
  
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
  
  expHall_log->SetVisAttributes (G4VisAttributes::Invisible);	// set hall volume invisible
  //always return the physical World
  return expHall_phys;

}

// ********************************************SciBooNE integration
///////////////////////////////////////////////////////////////
void WCLiteDetectorConstruction::DefineMRD(G4PVPlacement* expHall)
///////////////////////////////////////////////////////////////
{
  #include "DefineMRD.icc"
  
  // ======================================================================================================
  // Below are extensions to SciBooNE code that add further SDs for paddle cladding and PMTs at paddle ends
  // ======================================================================================================
  
  // Create sensitive detector for pmts that will detect photon hits at the ends of the paddles
  G4VSensitiveDetector* mrdpmtSD = new mrdPMTSD("MRDPMTSD"); 
  // Register detector with manager
  SDman->AddNewDetector(mrdpmtSD);
  // gets called manually and tracks get killed before they enter it; don't need to associate with any logical volumes
  
  // Define Paddle Boundaries for MRD PMT photon detection
  // =====================================================
  G4int numpaddlesperpanelh=13;
  G4int numpaddlesperpanelv=15;
  G4int numpanels = mrdmod->NLayer+1;
  
  // Define efficiency for reflection & detection between paddles and PMTs at their ends
  // ===================================================================================
  G4OpticalSurface* lgSurface_op = new G4OpticalSurface("lgopsurface",glisur, polished, dielectric_metal);  
  const G4int lgmptentries = 2;
  G4double ELGphoton[lgmptentries] = {2.038*eV, 4.144*eV};
  G4double lgsurf_EFF[lgmptentries]={0.95,0.95};
  G4double lgsurf_REFL[lgmptentries]={0.,0.};
  G4MaterialPropertiesTable* lgsurf_MPT = new G4MaterialPropertiesTable();
  lgsurf_MPT->AddProperty("EFFICIENCY",  ELGphoton, lgsurf_EFF,  lgmptentries);
  lgsurf_MPT->AddProperty("REFLECTIVITY",ELGphoton, lgsurf_REFL, lgmptentries);
  lgSurface_op->SetMaterialPropertiesTable(lgsurf_MPT);
  
  // must be dielectric_metal to invoke absorption/detection process - but is this overridden if both volumes have a ref index?
  // for dielectric_metal transmittance isn't possible, so either reflection or absorption with probability from mat. properties. 
  for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
      G4LogicalBorderSurface* lgSurface_log = new G4LogicalBorderSurface("lgborderlog",tapers_phys.at(i),lgs_phys.at(i),lgSurface_op);
  }
  
  // Define Mylar cladding on paddles
  // ================================
  G4OpticalSurface* scintSurface_op = new G4OpticalSurface("mylarSurface",glisur, ground, dielectric_metal);
  const G4int mylarmptentries = 2;
  G4double mylar_Energy[mylarmptentries] = {2.0*eV, 3.6*eV};
  G4double mylar_REFL[mylarmptentries] = {0.9,0.9};
  G4double mylar_EFFI[mylarmptentries] = {0.0, 0.0};
  G4MaterialPropertiesTable* MPTmylarSurface = new G4MaterialPropertiesTable();
  MPTmylarSurface->AddProperty("REFLECTIVITY",mylar_Energy,mylar_REFL,mylarmptentries);
  MPTmylarSurface->AddProperty("EFFICIENCY",mylar_Energy,mylar_EFFI,mylarmptentries);
  scintSurface_op->SetMaterialPropertiesTable(MPTmylarSurface);
  
  // Place cladding on MRD paddles
  for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
    G4LogicalBorderSurface* scintSurface_log
      = new G4LogicalBorderSurface("scintcladdinglog",lgs_phys.at(i),expHall,scintSurface_op);
  }
  
}
// *******************************************/SciBooNE integration

void WCLiteDetectorConstruction::ConstructMRD(G4LogicalVolume* expHall_log, G4VPhysicalVolume* expHall_phys){

  //===============================================================================================================
  // MRD DETECTOR DEFINITION
  //===============================================================================================================

  // Set up some sort of offset for the MRD from the tank. Leave the tank centred on the origin. 
  // N.B. Hall is 50*500*500m

  // need for the mother volume size
  extern G4double mrdZlen; 
//  extern G4double maxwidth;
//  extern G4double maxheight;
   
  // offset MRD by half of length of both so edges touch + 2cm offset, with 1.52m tank radius puts MRD at z = 1.54*m.
  G4double mrdZoffset = (2*tankouterRadius) + tankzoffset + (mrdZlen/2.) + 5*cm;
  G4cout<<"MRD z start: "<<(tankouterRadius + 2*cm)<<" and total length: "<<mrdZlen<<G4endl;

  extern G4Box* totMRD_box;
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
  
  // ADD VIS ATTRIBS
  //================
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
  
  // Create sensitive detector for pmts that will detect photon hits at the ends of the paddles
  G4VSensitiveDetector* mrdpmtSD = new mrdPMTSD("MRDPMTSD"); 
  // Register detector with manager
  SDman->AddNewDetector(mrdpmtSD);
  // gets called manually and tracks get killed before they enter it; don't need to associate with any logical volumes
  
  // Define Paddle Boundaries for MRD PMT photon detection
  // =====================================================
  extern G4int numpaddlesperpanelh;
  extern G4int numpaddlesperpanelv;
  extern G4int numpanels;
  extern std::vector<G4VPhysicalVolume*> tapers_phys;
  extern std::vector<G4VPhysicalVolume*> lgs_phys;
  
  // Define efficiency for reflection & detection between paddles and PMTs at their ends
  // ===================================================================================
  G4OpticalSurface* lgSurface_op = new G4OpticalSurface("lgopsurface",glisur, polished, dielectric_metal);  
  const G4int lgmptentries = 2;
  G4double ELGphoton[lgmptentries] = {2.038*eV, 4.144*eV};
  G4double lgsurf_EFF[lgmptentries]={0.95,0.95};
  G4double lgsurf_REFL[lgmptentries]={0.,0.};
  G4MaterialPropertiesTable* lgsurf_MPT = new G4MaterialPropertiesTable();
  lgsurf_MPT->AddProperty("EFFICIENCY",  ELGphoton, lgsurf_EFF,  lgmptentries);
  lgsurf_MPT->AddProperty("REFLECTIVITY",ELGphoton, lgsurf_REFL, lgmptentries);
  lgSurface_op->SetMaterialPropertiesTable(lgsurf_MPT);
  
  // must be dielectric_metal to invoke absorption/detection process - but is this overridden if both volumes have a ref index?
  // for dielectric_metal transmittance isn't possible, so either reflection or absorption with probability from mat. properties. 
  for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
  	G4LogicalBorderSurface* lgSurface_log = new G4LogicalBorderSurface("lgborderlog",tapers_phys.at(i),lgs_phys.at(i),lgSurface_op);
  }
  
  // Define Mylar cladding on paddles
  // ================================
  G4OpticalSurface* scintSurface_op = new G4OpticalSurface("mylarSurface",glisur, ground, dielectric_metal);
  const G4int mylarmptentries = 2;
  G4double mylar_Energy[mylarmptentries] = {2.0*eV, 3.6*eV};
  G4double mylar_REFL[mylarmptentries] = {0.9,0.9};
  G4double mylar_EFFI[mylarmptentries] = {0.0, 0.0};
  G4MaterialPropertiesTable* MPTmylarSurface = new G4MaterialPropertiesTable();
  MPTmylarSurface->AddProperty("REFLECTIVITY",mylar_Energy,mylar_REFL,mylarmptentries);
  MPTmylarSurface->AddProperty("EFFICIENCY",mylar_Energy,mylar_EFFI,mylarmptentries);
  scintSurface_op->SetMaterialPropertiesTable(MPTmylarSurface);
	
	// Place cladding on MRD paddles
  for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
  		G4LogicalBorderSurface* scintSurface_log
  		 = new G4LogicalBorderSurface("scintcladdinglog",lgs_phys.at(i),expHall_phys,scintSurface_op);
  }

}

void WCLiteDetectorConstruction::ConstructVETO(G4LogicalVolume* expHall_log, G4VPhysicalVolume* expHall_phys){
  //===============================================================================================================
  // VETO WALL DEFINITION
  //===============================================================================================================
  extern G4double vetoZlen;
  
  
  G4double vetoZoffset = -(0*tankouterRadius) + (vetoZlen/2.) + 2*cm;	// by definition of rob's geometry, FACC is ~z=0
  G4cout << "Veto z start: " << vetoZoffset-(vetoZlen/2.) << " and total length: "<< vetoZlen << G4endl;
  
  extern G4Box* totVeto_box;
  G4LogicalVolume* totVeto_log = new G4LogicalVolume(totVeto_box, G4Material::GetMaterial("Air"), "totVetolog",0,0,0);
  G4VPhysicalVolume* totVeto_phys = new G4PVPlacement(0,G4ThreeVector(0,-137.64875*mm,vetoZoffset),totVeto_log,"totVetoPhys",expHall_log,false,0); 
  
  
  // FACC VETO
  // =========
  PlaceVetoPaddles(totVeto_log); 
  
  // ADD VIS ATTRIBS
  //================
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(false);
  totVeto_log->SetVisAttributes(simpleBoxVisAtt);
  
  // Retrieve Logical volumes and associate with sensitive detectors
  //================================================================
  extern G4LogicalVolume* vetoPaddle_log;
  extern G4LogicalVolume* vetol2Paddle_log;
  
  // Get pointer to detector manager
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  // Sensitive detector for Veto paddles
  G4VSensitiveDetector* faccSD = new FACCSD("FrontAntiCoincidenceCounter");
  SDman->AddNewDetector(faccSD);
  vetoPaddle_log->SetSensitiveDetector(faccSD);
  vetol2Paddle_log->SetSensitiveDetector(faccSD);
  
  // Sensitive detector for Veto PMTs
  G4VSensitiveDetector* faccpmtSD = new faccPMTSD("FACCPMTSD"); 
  // Register detector with manager
  SDman->AddNewDetector(faccpmtSD);
  // gets called manually and tracks get killed before they enter it; don't need to associate with any logical volumes
  
  // Define Paddle Boundaries for MRD PMT photon detection
  // =====================================================
  extern G4int numvetopaddles;
  extern std::vector<G4VPhysicalVolume*> vetopaddles_phys;
  
  // Define 'Surface' volume boundaries for Veto PMT detection
  // =========================================================
  extern std::vector<G4VPhysicalVolume*> vetosurfaces_phys;
  
  // Define efficiency for reflection & detection between paddles and PMTs at their ends
  // ===================================================================================
  G4OpticalSurface* vetoSensSurface_op = new G4OpticalSurface("vetosenssurfaceop",glisur, polished, dielectric_metal);
  const G4int vetosensmptentries = 2;
  G4double EVetophoton[vetosensmptentries] = {2.038*eV, 4.144*eV};
  G4double vetosenssurf_EFF[vetosensmptentries]={0.95,0.95};	// should be 5% efficient?!!
  G4double vetosenssurf_REFL[vetosensmptentries]={0.,0.};
  G4MaterialPropertiesTable* vetosenssurf_MPT = new G4MaterialPropertiesTable();
  vetosenssurf_MPT->AddProperty("EFFICIENCY",  EVetophoton, vetosenssurf_EFF,  vetosensmptentries);
  vetosenssurf_MPT->AddProperty("REFLECTIVITY",EVetophoton, vetosenssurf_REFL, vetosensmptentries);
  vetoSensSurface_op->SetMaterialPropertiesTable(vetosenssurf_MPT);
  
  // must be dielectric_metal to invoke absorption/detection process - but is this overridden if both volumes have a ref index?
  // for dielectric_metal transmittance isn't possible, so either reflection or absorption with probability from mat. properties. 
  for(G4int i=0;i<(numvetopaddles);i++){
  	G4LogicalBorderSurface* vetoBorderSurface_log = new G4LogicalBorderSurface("vetoborderlog",vetopaddles_phys.at(i),vetosurfaces_phys.at(i),vetoSensSurface_op);
  }
  
  // Define Mylar cladding on paddles
  // ================================
  G4OpticalSurface* scintSurface_op = new G4OpticalSurface("mylarSurface",glisur, ground, dielectric_metal);
  const G4int mylarmptentries = 2;
  G4double mylar_Energy[mylarmptentries] = {2.0*eV, 3.6*eV};
  G4double mylar_REFL[mylarmptentries] = {0.9,0.9};
  G4double mylar_EFFI[mylarmptentries] = {0.0, 0.0};
  G4MaterialPropertiesTable* MPTmylarSurface = new G4MaterialPropertiesTable();
  MPTmylarSurface->AddProperty("REFLECTIVITY",mylar_Energy,mylar_REFL,mylarmptentries);
  MPTmylarSurface->AddProperty("EFFICIENCY",mylar_Energy,mylar_EFFI,mylarmptentries);
  scintSurface_op->SetMaterialPropertiesTable(MPTmylarSurface);
  
  // Place cladding on Veto paddles
  for(G4int i=0;i<(numvetopaddles);i++){
  		G4LogicalBorderSurface* vetoSurface_log
  		 = new G4LogicalBorderSurface("vetocladdinglog",vetopaddles_phys.at(i),expHall_phys,scintSurface_op);
  }
  
}

void WCLiteDetectorConstruction::ConstructNCV(G4LogicalVolume* waterTank_log){
  //===============================================================================================================
  //NEUTRON CAPTURE VOLUME DEFINITION
  //===============================================================================================================
  // NCV is defined as two cylinders, one filled with liquid, the other with acrylic. Two 'caps' are added to the 
  // top of the NCV as well. 
  // The metal structure will be added as well (from pictures) since it will induce n-Fe captures
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
    = new G4VisAttributes(G4Color(0.1,0.,1.0,1));
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
    = new G4VisAttributes(G4Color(0.1,0.,1.0,1));
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
      = new G4VisAttributes(G4Color(0.1,1.0,0.,1));
    NCVliquid_log -> SetVisAttributes(NCVliquid_vis);
  
    G4VPhysicalVolume* NCVliquid_phys
     = new G4PVPlacement(0,  // no rotation
  			  G4ThreeVector(), // centered in the vessel
  			  NCVliquid_log,
			  "NCVliquid_phys",
  			  NCVvessel_log,          // mother
  			  false,
  			  0);
  			  
  // Make liquid scintillator a sensitive detector
  //----------------------------------------------
  // Get pointer to detector manager
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  // Create a new instance of NCV sensitive detector
  G4VSensitiveDetector* ncvSD = new NCVSD("NeutronCaptureVolume"); 
  // Register detector with manager
  SDman->AddNewDetector(ncvSD);
  // Attach detector to liquid scintillator volume
  NCVliquid_log->SetSensitiveDetector(ncvSD);

}
