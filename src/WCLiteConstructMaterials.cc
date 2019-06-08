#include "WCLiteDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4NistManager.hh" // marcus: added
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

void WCLiteDetectorConstruction::ConstructMaterials()
{
  //****Materials Definitions****

  G4double density;
  G4double a;
  G4int z,n;
  G4int nelements;
  G4double fracMass;

  //---Vacuum

  density     = universe_mean_density;              //from PhysicalConstants.h
  a = 1.01*g/mole;
//  G4double pressure    = 1.e-19*pascal;
//  G4double temperature = 0.1*kelvin;
//  G4Material* Vacuum = new G4Material("Vacuum", 1., a, density, kStateGas,temperature,pressure);

  //LS-Water 
  G4Element* C = new G4Element("Carbon","C", z=6, a=12.01*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* S = new G4Element("Sulphur"  , "S", z=16 , a=32.064*g/mole);
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
 
  G4Material* Scintillator = new G4Material("Scintillator", density=0.988684*g/cm3, nelements=5);
  Scintillator->AddElement(C, fracMass=0.031042);
  Scintillator->AddElement(H, fracMass=0.659);
  Scintillator->AddElement(O, fracMass=0.309);
  Scintillator->AddElement(S, fracMass=0.0009);
  Scintillator->AddElement(N, fracMass=0.000058);
  
  //---Water
  
  a = 1.01*g/mole;
  G4Element* elH 
    = new G4Element("Hydrogen","H", 1,a);
  
  a = 16.00*g/mole;
  G4Element* elO 
    = new G4Element("Oxygen","O", 8,a);
  
  density = 1.00*g/cm3;
  G4Material* Water 
    = new G4Material("Water",density,2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);
  
   //---Steel
  
  a= 12.01*g/mole;
  G4Element* elC 
    = new G4Element("Carbon","C", 6,a);
  
  a = 55.85*g/mole;
  G4Element* elFe
    = new G4Element("Iron","Fe", 26,a);
  
  density = 7.8*g/cm3;
  Steel = new G4Material("Steel",density,2);
  Steel->AddElement(elC, 1.*perCent);
  Steel->AddElement(elFe, 99.*perCent);
  
  //---Stainless steel 304L (simple example, particular alloy can be different)
  
  a = 51.9961*g/mole;
  G4Element* elCr = new G4Element("Chromium", "Cr", 24., a);
  
  a = 58.6934*g/mole;
  G4Element* elNi = new G4Element("Nickel", "Ni", 28., a);

  a = 54.938*g/mole;
  G4Element* elMn = new G4Element("Manganese", "Mn", 25., a); 

  a = 30.974*g/mole;
  G4Element* elP = new G4Element("Phosphore", "P", 15., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element("Silicon", "Si", 14., a); 

  a = 32.065*g/mole;
  G4Element* elS = new G4Element("Sulphur", "S", 16., a);

  density = 7.81*g/cm3;
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density, 8);
 
  StainlessSteel->AddElement(elFe, 70.44*perCent);
  StainlessSteel->AddElement(elCr, 18*perCent);
  StainlessSteel->AddElement(elC,  0.08*perCent);
  StainlessSteel->AddElement(elNi, 8*perCent); 
  StainlessSteel->AddElement(elP, 0.45*perCent);
  StainlessSteel->AddElement(elSi,  1*perCent);
  StainlessSteel->AddElement(elMn, 2*perCent);
  StainlessSteel->AddElement(elS, 0.03*perCent);
  
  G4MaterialPropertiesTable *mpt = new G4MaterialPropertiesTable();
   
  const G4int nEntries = 2;
  G4double photonEnergy[nEntries] = {1.*eV , 7.*eV};
  
  //G4double rindex_Steel[nEntries] = {1.462 , 1.462}; // No I haven't gone mad
  G4double abslength_Steel[nEntries] = {.001*mm , .001*mm};
  //mpt->AddProperty("RINDEX", photonEnergy, rindex_Steel, nEntries);
  mpt->AddProperty("ABSLENGTH", photonEnergy, abslength_Steel, nEntries);
  
  StainlessSteel->SetMaterialPropertiesTable(mpt);
  
// **** NCV additions *****
  // Definition of gadolinium
  // natural gadolinium
  G4Isotope* Gd152 = new G4Isotope("Gd152", z= 64, n= 152, a= 152.0*g/mole);
  G4Isotope* Gd154 = new G4Isotope("Gd154", z= 64, n= 154, a= 154.0*g/mole);
  G4Isotope* Gd155 = new G4Isotope("Gd155", z= 64, n= 155, a= 155.0*g/mole);
  G4Isotope* Gd156 = new G4Isotope("Gd156", z= 64, n= 156, a= 156.0*g/mole);
  G4Isotope* Gd157 = new G4Isotope("Gd157", z= 64, n= 157, a= 157.0*g/mole);
  G4Isotope* Gd158 = new G4Isotope("Gd158", z= 64, n= 158, a= 158.0*g/mole);
  G4Isotope* Gd160 = new G4Isotope("Gd160", z= 64, n= 160, a= 160.0*g/mole);

  G4Element* Gd = new G4Element("Gadolinium", "Gd",7);
  Gd->AddIsotope(Gd152,  0.2*perCent); // beware : it is abundance, not fraction mass
  Gd->AddIsotope(Gd154,  2.2*perCent);
  Gd->AddIsotope(Gd155, 14.9*perCent);
  Gd->AddIsotope(Gd156, 20.6*perCent);
  Gd->AddIsotope(Gd157, 15.7*perCent);
  Gd->AddIsotope(Gd158, 24.7*perCent);
  Gd->AddIsotope(Gd160, 21.7*perCent);
  
  // --- Mineral Oil  (CH2)n ------
  density = 0.77*g/cm3;
  nelements = 2;
  G4Material* mineral_oil 
   = new G4Material("MineralOil", density, nelements);
  mineral_oil -> AddElement(C, 1);
  mineral_oil -> AddElement(H, 2);
  
// --- Pseudo-cumene (C9 H12) also called 1,2,4-Trimethybenzene
  density = 0.8758 *g/cm3;  // at T=20 deg C
  nelements=2;
  G4Material* pseudocumene 
    = new G4Material("pseudocumene",density,nelements);
  pseudocumene->AddElement(C, 9);
  pseudocumene->AddElement(H, 12); 
   
// --- PPO (C15 H11 N 0) -- also called DPO, 2,5-diphenyloxazole
  density = 1.0 *g/cm3;  // ??? at T=?
  G4Material* PPO 
     = new G4Material("PPO",density,nelements=4);
  PPO->AddElement(C, 15);
  PPO->AddElement(H, 11);
  PPO->AddElement(N, 1);
  PPO->AddElement(O, 1);
   
// --- Bis-MSB ------
  density= 1.3*g/cm3;
  G4Material* BisMSB 
    = new G4Material("Bis-MSB",density, nelements= 2);  // density unknown
  BisMSB -> AddElement(C, 24);
  BisMSB -> AddElement(H, 2);
  
// --- NCVliquid - EJ-335 (0.25% Gd), composition taken from the simulation of the Nucifer experiment ---
// --- Concentrations (IN VOLUME): 60% mineral oil, 40% pseudocumene, 6e-3 g/cm3 of PPO, 2e-5 g/cm3 of bis-MSB
// --- Mail from Eljen -> 45% oil, 45% PC, the rest is unknown 
// --- These volume concentrations need to be multiplied by the density to be handled as mass concentrations 
  density = 0.89*g/cm3;
  G4double PPO_fraction= 6.*g/(1e3*cm3*density);
  G4double BisMSB_fraction= 2.*g/(1e5*cm3*density);
  G4double Gd_fraction= 0.25*perCent;
  G4Material* NCVliquid = new G4Material("NCVliquid",density,5);
  NCVliquid->AddMaterial(mineral_oil, 60.*perCent / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) ); G4cout << 60.*perCent / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) << G4endl;
  NCVliquid->AddMaterial(pseudocumene, 40.*perCent / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) ); G4cout << 40.*perCent / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) << G4endl;
  NCVliquid->AddMaterial(PPO, PPO_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) ); G4cout << PPO_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) << G4endl;
  NCVliquid->AddMaterial(BisMSB, BisMSB_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) ); G4cout << BisMSB_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) << G4endl;
  NCVliquid->AddElement(Gd, Gd_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) ); G4cout << Gd_fraction / (1.0 + PPO_fraction + BisMSB_fraction + Gd_fraction) << G4endl;
  
   //               H H 
  // --- Acrylic  -C-C- --------------------
  //               H COOCH3
  density = 1.14*g/cm3;
  nelements = 3;
  G4Material* Acrylic = new G4Material("Acrylic", density, nelements);
  Acrylic -> AddElement(elH, 6);
  Acrylic -> AddElement(elC, 4);
  Acrylic -> AddElement(elO, 2);
 
// **** END NCV additions ****

// **** REPLACED BY NCV additions ****
// a = 157.25*g/mole;
// G4Element* Gd = new G4Element("Gadolinium","Gd", 64,a);
// **** END REPLACED BY NCV additions ****
    
    //---Gd doped Water

  density = 1.00*g/cm3;
  G4Material* DopedWater
     = new G4Material("Doped Water",density,2);
  DopedWater->AddMaterial(Water,99.9*perCent);
  DopedWater->AddElement(Gd,0.1*perCent);


//---Gd doped Scintillator

  density = 1.00*g/cm3;
  G4Material* DopedScintillator
     = new G4Material("Doped Scintillator",density,2);
  DopedScintillator->AddMaterial(Scintillator,99.5*perCent);
  DopedScintillator->AddElement(Gd,0.5*perCent);



  //---Ice 
  
  density = 0.92*g/cm3;//Ice
  G4Material* Ice = new G4Material("Ice",density,2);
  Ice->AddElement(elH, 2);
  Ice->AddElement(elO, 1);

  //---Solid Dry Ice

  density = 1.563*g/cm3;
  G4Material* DryIce = new G4Material("SolidDryIce", density, 2);
  DryIce->AddElement(elC, 1);
  DryIce->AddElement(elO, 2);

  //---Air
  
  a = 14.01*g/mole;
  G4Element* elN 
    = new G4Element("Nitrogen","N", 7,a);
  
  density = 1.290*mg/cm3;
  Air = new G4Material("Air",density,2);
  Air->AddElement(elN, 70.*perCent);
  Air->AddElement(elO, 30.*perCent);
  
  //---Plastic
  
  density = 0.95*g/cm3;
  G4Material* Plastic
    = new G4Material("Plastic",density,2);
  Plastic->AddElement(elC, 1);
  Plastic->AddElement(elH, 2);

  //---Aluminum (element and material)

  a = 26.98*g/mole;
  G4Element* elAl = new G4Element("Aluminum", "Al", 13, a);  

  density = 2.7*g/cm3;
  G4Material* Aluminum
    = new G4Material("Aluminum",density,1);
  Aluminum->AddElement(elAl, 1);

  //---Black sheet

  density = 0.95*g/cm3;
  G4Material* Blacksheet
    = new G4Material("Blacksheet",density,2);
  Blacksheet->AddElement(elC, 1);
  Blacksheet->AddElement(elH, 2);

  //---Glass
 
  density = 2.20*g/cm3;
  G4Material* SiO2 = new G4Material("SiO2",density,2);
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO , 2);

  a = 10.81*g/mole;
  G4Element* elB = new G4Element("Boron", "B", 5, a);  

  density = 2.46*g/cm3;
  G4Material* B2O3 = new G4Material("B2O3",density,2);
  B2O3->AddElement(elB, 2);
  B2O3->AddElement(elO, 3);

  a = 22.99*g/mole;
  G4Element* elNa = new G4Element("Sodium", "Na", 11, a);  

  density = 2.27*g/cm3;
  G4Material* Na2O = new G4Material("Na2O",density,2);
  Na2O->AddElement(elNa, 2);
  Na2O->AddElement(elO, 1);

  density = 4.00*g/cm3;
  G4Material* Al2O3 = new G4Material("Al2O3",density,2);
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);

//   G4Material* blackAcryl
//     = new G4Material("blackAcryl", density, 3);
//   blackAcryl -> AddElement(elH, 6);
//   blackAcryl -> AddElement(elC, 4);
//   blackAcryl -> AddElement(elO, 2);

  density = 2.23*g/cm3;
  G4Material* Glass
    = new G4Material("Glass",density,4);
  //G4Material* Glass
  //= new G4Material("Glass",density,8);  //Put in 8 materials later
  
  Glass->AddMaterial(SiO2, 80.6*perCent);
  Glass->AddMaterial(B2O3, 13.0*perCent);
  Glass->AddMaterial(Na2O, 4.0*perCent);
  Glass->AddMaterial(Al2O3, 2.4*perCent);
  //Glass->AddMaterial(Al2O3, 2.3*perCent);  
  //Put in 2.3 percent if the other 4 materials = 0.1 percent
  
  
  //*****************SciBooNE integration
  
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  
  // Polystyrene scintillator for SciBar (scintillator, WLS fiber) and MRD
  //G4Material* Scinti = man->FindOrBuildMaterial("G4_POLYSTYRENE");	// not found?
  
  density= 1.021 *g/cm3;
  Scinti = new G4Material("Scinti", density, 2);
  Scinti-> AddElement(elC, 8);	// marcus: elC-> C, natoms=8 -> 8
  Scinti-> AddElement(elH, 8);	// marcus: elH-> H, natoms=8 -> 8
  
  
  // Iron for MRD
  a= 55.847 *g/mole;
  density= 7.841 *g/cm3;
  MRDIron = new G4Material("Iron", 26., a, density);
  
  // Aluminum for MRD support structure
  Al  = man->FindOrBuildMaterial("G4_Al");
  
  // Mylar for MRD scintillator paddle cladding
  // from http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt
  G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");
  /*
  density = 1.397*g/cm3;
  G4Material* Mylar = new G4Material("Mylar", density, 3);
  Mylar->AddElement(elO,2);
  Mylar->AddElement(elC,5);
  Mylar->AddElement(elH,4);
  */
  
  G4MaterialPropertiesTable* mptMylar = new G4MaterialPropertiesTable();
  G4double AbsorptionLengthMylar = 100.*cm;
  G4double refractiveIndexMylar = 1.640;
  mptMylar->AddConstProperty("RINDEX",refractiveIndexMylar);
  mptMylar->AddConstProperty("ABSLENGTH",AbsorptionLengthMylar);
  Mylar->SetMaterialPropertiesTable(mptMylar);
  
  //****************/SciBooNE integration
  
  // Sciboone code does not use materials property tables to define scintillation properties??
  // these taken from: http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2013/06/14/11.28-71263-XeDetectorConstruction.cc
  // or http://g4course2010.web.cern.ch/g4course2010/task6/DetectorConstruction_8cc-source.html
  const G4int WLS_NUMENTRIES = 4;
  G4double WLS_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV}; // << WLS energy???
  G4double RIndexPstyrene[WLS_NUMENTRIES]={ 1.5, 1.5, 1.5, 1.5};
  //G4double Absorption1[WLS_NUMENTRIES]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  G4double Absorption1[WLS_NUMENTRIES]={124.*cm,124.*cm,124.*cm,124.*cm};	// average abs length from mrdmodule.txt
  G4double ScintilFast[WLS_NUMENTRIES]={0.00, 0.00, 1.00, 1.00};
  G4MaterialPropertiesTable* MPTPStyrene = new G4MaterialPropertiesTable();
  MPTPStyrene->AddProperty("RINDEX",WLS_Energy,RIndexPstyrene,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("ABSLENGTH",WLS_Energy,Absorption1,WLS_NUMENTRIES);
  MPTPStyrene->AddProperty("FASTCOMPONENT",WLS_Energy, ScintilFast, WLS_NUMENTRIES);
  // also available: slow component. not needed? Organic scintillators tend to have fast response. Need to also specify YIELDRATIO if
  // specifying both FAST and SLOW components. 
  MPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  MPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  Scinti->SetMaterialPropertiesTable(MPTPStyrene);
  // Set the Birks Constant for the MRD scintillator layers
  Scinti->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  
  /* if scintillation depends on particle type, see http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch05s02.html#sect.PhysProc.Photo */
  
  /* Do we need to add a wavelength shifter (WLS)? see above link if so.
     Do we need to specify scattering - RAYLEIGH and/or MIE? see further down in this doc for templates. Some numbers entered 
     but not used, described as "utter fiction" - more suitable values? */

// --------------------------------------------------------------
// Materials for FACC scintillation paddles including WLS fibres
// taken from /examples/extended/optical/help, but properties have been verified as correct for CDF paddles
// --------------------------------------------------------------


//  //--------------------------------------------------
//  // WLSfiber PMMA
//  //--------------------------------------------------

//  density = 1.190*g/cm3;
//  std::vector<G4int> natoms;
//  std::vector<G4double> fractionMass;
//  std::vector<G4String> elements;

//  elements.push_back("C");     natoms.push_back(5);
//  elements.push_back("H");     natoms.push_back(8);
//  elements.push_back("O");     natoms.push_back(2);

//  G4Material* PMMA = fNistMan-> ConstructNewMaterial("PMMA", elements, natoms, density);

//  elements.clear();
//  natoms.clear();

//  //--------------------------------------------------
//  // Cladding (polyethylene)
//  //--------------------------------------------------

//  elements.push_back("C");     natoms.push_back(2);
//  elements.push_back("H");     natoms.push_back(4);

//  density = 1.200*g/cm3;

//  fPethylene = fNistMan->
//          ConstructNewMaterial("Pethylene", elements, natoms, density);

//  elements.clear();
//  natoms.clear();

//  //--------------------------------------------------
//  // Double Cladding (fluorinated polyethylene)
//  //--------------------------------------------------

//  elements.push_back("C");     natoms.push_back(2);
//  elements.push_back("H");     natoms.push_back(4);

//  density = 1.400*g/cm3;

//  fFPethylene = fNistMan->
//          ConstructNewMaterial("FPethylene", elements, natoms, density);

//  elements.clear();
//  natoms.clear();

//  //--------------------------------------------------
//  // Polystyrene
//  //--------------------------------------------------
// 
//  elements.push_back("C");     natoms.push_back(8);
//  elements.push_back("H");     natoms.push_back(8);

//  density = 1.050*g/cm3;

//  fPolystyrene = fNistMan->
//          ConstructNewMaterial("Polystyrene", elements, natoms, density);

//  elements.clear();
//  natoms.clear();


// -------------------------------------------------------------
// Generate & Add Material Properties Table
// -------------------------------------------------------------

/*
  const G4int NUMENTRIES = 32;
  
  G4double PPCKOV[NUMENTRIES] =
    { 2.034E-9*GeV, 2.068E-9*GeV, 2.103E-9*GeV, 2.139E-9*GeV,
      2.177E-9*GeV, 2.216E-9*GeV, 2.256E-9*GeV, 2.298E-9*GeV,
      2.341E-9*GeV, 2.386E-9*GeV, 2.433E-9*GeV, 2.481E-9*GeV,
      2.532E-9*GeV, 2.585E-9*GeV, 2.640E-9*GeV, 2.697E-9*GeV,
      2.757E-9*GeV, 2.820E-9*GeV, 2.885E-9*GeV, 2.954E-9*GeV,
      3.026E-9*GeV, 3.102E-9*GeV, 3.181E-9*GeV, 3.265E-9*GeV,
      3.353E-9*GeV, 3.446E-9*GeV, 3.545E-9*GeV, 3.649E-9*GeV,
      3.760E-9*GeV, 3.877E-9*GeV, 4.002E-9*GeV, 4.136E-9*GeV };
*/


  // default values
  /*
  G4double RINDEX1[NUMENTRIES] =
    { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
      1.346,  1.3465, 1.347,  1.3475, 1.348,
      1.3485, 1.3492, 1.35,   1.3505, 1.351,
      1.3518, 1.3522, 1.3530, 1.3535, 1.354,
      1.3545, 1.355,  1.3555, 1.356,  1.3568,
      1.3572, 1.358,  1.3585, 1.359,  1.3595,
      1.36,   1.3608};
  */
  
  // M Fechner,  values from refsg.F in apdetsim
  /*  G4double RINDEX1[NUMENTRIES] = 
    { 1.33332, 1.333364, 1.333396, 1.3343, 1.33465,
      1.33502, 1.3354, 1.33579, 1.3362, 1.33663, 1.33709,
      1.33756, 1.33806, 1.3386, 1.33915, 1.33974,
      1.34037, 1.34105, 1.34176, 1.34253, 1.34336,
      1.34425, 1.34521, 1.34626, 1.3474, 1.34864,
      1.35002, 1.35153, 1.35321, 1.35507, 1.35717, 1.35955 };
  */

   //From SFDETSIM water absorption
   const G4int NUMENTRIES_water=60;

   G4double ENERGY_water[NUMENTRIES_water] =
     { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
       1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
       1.74647e-09*GeV, 1.77142e-09*GeV,1.7971e-09*GeV, 1.82352e-09*GeV, 
       1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
       1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
       2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
       2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
       2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
       2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
       2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
       3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
       3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
       3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
       4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
       5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };



      // Air
   G4double RINDEX_air[NUMENTRIES_water] = 
     { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
       

   // M Fechner : new ; define the water refraction index using refsg.F 
   //from skdetsim using the whole range.   
    G4double RINDEX1[NUMENTRIES_water] = 
     {1.32885, 1.32906, 1.32927, 1.32948, 1.3297, 1.32992, 1.33014, 
      1.33037, 1.3306, 1.33084, 1.33109, 1.33134, 1.3316, 1.33186, 1.33213,
      1.33241, 1.3327, 1.33299, 1.33329, 1.33361, 1.33393, 1.33427, 1.33462,
      1.33498, 1.33536, 1.33576, 1.33617, 1.3366, 1.33705, 1.33753, 1.33803,
      1.33855, 1.33911, 1.3397, 1.34033, 1.341, 1.34172, 1.34248, 1.34331,
      1.34419, 1.34515, 1.3462, 1.34733, 1.34858, 1.34994, 1.35145, 1.35312,
      1.35498, 1.35707, 1.35943, 1.36211, 1.36518, 1.36872, 1.37287, 1.37776,
      1.38362, 1.39074, 1.39956, 1.41075, 1.42535};
   

   
   G4double ABSORPTION_water[NUMENTRIES_water] = //old
     { 22.8154*cm, 28.6144*cm, 35.9923*cm, 45.4086*cm, 57.4650*cm,
       72.9526*cm, 92.9156*cm, 118.737*cm, 152.255*cm, 195.925*cm,
       202.429*cm, 224.719*cm, 236.407*cm, 245.700*cm, 289.017*cm,
       305.810*cm, 316.456*cm, 326.797*cm, 347.222*cm, 414.938*cm,
       636.943*cm, 934.579*cm, 1245.33*cm, 1402.52*cm, 1550.39*cm,
       1745.20*cm, 1883.24*cm, 2016.13*cm, 2442.18*cm, 3831.28*cm,
       4652.89*cm, 5577.04*cm, 6567.08*cm, 7559.88*cm, 8470.06*cm,
       9205.54*cm, 9690.95*cm, 9888.36*cm, 9804.50*cm, 9482.17*cm,
       8982.77*cm, 8369.39*cm, 7680.31*cm, 6902.11*cm, 6183.84*cm,
       5522.27*cm, 4914.33*cm, 4357.09*cm, 3847.72*cm, 3383.51*cm,
       2961.81*cm, 2580.08*cm, 2235.83*cm, 1926.66*cm, 1650.21*cm,
       1404.21*cm, 1186.44*cm, 994.742*cm, 827.027*cm, 681.278*cm};
    /*
   G4double ABSORPTION_water[NUMENTRIES_water] = //new
     {25.3504*cm, 31.7938*cm, 39.9915*cm, 50.454*cm, 63.85*cm, 
      81.0584*cm, 103.24*cm, 131.93*cm, 169.172*cm, 217.694*cm, 
      224.921*cm, 249.688*cm, 262.674*cm, 273*cm, 321.13*cm, 339.789*cm,
      351.617*cm, 363.108*cm, 385.802*cm, 461.042*cm, 707.714*cm, 
      1038.42*cm, 1383.7*cm, 1558.36*cm, 1722.65*cm, 1939.11*cm, 
      2092.49*cm, 2240.14*cm, 2962.96*cm, 4967.03*cm, 6368.58*cm, 
      8207.56*cm, 10634.2*cm, 13855.3*cm, 18157.3*cm, 23940.2*cm, 
      31766*cm, 42431.6*cm, 57074.9*cm, 77335.9*cm, 105598*cm, 
      145361*cm, 192434*cm, 183898*cm, 176087*cm, 168913*cm, 162301*cm, 
      156187*cm, 150516*cm, 145243*cm, 140327*cm, 135733*cm, 131430*cm, 
      127392*cm, 123594*cm, 120016*cm, 116640*cm, 113448*cm, 110426*cm, 
      107562*cm};
  */
   // M Fechner: Rayleigh scattering -- as of version 4.6.2 of GEANT,
   // one may use one's own Rayleigh scattering lengths (the buffer is no
   // longer overwritten for "water", see 4.6.2 release notes)

   // RAYFF = 1/ARAS, for those of you who know SKdetsim...
   // actually that's not quite right because the scattering models
   // are different; in G4 there is no scattering depolarization
   // std value at SK = 0.6. But Mie scattering is implemented
   // in SKdetsim and not in G4

   
  // april 2005 : reduced reflections, let's increase scattering...
  // sep 09: for the large detector like superK the old values are muc better
      G4double RAYFF = 1.0/1.65;  //old
   // G4double RAYFF = 1.0/1.5;  new

   G4double RAYLEIGH_water[NUMENTRIES_water] = {
     167024.4*cm*RAYFF, 158726.7*cm*RAYFF, 150742*cm*RAYFF,
     143062.5*cm*RAYFF, 135680.2*cm*RAYFF, 128587.4*cm*RAYFF,
     121776.3*cm*RAYFF, 115239.5*cm*RAYFF, 108969.5*cm*RAYFF,
     102958.8*cm*RAYFF, 97200.35*cm*RAYFF, 91686.86*cm*RAYFF,
     86411.33*cm*RAYFF, 81366.79*cm*RAYFF, 76546.42*cm*RAYFF,
     71943.46*cm*RAYFF, 67551.29*cm*RAYFF, 63363.36*cm*RAYFF,
     59373.25*cm*RAYFF, 55574.61*cm*RAYFF, 51961.24*cm*RAYFF,
     48527.00*cm*RAYFF, 45265.87*cm*RAYFF, 42171.94*cm*RAYFF,
     39239.39*cm*RAYFF, 36462.50*cm*RAYFF, 33835.68*cm*RAYFF,
     31353.41*cm*RAYFF, 29010.30*cm*RAYFF, 26801.03*cm*RAYFF,
     24720.42*cm*RAYFF, 22763.36*cm*RAYFF, 20924.88*cm*RAYFF,
     19200.07*cm*RAYFF, 17584.16*cm*RAYFF, 16072.45*cm*RAYFF,
     14660.38*cm*RAYFF, 13343.46*cm*RAYFF, 12117.33*cm*RAYFF,
     10977.70*cm*RAYFF, 9920.416*cm*RAYFF, 8941.407*cm*RAYFF,
     8036.711*cm*RAYFF, 7202.470*cm*RAYFF, 6434.927*cm*RAYFF,
     5730.429*cm*RAYFF, 5085.425*cm*RAYFF, 4496.467*cm*RAYFF,
     3960.210*cm*RAYFF, 3473.413*cm*RAYFF, 3032.937*cm*RAYFF,
     2635.746*cm*RAYFF, 2278.907*cm*RAYFF, 1959.588*cm*RAYFF,
     1675.064*cm*RAYFF, 1422.710*cm*RAYFF, 1200.004*cm*RAYFF,
     1004.528*cm*RAYFF, 833.9666*cm*RAYFF, 686.1063*cm*RAYFF
   };

   //From SFDETSIM
   /*
   G4double RINDEX_glass[NUMENTRIES] =
     { 1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600 };
   */
   // M Fechner : unphysical, I want to reduce reflections
   // back to the old value 1.55

   G4double RINDEX_glass[NUMENTRIES_water] =
     { 1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600,
       1.600, 1.600 }; 

   //G4double RINDEX_blacksheet[NUMENTRIES] =
   //{ 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500,
   //  2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500,
   //  2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500,
   //  2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500,
   //  2.500, 2.500, 2.500, 2.500 };
   
   
   //G4double ABSORPTION1[NUMENTRIES] =
   //{344.8*cm,  408.2*cm,  632.9*cm,  917.4*cm, 1234.6*cm, 1388.9*cm,
   // 1515.2*cm, 1724.1*cm, 1886.8*cm, 2000.0*cm, 2631.6*cm, 3571.4*cm,
   // 4545.5*cm, 4761.9*cm, 5263.2*cm, 5263.2*cm, 5555.6*cm, 5263.2*cm,
   // 5263.2*cm, 4761.9*cm, 4545.5*cm, 4166.7*cm, 3703.7*cm, 3333.3*cm,
   // 3000.0*cm, 2850.0*cm, 2700.0*cm, 2450.0*cm, 2200.0*cm, 1950.0*cm,
   // 1750.0*cm, 1450.0*cm };
   
   /*   
   G4double ABSORPTION_glass[NUMENTRIES] = 
     { 100.0*cm, 110.0*cm, 120.0*cm, 130.0*cm, 140.0*cm, 150.0*cm, 160.0*cm,
       165.0*cm, 170.0*cm, 175.0*cm, 180.0*cm, 185.0*cm, 190.0*cm, 195.0*cm,
       200.0*cm, 200.0*cm, 200.0*cm, 200.0*cm, 200.0*cm, 195.0*cm, 190.0*cm,
       185.0*cm, 180.0*cm, 175.0*cm, 170.0*cm, 160.0*cm, 150.0*cm, 140.0*cm,
       130.0*cm, 120.0*cm, 110.0*cm, 100.0*cm };
   */
   // M Fechner : the quantum efficiency already takes glass abs into account

   G4double ABSORPTION_glass[NUMENTRIES_water]= 
     { 1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm, 1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,1.0e9*cm,
       1.0e9*cm, 1.0e9*cm };
   
   G4double BLACKABS_blacksheet[NUMENTRIES_water] =
     { 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 
       1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm,
       1.0e-9*cm, 1.0e-9*cm};
   
   
   //The following reflectivity for blacksheet is obtained from skdetsim
   //There is a fudge factor of 2.7 for the reflectivity of blacksheet
   //depending on whether SK1 or SK2 simulation is used.  
   //The BlackSheetFudgeFactor is set to true if you want to use the 
   //SK2 values, false if not.
   G4double SK1SK2FF = 1.0;
   //G4bool BlackSheetFudgeFactor=false;
   G4bool BlackSheetFudgeFactor=true;
   //   if (BlackSheetFudgeFactor) SK1SK2FF=SK1SK2FF*2.7;
   if (BlackSheetFudgeFactor) SK1SK2FF=SK1SK2FF*1.55;

   /*
   G4double REFLECTIVITY_blacksheet[NUMENTRIES] =
     { 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.057*SK1SK2FF, 0.059*SK1SK2FF, 0.060*SK1SK2FF, 
       0.059*SK1SK2FF, 0.058*SK1SK2FF, 0.057*SK1SK2FF, 0.055*SK1SK2FF, 
       0.050*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF };
   */
   
   G4double REFLECTIVITY_blacksheet[NUMENTRIES_water] =
     { 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
       0.055*SK1SK2FF, 0.057*SK1SK2FF, 0.059*SK1SK2FF, 0.060*SK1SK2FF, 
       0.059*SK1SK2FF, 0.058*SK1SK2FF, 0.057*SK1SK2FF, 0.055*SK1SK2FF, 
       0.050*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
       0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF ,
       0.045*SK1SK2FF, 0.045*SK1SK2FF };


   //utter fiction at this stage
   G4double EFFICIENCY[NUMENTRIES_water] =
     { 0.001*m };
      
   //utter fiction at this stage, does not matter
   G4double RAYLEIGH_air[NUMENTRIES_water] =
     { 0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,
       0.001*m,0.001*m,0.001*m,0.001*m,0.001*m,0.001*m};
      
   //Not used yet, fictional values
   //G4double SPECULARLOBECONSTANT1[NUMENTRIES] =
   //{ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001 };
   
   //Not used yet, fictional values
   //G4double SPECULARSPIKECONSTANT1[NUMENTRIES] =
   //{ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001 };
   
   //Not used yet, fictional values
   //G4double BACKSCATTERCONSTANT1[NUMENTRIES] =
   //{ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
   //  0.001, 0.001, 0.001, 0.001 };
   
//   G4double EFFICIENCY_blacksheet[NUMENTRIES_water] =
//    { 0.0 };

   //	------------- Surfaces --------------

   /*
   OpWaterBSSurface =
     new G4OpticalSurface("WaterBSCellSurface");

   OpWaterBSSurface->SetType(dielectric_dielectric);
   OpWaterBSSurface->SetModel(unified);
   OpWaterBSSurface->SetFinish(groundfrontpainted);
   OpWaterBSSurface->SetSigmaAlpha(0.1);

   const G4int NUM = 2;
   //   G4double PP[NUM] =
   //{ 2.038E-9*GeV, 4.144E-9*GeV };

   G4double PP[NUM] = { 1.4E-9*GeV,6.2E-9*GeV};
   G4double RINDEX_blacksheet[NUM] =
     { 1.6, 1.6 };

   G4double SPECULARLOBECONSTANT[NUM] =
     { 0.3, 0.3 };
   G4double SPECULARSPIKECONSTANT[NUM] =
     { 0.2, 0.2 };
   G4double BACKSCATTERCONSTANT[NUM] =
     { 0.2, 0.2 };

   OpGlassCathodeSurface =
     new G4OpticalSurface("GlassCathodeSurface");
   OpGlassCathodeSurface->SetType(dielectric_dielectric);
   OpGlassCathodeSurface->SetModel(unified);
   //   OpGlassCathodeSurface->SetFinish(groundbackpainted);
   OpGlassCathodeSurface->SetFinish(polished);
   //OpGlassCathodeSurface->SetSigmaAlpha(0.002);
   // was 1.0
   // totally unphysical anyway 
   G4double RINDEX_cathode[NUM] =
     { 1.0, 1.0 };
      
   G4double SPECULARLOBECONSTANT_glasscath[NUM] =
     { 1.0, 1.0 };
     // { 0.3, 0.3 };
   G4double SPECULARSPIKECONSTANT_glasscath[NUM] =
     { 0.0, 0.0 };
     //     { 0.2, 0.2 };
   G4double BACKSCATTERCONSTANT_glasscath[NUM] =
     {0.0, 0.0};
   //     { 0.2, 0.2 };
   
   G4double REFLECTIVITY_glasscath[NUM] =
     { 0.0, 0.0 };
   G4double EFFICIENCY_glasscath[NUM] =
     { 0.0, 0.0 };
   */

   G4MaterialPropertiesTable *myMPT1 = new G4MaterialPropertiesTable();
   // M Fechner : new   ; wider range for lambda
   myMPT1->AddProperty("RINDEX", ENERGY_water, RINDEX1, NUMENTRIES_water);
   myMPT1->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
//   myMPT1->AddConstProperty("SCINTILLATIONYIELD",0./MeV);
//   myMPT1->AddConstProperty("RESOLUTIONSCALE",.0);
   // M Fechner: new, don't let G4 compute it.
   myMPT1->AddProperty("RAYLEIGH",ENERGY_water,RAYLEIGH_water,NUMENTRIES_water);
   Water->SetMaterialPropertiesTable(myMPT1);
   //Gd doped water has the same optical properties as pure water
   DopedWater->SetMaterialPropertiesTable(myMPT1);
   //myMPT1->DumpTable();
   



 // DEFINE SCINTILLATOR MATERIAL PROPERTY

   const G4int NumEntries = 60;
   
   G4double Scint_Energy[NumEntries] = {2.65086E-9*GeV,2.70739E-9*GeV,2.76153E-9*GeV,2.80895E-9*GeV,2.81081E-9*GeV,2.84367E-9*GeV,2.86017E-9*GeV,2.8767E-9*GeV,2.8942E-9*GeV,2.89602E-9*GeV,2.91479E-9*GeV,2.94447E-9*GeV,2.95178E-9*GeV,2.95449E-9*GeV,2.96948E-9*GeV,2.98294E-9*GeV,2.99728E-9*GeV,3.00297E-9*GeV,3.01188E-9*GeV,3.01761E-9*GeV,3.02239E-9*GeV,3.06615E-9*GeV,3.07136E-9*GeV,3.07403E-9*GeV,3.07766E-9*GeV,3.08913E-9*GeV,3.09219E-9*GeV,3.1048E-9*GeV,3.14546E-9*GeV,3.15049E-9*GeV,3.17008E-9*GeV,3.17567E-9*GeV,3.19696E-9*GeV,3.20206E-9*GeV,3.20943E-9*GeV,3.21148E-9*GeV,3.21352E-9*GeV,3.2208E-9*GeV,3.2318E-9*GeV,3.23576E-9*GeV,3.24292E-9*GeV,3.25304E-9*GeV,3.28128E-9*GeV,3.28671E-9*GeV,3.2919E-9*GeV,3.29532E-9*GeV,3.30592E-9*GeV,3.31226E-9*GeV,3.31979E-9*GeV,3.32572E-9*GeV,3.3332E-9*GeV,3.33545E-9*GeV,3.3435E-9*GeV,3.36639E-9*GeV,3.41137E-9*GeV,3.41482E-9*GeV,3.45157E-9*GeV,3.45439E-9*GeV,3.46271E-9*GeV,4.66655E-9*GeV};
   
//   G4double RI_Scint[NumEntries] = {1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,1.3492,};
   
   
   G4double Scint_WB_LS[NumEntries] = {.07208,.08919,.12475,.15227,.15268,.17306,.18597,.21360,.22732,.21627,.24965,.28076,.31149,.27979,.32093,.33267,.38079,.37930,.37735,.39199,.41786,.49865,.49669,.49669,.53000,.54034,.56376,.57231,.69088,.71348,.73073,.72491,.77569,.82767,.82767,.82978,.84346,.84021,.87063,.86140,.86035,.89661,.91665,.91665,.95458,.94635,.93879,.94448,.92571,.94167,.96505,.96505,.97168,1.00000,.87935,.86627,.49176,.43960,.39031,.01287};
   
   
   G4MaterialPropertiesTable* myMPT7 = new G4MaterialPropertiesTable();
 

   myMPT7->AddProperty("RINDEX", ENERGY_water, RINDEX1, NUMENTRIES_water);
   myMPT7->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
   myMPT7->AddProperty("RAYLEIGH",ENERGY_water,RAYLEIGH_water,NUMENTRIES_water);
   /*
   myMPT7->AddProperty("RINDEX",Scint_Energy, RI_Scint,NumEntries);
   myMPT7->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_water, NUMENTRIES_water);
   myMPT7->AddProperty("RAYLEIGH",ENERGY_water,RAYLEIGH_water,NUMENTRIES_water);
   */
   myMPT7->AddProperty("FASTCOMPONENT",Scint_Energy, Scint_WB_LS,NumEntries);
   myMPT7->AddProperty("SLOWCOMPONENT",Scint_Energy, Scint_WB_LS,NumEntries);
   //myMPT7->AddProperty("MIEHG",ENERGY_water,MIE_water,NUMENTRIES_water);
   //myMPT7->AddConstProperty("MIEHG_FORWARD",MIE_water_const[0]);
   //myMPT7->AddConstProperty("MIEHG_BACKWARD",MIE_water_const[1]);
   //myMPT7->AddConstProperty("MIEHG_FORWARD_RATIO",MIE_water_const[2]);
   
   myMPT7->AddConstProperty("SCINTILLATIONYIELD",1800./MeV);
   myMPT7->AddConstProperty("RESOLUTIONSCALE",.01);
   myMPT7->AddConstProperty("FASTTIMECONSTANT", 1.23*ns);
   myMPT7->AddConstProperty("SLOWTIMECONSTANT",9.26*ns);
   myMPT7->AddConstProperty("YIELDRATIO",0.26);
   
   Scintillator->SetMaterialPropertiesTable(myMPT7);
   Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
   

// could be useful to have optical properties of air for bubbles
   G4MaterialPropertiesTable *myMPT2 = new G4MaterialPropertiesTable();
   myMPT2->AddProperty("RINDEX", ENERGY_water, RINDEX_air, NUMENTRIES_water);
   // M Fechner : what is that ?????
   myMPT2->AddProperty("ABSLENGTH", ENERGY_water, BLACKABS_blacksheet, NUMENTRIES_water);
   myMPT2->AddProperty("RAYLEIGH",ENERGY_water, RAYLEIGH_air, NUMENTRIES_water);
//   myMPT2->AddConstProperty("SCINTILLATIONYIELD",0./MeV);	 shouldn't be necessary if it's not defined as a scintillator
//   myMPT2->AddConstProperty("RESOLUTIONSCALE",.0);
   Air->SetMaterialPropertiesTable(myMPT2);
   
   G4MaterialPropertiesTable *myMPT3 = new G4MaterialPropertiesTable();
   myMPT3->AddProperty("ABSLENGTH", ENERGY_water, BLACKABS_blacksheet, NUMENTRIES_water);
   myMPT3->AddProperty("REFLECTIVITY", ENERGY_water, REFLECTIVITY_blacksheet, NUMENTRIES_water);
   myMPT3->AddProperty("EFFICIENCY",   ENERGY_water, EFFICIENCY, NUMENTRIES_water);
   Plastic->SetMaterialPropertiesTable(myMPT3);
   
   G4MaterialPropertiesTable *myMPT4 = new G4MaterialPropertiesTable();
   myMPT4->AddProperty("ABSLENGTH", ENERGY_water, BLACKABS_blacksheet, NUMENTRIES_water);
   Blacksheet->SetMaterialPropertiesTable(myMPT4);
   
   G4MaterialPropertiesTable *myMPT5 = new G4MaterialPropertiesTable();
   myMPT5->AddProperty("RINDEX", ENERGY_water, RINDEX_glass, NUMENTRIES_water);
   myMPT5->AddProperty("ABSLENGTH",ENERGY_water, ABSORPTION_glass, NUMENTRIES_water);
   Glass->SetMaterialPropertiesTable(myMPT5);
   
 
 
   //	------------- Surfaces --------------

   // Blacksheet
   /*
   G4MaterialPropertiesTable *myST1 = new G4MaterialPropertiesTable();
   myST1->AddProperty("RINDEX", ENERGY_water, RINDEX_blacksheet, NUMENTRIES_water);
   myST1->AddProperty("SPECULARLOBECONSTANT", PP, SPECULARLOBECONSTANT, NUM);
   myST1->AddProperty("SPECULARSPIKECONSTANT", PP, SPECULARSPIKECONSTANT, NUM);
   myST1->AddProperty("BACKSCATTERCONSTANT", PP, BACKSCATTERCONSTANT, NUM);
   myST1->AddProperty("REFLECTIVITY", ENERGY_water, REFLECTIVITY_blacksheet, NUMENTRIES_water);
   myST1->AddProperty("EFFICIENCY", ENERGY_water, EFFICIENCY_blacksheet, NUMENTRIES_water);
   OpWaterBSSurface->SetMaterialPropertiesTable(myST1);
   */

   //Glass to Cathode surface inside PMTs
   /*
   G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();
   myST2->AddProperty("RINDEX", PP, RINDEX_cathode, NUM);
   //   myST2->AddProperty("SPECULARLOBECONSTANT", PP, SPECULARLOBECONSTANT_glasscath, NUM);
   //   myST2->AddProperty("SPECULARSPIKECONSTANT", PP, SPECULARSPIKECONSTANT_glasscath, NUM);
   //   myST2->AddProperty("BACKSCATTERCONSTANT", PP, BACKSCATTERCONSTANT_glasscath, NUM);
   myST2->AddProperty("REFLECTIVITY", PP, REFLECTIVITY_glasscath, NUM);
   myST2->AddProperty("EFFICIENCY", PP, EFFICIENCY_glasscath, NUM);
   //myST2->AddProperty("ABSLENGTH", PP, abslength_paint , NUM);
   OpGlassCathodeSurface->SetMaterialPropertiesTable(myST2);
   */

}
