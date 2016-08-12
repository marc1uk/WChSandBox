
#ifndef MRDDetectorConstruction_h
#define MRDDetectorConstruction_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4AssemblyVolume.hh"
	
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ComputePaddleTransformation (const G4int copyNo, G4VPhysicalVolume* physVol);
void ComputeTaperTransformation (const G4int copyNo, G4VPhysicalVolume* physVol, G4bool lgs);
void ComputeSteelTransformation (const G4int copyNo, G4VPhysicalVolume* physVol);
void ComputeVetoPaddleTransformation (const G4int copyNo, G4VPhysicalVolume* physVol, G4bool surf);

void PlacePaddles(G4LogicalVolume* totMRD_log);
void PlaceTapers(G4LogicalVolume* totMRD_log);
void PlaceLGs(G4LogicalVolume* totMRD_log);
void PlaceSteels(G4LogicalVolume* totMRD_log);
void makeAlu(G4AssemblyVolume* totMRD);
void PlaceVetoPaddles(G4LogicalVolume* totVeto_log);



#endif
