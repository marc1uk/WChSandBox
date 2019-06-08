// ====================================================================
//   MRDDetectorConstruction.cc
//
//   17/11/15 M. O'Flaherty
// ====================================================================
//===============================================================================================================
  //MRD DETECTOR DEFINITION
//===============================================================================================================

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
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
#include "MRDDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"

G4double Xposition=0, Yposition=0, Zposition=0;		// used for positioning parameterisations.
G4int numpaddlesperpanelh=26;								// paddles per h scintillator panel
G4int numpaddlesperpanelv=30;								// paddles per v scintillator panel
G4int numpanels=12;													// scintillator panels
G4int numrpcs=10;														// rpc panels
G4int numplates=12;													// steel plates
G4int numalustructs=13;											// number of supporting structs. We may be dropping one as we have fewer scintillators?
G4int numvetopaddles=26;										// number of scintillator paddles in the FACC; 13 panels in 2 layers
G4int vetopaddlesperpanel=13;								// number of scintillator paddles in each FACC panel

G4double scintbordergap=0.3*cm;	// gap between each scintillator (cm) (border to account for cladding etc)
G4double steelscintgap=0.5*cm;	// gap between steel and scintillator
G4double scintalugap=0.2*cm;			// gap between scintillator and alu struct
G4double alusteelgap=2.0*cm; 		// gap between alu struct and subsequent steel of next layer
G4double layergap = steelscintgap + scintalugap + alusteelgap;	// total gaps of layers

G4double steelfullxlen = 305*cm;
G4double steelfullylen = 274*cm;
G4double steelfullzlen = 5*cm;

G4double scintfullxlen = 20*cm;
G4double scintfullzlen= 0.6*cm;
G4double scinthfullylen = 147.2*cm; //155cm - 7.8cm tapered section		147.1 according to sciboone gdml export		//swapped???
G4double scintvfullylen= 130.2*cm;  //138cm - 7.8cm tapered section		129.4 according to sciboone gdml export

G4double scinttapfullwidth = 17.1*cm; 			// width of the tapering part of scintillator paddles at the narrow end
G4double scinttapfullheight = 7.8*cm; 			// z length of tapering part of scint paddles.

G4double scintlgfullwidth = 5.08*cm; 			// tapered light guides at narrow end
G4double scintlgfullheight = 33.3*cm; 			// 

G4double alufullxlen = steelfullxlen+15*cm;	// outer thicknesses - total frame dims are about those of the steel plate
G4double alufullylen = steelfullylen+15*cm;	// basically making this up
G4double alufullzlen = 3.81*cm;							// from eye and a skim of the mrdmodule.txt file, i'm guessing depth is ~0.75 inches (1.9cm)
G4double alufullxthickness = 2.54*cm;				// as above, guessing frame to be 1 inch box cross-section
G4double alufullythickness = 2.54*cm;
G4double windowwidth = (steelfullxlen-(4*alufullxthickness))/3;	// (full length - 4 beams) / 3 windows
G4double windowheight= (steelfullylen-(4*alufullythickness))/3;
	
G4double mrdZlen = numplates*steelfullzlen + (numpanels+1)*scintfullzlen + numalustructs*alufullzlen + numpanels*layergap + scintalugap; 
// add another panel to full length because last alu struct is placed back assuming it's there. Maybe need to change... 

G4double vetopaddlefullxlen = 320.0*cm;
G4double vetopaddlefullylen = 30.5*cm;
G4double vetopaddlefullzlen = 2.0*cm;
G4double vetolayer2offset = 1.27*cm;
G4double vetopaddlegap = 0.2*cm;
G4double nothickness = 0.01*cm;

G4double vetoZlen = 2*vetopaddlefullzlen+vetopaddlegap;

extern G4double tankouterRadius;

// following measurements taken with tape measure; can be used to derive gaps. 
// veto wall full height measured to be 158.75" = 403.2cm / 13 paddles = 31.1cm each in y height
// paddle full x length measured to be 126.5" = 321.3cm in square section, followed by taper
// paddle full z length measured to be 1" = 2.54cm;
// taper length measured to be 14.5" = 36.8cm
// taper narrow width measured to be 2" = 5.1cm (diameter of PMT)
// taper depth at narrow end is difficult to measure due to interface between PMT and paddle/LG - assuming same
// layer offset is 0.5" = 1.27cm

/*
Allowing for a typical deterioration rate of 5–10% per year, full efficiency should be retained more than 10 years and it will be over than useful lifetime of CDF. .. The technique relies on a wavelength shifter fibers to extract the light from the longer side of the scintillator bar. ... the scintillator used to construct the counters (UPS 923A) is a polystyrene-based plastic...The results of quality control tests performed at JINR show that the average ligh output ranges between 21ph.e./MIP (for the longest counters) for muons traversing the counters transversely at the furthest ends from the photomultipliers.
*/

G4double widths[] = {2*(scinthfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)),((numpaddlesperpanelv/2)*(scintfullxlen+scintbordergap))};	
// 2* and y dim because we stack 2 rotated h scint paddles end-to-end. 
//
G4double heights[] = {2*(scintvfullylen+scinttapfullheight+scintlgfullheight+(scintbordergap/2)),((numpaddlesperpanelh/2)*(scintfullxlen+scintbordergap))};
//
G4double maxwidth = *std::max_element(widths,widths+(sizeof(widths)/sizeof(widths[0])))+0.1*cm;
G4double maxheight = *std::max_element(heights,heights+(sizeof(heights)/sizeof(heights[0])))+0.1*cm;
        
// Define solids 
//==============  
// G4Box* variableName = new G4Box("SolidName", x_halflength, y_halflength, z_halflength);
// G4Trd("SolidName", x_halflength1, x_halflength2, y_halflength1, y_halflength2, z_halflength); 
// 2 x and y dimensions given - define dims of the cross-sections at the two z ends. 

// Paddles - h and v
G4Box* sciMRDhpaddle_box = new G4Box("scintHpaddle",scintfullxlen/2,scinthfullylen/2,scintfullzlen/2);
G4Box* sciMRDvpaddle_box = new G4Box("scintVpaddle",scintfullxlen/2,scintvfullylen/2,scintfullzlen/2);

// Paddle Tapered ends
G4Trd* mrdScintTap_box = new G4Trd("mrdScintTap_box", scintfullxlen/2, scinttapfullwidth/2, scintfullzlen/2, scintfullzlen/2, scinttapfullheight/2);	

// Tapered Light Guides - same code as tapered ends. re-use rotations matrices as they're the same.
G4Trd* mrdLG_box = new G4Trd("mrdLG_box", scinttapfullwidth/2, scintlgfullwidth/2, scintfullzlen/2, scintfullzlen/2, scintlgfullheight/2);	

// Steel plates
G4Box* steelMRDplate_box = new G4Box("steelPlate",steelfullxlen/2,steelfullylen/2,steelfullzlen/2);

//The alu support structure is roughly the external size of the steel plates...
G4Box* aluMRDstruc_box = new G4Box("outer_Box", alufullxlen/2, alufullylen/2, alufullzlen/2);

// ...with a 3x3 grid of ~even sized holes  
G4Box* aluMRDwindow_box = new G4Box("inner_Box", windowwidth/2, windowheight/2, alufullzlen/2);

G4Box* totMRD_box = new G4Box("totMRD",(maxwidth/2),(maxheight/2),mrdZlen/2);

G4Box* vetoPaddle_box = new G4Box("vetoPaddle_box",vetopaddlefullxlen/2, vetopaddlefullylen/2, vetopaddlefullzlen/2);
G4Box* vetoSurface_box = new G4Box("vetoSurface_box",vetopaddlefullxlen/2,nothickness,vetopaddlefullzlen/2);

G4Box* totVeto_box = new G4Box("totVeto_box", (vetopaddlefullxlen/2), ((vetopaddlefullylen/2)+vetopaddlegap)*(vetopaddlesperpanel)+(vetolayer2offset/2)+nothickness*2, (vetopaddlefullzlen+(vetopaddlegap/2)));
	
// Define logical volumes 
//=======================
// G4LogicalVolume* variableName = new G4LogicalVolume(solidVariableName, Material, "logicalVolName");
// Can't do this outside a function - the materials aren't defined yet.

G4LogicalVolume* hpaddle_log;
			
G4LogicalVolume* vpaddle_log;

G4LogicalVolume* taper_log;

G4LogicalVolume* lg_log;

G4LogicalVolume* steel_log;

G4LogicalVolume* vetoPaddle_log;
G4LogicalVolume* vetol2Paddle_log;
G4LogicalVolume* vetoSurface_log;

// Physical Volumes
// ================
// to place optical surface between tapers and LGs for optical detection we need pointers to the physical volumes
// declare containers for the pointers here, they'll be filled by the placement function.
std::vector<G4VPhysicalVolume*> paddles_phys;
std::vector<G4VPhysicalVolume*> tapers_phys;
std::vector<G4VPhysicalVolume*> lgs_phys;
std::vector<G4VPhysicalVolume*> vetopaddles_phys;
std::vector<G4VPhysicalVolume*> vetosurfaces_phys;

// Define rotation matrices 
//=========================
// rotated and unrotated scintillator paddles. Could do away with one by changing dims but hey.
G4RotationMatrix* noRot = new G4RotationMatrix();											// null rotation pointer
G4RotationMatrix* rotatedmatx = new G4RotationMatrix(0,0,90*deg);			// horizontal config for scint paddles

// Note trapezium narrows along z axis, so need to define a rotation of 90deg about the y(??)-axis to bring taper into the x-y plane. 
G4RotationMatrix* upmtx = new G4RotationMatrix(180*deg,90*deg,0*deg);	
G4RotationMatrix* downmtx = new G4RotationMatrix(0*deg,90*deg,0*deg);
G4RotationMatrix* rightmtx = new G4RotationMatrix(90*deg,90*deg,0*deg);	
G4RotationMatrix* leftmtx = new G4RotationMatrix(-90*deg,90*deg,0*deg);

// Define Visualisation Attributes 
//================================
G4VisAttributes* scinthatts = new G4VisAttributes(G4Colour(0.5,0.0,1.0));     // pink
G4VisAttributes* scintvatts = new G4VisAttributes(G4Colour(1.0,0.0,0.5));     // purple
G4VisAttributes* steelatts = new G4VisAttributes(G4Colour(0.0,1.0,1.0));      // blue?
G4VisAttributes* scinttapatts = new G4VisAttributes(G4Colour(0.6, 1.0, 0.8)); // light-green
G4VisAttributes* scintlgatts = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));	// grey

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Done declaring global variables. Now declare functions for how each instance is placed
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Define Positioning of Steel Plates
//===================================
void ComputeSteelTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) {
		Xposition=0, Yposition=0, Zposition=0;
			
		Zposition=(copyNo)*(steelfullzlen + scintfullzlen + alufullzlen + layergap);	// layer width offset is always constant
		//Zposition=Zposition + (scintfullzlen + scintalugap + alufullzlen + alusteelgap); 							// offset of first layer 
		// no z offset of steel: front of first layer is front face of MRD.
		Zposition=Zposition + (steelfullzlen/2);													// offset by half depth so we are placing front face not centre
		Zposition=Zposition - (mrdZlen/2);																// offset by half total length to shift to front.

		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(0);
		physVol->GetLogicalVolume()->SetVisAttributes(steelatts);
	}
	
// Define Positioning of Scintillator Paddles
//===========================================
void ComputePaddleTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) {

		Xposition=0, Yposition=0, Zposition=0;
		G4RotationMatrix* rotmtx;
		G4int panelpairnum = floor(copyNo/(numpaddlesperpanelv+numpaddlesperpanelh));	// which pair of panels
		G4int panelnumrem = copyNo - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);	// copy num within a pair of panels
		G4int panelnum;
		G4int paddlenum;
		G4bool isvpaddle=false, ishpaddle=false;
		if(panelnumrem>numpaddlesperpanelv){																	// first layer is a VERTICAL layer (leading horizontal layer removed)
			panelnum = (panelpairnum*2) +1;
			ishpaddle = true;
			paddlenum = copyNo%numpaddlesperpanelh;
		} else {
			panelnum = (panelpairnum*2);
			isvpaddle = true;
			paddlenum = copyNo%numpaddlesperpanelv;
		}
		G4int pairnum = floor(paddlenum/2);																		// Paddles 0&1 are a pair; 
																																					// then X offset is the same for every pair
		Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen + layergap);
																																					// layer width offset is always constant (except first)
		//if(panelnum==0){Zposition = Zposition + alufullzlen + scintalugap;}								// first layer intrudes into 'layer offset' 
		Zposition = Zposition + steelfullzlen + steelscintgap;								// scint follows closely behind first steel
		Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we place front face not centre
		Zposition = Zposition - (mrdZlen/2);																	// offset by half total length to shift to front
		
		if(paddlenum==0){G4cout<<"scint layer placed at z = " << (Zposition + (mrdZlen/2) - (scintfullzlen/2) + tankouterRadius + 2*cm)<< " cm to "<< (Zposition + (mrdZlen/2) + (scintfullzlen/2) + tankouterRadius + 2*cm) << " cm." << G4endl;}
		
		// Y position is offset by half length, so that one end is at the origin. This needs to be the correct half-length
		// for the appropriate paddle (H or V type). Offsets are defined in the MOTHER ref frame. 
		// Then rotate paddles by 90deg in Y axis in every other panel. 
		if (ishpaddle){
			// horizontal panel
			if (paddlenum%2==0){
				Xposition=((scinthfullylen+scintbordergap)/2); 										// offset by +half length so one end is at x=0
			} else {
				Xposition=-((scinthfullylen+scintbordergap)/2);										// offset by -half length so one end is at x=0
			}
			Yposition = pairnum*(scintfullxlen+scintbordergap); 								// individual offset by pair number
			Yposition = Yposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelh)+(scintfullxlen/2);
			// shift whole set by 1/2 total X extent to shift center back to X=0: HalfLength cancels doubed num of paddles
			rotmtx=rotatedmatx;
		} else {
			// vertical panel
			if (paddlenum%2==0){
				Yposition=((scintvfullylen+scintbordergap)/2); 
			} else {
				Yposition=-((scintvfullylen+scintbordergap)/2);
			}
			// for vertical panels need to shift Y for each pair
			Xposition = pairnum*(scintfullxlen+scintbordergap); 	// individual offset by pair number
			// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
			Xposition = Xposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelv)+(scintfullxlen/2); 
			
			rotmtx=0;							// don't rotate vertical panels
		}
		
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(rotmtx);
//	physVol->GetLogicalVolume()->SetVisAttributes(scintvatts);	//can set visualisation attributes like this

}

// Define Positioning of Trapezoidal Taper ends of MRD paddles
// ===========================================================
// Combined with trapezoidal light-guide tapers as the code is 90% the same
void ComputeTaperTransformation (const G4int copyNo, G4VPhysicalVolume* physVol, G4bool lgs) {

		Xposition=0, Yposition=0, Zposition=0;
		G4RotationMatrix* rotmtx;
		G4int panelpairnum = floor(copyNo/(numpaddlesperpanelv+numpaddlesperpanelh));	// which pair of panels
		G4int panelnumrem = copyNo - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);	// copy num within a pair of panels
		G4int panelnum;
		G4int paddlenum;
		G4bool isvpaddle=false, ishpaddle=false;
		if(panelnumrem>numpaddlesperpanelv){																	// first layer is a VERTICAL layer (leading horizontal layer removed)
			panelnum = (panelpairnum*2) +1;
			ishpaddle = true;
			paddlenum = copyNo%numpaddlesperpanelh;
		} else {
			panelnum = (panelpairnum*2);
			isvpaddle = true;
			paddlenum = copyNo%numpaddlesperpanelv;
		}
		G4int pairnum = floor(paddlenum/2);											// LGs 0,1 are a vertical pair; X offset is the same for every pair
		
		// exact same z position calculations as paddles
		Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen + layergap);
																																					// layer width offset is always constant (except first)
		//if(panelnum==0){Zposition = Zposition + alufullzlen + scintalugap;}		// first layer intrudes into 'layer offset' 
		Zposition = Zposition + steelfullzlen + steelscintgap;								// scint follows closely behind first steel
		Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we place front face not centre
		Zposition=Zposition - (mrdZlen/2);																		// offset by half total length to shift to front
	
		// Y offset is the full length of the paddle plus half length of LG. Paddle length needs to be the correct type
		// (H or V type). Same rotation as paddles. 
		if (ishpaddle){
			// horizontal panel
			if(lgs){
			  Xposition=scinthfullylen+scinttapfullheight+(scintlgfullheight/2)+(scintbordergap/2);
			} else {
				Xposition=scinthfullylen+(scinttapfullheight/2)+(scintbordergap/2);
			}

			if (paddlenum%2==0){
				rotmtx=rightmtx;
			} else {
				Xposition=-Xposition;
				rotmtx=leftmtx;
			}
			// Y offset exactly the same as paddles
			Yposition = pairnum*(scintfullxlen+scintbordergap);
			Yposition = Yposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelh)+(scintfullxlen/2); 
		} else {
		// vertical panel
			if(lgs){
				Yposition=scintvfullylen+scinttapfullheight+(scintlgfullheight/2)+(scintbordergap/2);
			} else {
				Yposition=scintvfullylen+(scinttapfullheight/2)+(scintbordergap/2);
			}
			if (paddlenum%2==0){
				rotmtx=upmtx; 
			} else {
				Yposition=-Yposition;
				rotmtx=downmtx;
			}
			Xposition = pairnum*(scintfullxlen+scintbordergap);
			Xposition = Xposition - 0.5*(((scintfullxlen+scintbordergap)/2)*numpaddlesperpanelv)+(scintfullxlen/2); 
		}

		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetRotation(rotmtx);
		physVol->SetTranslation(origin);
		if(lgs){physVol->GetLogicalVolume()->SetVisAttributes(scintlgatts);}
		else {physVol->GetLogicalVolume()->SetVisAttributes(scinttapatts);}
	}


// Define Positioning of Veto Paddles
//===========================================
void ComputeVetoPaddleTransformation (const G4int copyNo, G4VPhysicalVolume* physVol, G4bool surf) {

		Xposition=0, Yposition=0, Zposition=0;
		G4RotationMatrix* rotmtx=0;																								// No rotations.
		G4int panelnum = floor(copyNo/vetopaddlesperpanel);												// numbering from 0
		G4int paddlenum = copyNo%vetopaddlesperpanel; 														// numering from 0 within a panel
		Zposition = panelnum*(vetopaddlefullzlen+vetopaddlegap); 									// layer width offset
		Zposition = Zposition + (vetopaddlefullzlen/2);														// offset by half depth to place front face not centre
		Zposition = Zposition - (vetoZlen/2);																			// offset by half total length to shift to front
			
		// AFAIK paddle dims are the same in both layers.
		Xposition = 0;																														// all paddles are centered on the beam
		Yposition = paddlenum*(vetopaddlefullylen+vetopaddlegap); 								// individual offset by paddle number
		// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
		Yposition = Yposition - 0.5*((vetopaddlefullylen+vetopaddlegap)*vetopaddlesperpanel)+(vetopaddlefullylen/2);
		if (panelnum==0){
 			Yposition = Yposition + vetolayer2offset;																// slight shift of second layer WRT first
 		}
		if(surf){
			Yposition = Yposition + (vetopaddlefullylen/2) + nothickness/2;							// shift empty volume for surface purposes
		}
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(rotmtx);
		// set vis attributes when placing, since they share a logical volume
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Done defining functions. Now do actual generation and placement of physical volumes
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Code to do generation and placement of scintillator paddles
// ===========================================================
void PlacePaddles(G4LogicalVolume* totMRD_log){

	hpaddle_log = new G4LogicalVolume(sciMRDhpaddle_box, G4Material::GetMaterial("Scinti"), "vpaddle_log");
	vpaddle_log = new G4LogicalVolume(sciMRDvpaddle_box, G4Material::GetMaterial("Scinti"), "hpaddle_log");
	hpaddle_log->SetVisAttributes(scinthatts);
	vpaddle_log->SetVisAttributes(scintvatts);
	G4LogicalVolume* paddle_log;
	G4VPhysicalVolume* paddle_phys;
	for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
			
			G4int panelpairnum = floor(i/(numpaddlesperpanelv+numpaddlesperpanelh));					// which pair of panels
			G4int panelnumrem = i - panelpairnum*(numpaddlesperpanelv+numpaddlesperpanelh);		// copy num within a pair of panels
			G4bool isvpaddle=false, ishpaddle=false;
			if(panelnumrem>numpaddlesperpanelv){																	// first layer is a VERTICAL layer (leading horizontal layer removed)
				paddle_log = hpaddle_log;
			} else {
				paddle_log = vpaddle_log;
			}
			
			paddle_phys = new G4PVPlacement(noRot,G4ThreeVector(), paddle_log, "paddle_phys", totMRD_log, false, i);
			ComputePaddleTransformation(i, paddle_phys);
			paddles_phys.push_back(paddle_phys);
	}
}

// Code to do generation and placement of scintillator taper ends
// ==============================================================
void PlaceTapers(G4LogicalVolume* totMRD_log){

	taper_log = new G4LogicalVolume(mrdScintTap_box, G4Material::GetMaterial("Scinti"), "taper_log");			
	G4VPhysicalVolume* taper_phys;
	for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
			taper_phys = new G4PVPlacement(noRot,G4ThreeVector(), taper_log, "taper_phys", totMRD_log, false, i);
			ComputeTaperTransformation(i, taper_phys, false);
			tapers_phys.push_back(taper_phys);
	}
}

// Code to do generation and placement of glass light-guides
// =========================================================
void PlaceLGs(G4LogicalVolume* totMRD_log){

	lg_log = new G4LogicalVolume(mrdLG_box, G4Material::GetMaterial("Glass"), "lg_log");
	G4VPhysicalVolume* lg_phys;
	for(G4int i=0;i<((numpaddlesperpanelv+numpaddlesperpanelh)*(numpanels/2));i++){
			lg_phys = new G4PVPlacement(noRot,G4ThreeVector(), lg_log, "MRDLG_PV", totMRD_log, false, i);
			ComputeTaperTransformation(i, lg_phys, true);
			lgs_phys.push_back(lg_phys);
	}
}

// Code to do generation and placement of steel plates
// ===================================================
void PlaceSteels(G4LogicalVolume* totMRD_log){

	steel_log = new G4LogicalVolume(steelMRDplate_box, G4Material::GetMaterial("Steel"), "steel_log");
	G4VPhysicalVolume* steel_phys;
	for(G4int i=0;i<(numplates);i++){		// first layer of steel has been removed
			steel_phys = new G4PVPlacement(noRot,G4ThreeVector(), steel_log, "steel_phys", totMRD_log, false, i);
			ComputeSteelTransformation(i, steel_phys);
	}
}

// Define Alu Support Structure
//=============================

void makeAlu(G4AssemblyVolume* aluMRDassembly){
	// this code must be in a function as it uses for loops.

	//N.B. subtraction solids can only use CGS solids or the output of boolean operations, not assemblies or other stuff
	G4RotationMatrix  Ra (0,0,0);
	G4ThreeVector  Ta (0,0,0);
	Ta.set(windowwidth+alufullxthickness,windowheight+alufullythickness,0);
	G4Transform3D transform = G4Transform3D(Ra,Ta);
	G4SubtractionSolid* aluMRDstruc_sol = new G4SubtractionSolid("aluStruct",aluMRDstruc_box,aluMRDwindow_box,transform);
	G4SubtractionSolid* aluMRDstruc_sol2;
	
	for (G4int row=-1;row<2;row++){
		for (G4int col=-1;col<2;col++){
			G4double xoffset=col*(alufullxthickness+windowwidth);
			G4double yoffset=row*(alufullythickness+windowheight);
			Ta.set(xoffset,yoffset,0);                        
			transform = G4Transform3D(Ra,Ta);
			// using method with G4Transform3D means we can change the transform after as it is passed byval.
			// If we used a passive transform we would need to maintain the G4Transform3D. 
			aluMRDstruc_sol2 = new G4SubtractionSolid("aluStruct",aluMRDstruc_sol,aluMRDwindow_box,transform);
			//delete aluMRDstruc_sol; -- can't do this! seems derived boolean solid requires it's originals are kept alive!
			aluMRDstruc_sol = aluMRDstruc_sol2;
		}
	}

	G4LogicalVolume* aluMRDstruc_log 
	= new G4LogicalVolume(aluMRDstruc_sol, G4Material::GetMaterial("Aluminum"), "aluStruct",0,0,0) ;

	Ra.set(0,0,0);
	Ta.set(0,0,0);
	for (G4int structnum=0;structnum<(numalustructs-1);structnum++){		// first alu layer removed
		G4double zpos = structnum*(steelfullzlen + alufullzlen + scintfullzlen + layergap); // layer width offset is always constant
		//if(structnum>0){zpos = zpos + scintfullzlen + scintalugap;}	// all layers > 1 have additional scint offset
		zpos+= steelfullzlen + steelscintgap + scintfullzlen + scintalugap;
		zpos+= (alufullzlen/2); 																		// offset by half depth to place front face not centre
		zpos+= -(mrdZlen/2);																				// offset by half total length to shift to front
		Ta.set(0,0,zpos); 
		aluMRDassembly->AddPlacedVolume(aluMRDstruc_log,Ta,&Ra);
	}
	//aluMRDassembly->MakeImprint(expHall->GetLogicalVolume(), Tm,&Rm);  //placement done in DetectorConstruction
}

// Code to do generation and placement of veto paddles
// ===================================================
void PlaceVetoPaddles(G4LogicalVolume* totVeto_log){

	// using same scintillator as MRD paddles
	vetoPaddle_log = new G4LogicalVolume(vetoPaddle_box, G4Material::GetMaterial("Scinti"), "vetoPaddle_log");
	vetol2Paddle_log = new G4LogicalVolume(vetoPaddle_box, G4Material::GetMaterial("Scinti"), "vetol2Paddle_log");
	vetoPaddle_log->SetVisAttributes(scinthatts);
	vetol2Paddle_log->SetVisAttributes(scintvatts);
	G4LogicalVolume* paddle_log;
	// gotta make stupid little physical volumes so we can put a surface in. It shouldn't take up any space. 
	vetoSurface_log = new G4LogicalVolume(vetoSurface_box, G4Material::GetMaterial("Air"),"vetoSurface_log");
	G4VisAttributes* surfacevisatts= new G4VisAttributes(G4Colour(0.6, 1.0, 0.8));
	vetoSurface_log->SetVisAttributes(surfacevisatts);
	G4VPhysicalVolume* vetoPaddle_phys;
	G4VPhysicalVolume* vetoSurface_phys;
	
	for(G4int i=0;i<numvetopaddles;i++){
			G4int panelnum = floor(i/vetopaddlesperpanel);
			if(panelnum==0){
				paddle_log=vetoPaddle_log;
			} else {
				paddle_log=vetol2Paddle_log;
			}
			vetoPaddle_phys = new G4PVPlacement(noRot,G4ThreeVector(), paddle_log, "vetoPaddle_phys", totVeto_log, false, i);
			ComputeVetoPaddleTransformation(i, vetoPaddle_phys,false);
			vetopaddles_phys.push_back(vetoPaddle_phys);
			
			vetoSurface_phys = new G4PVPlacement(noRot,G4ThreeVector(), vetoSurface_log, "vetoSurface_phys", totVeto_log, false, i);
			ComputeVetoPaddleTransformation(i, vetoSurface_phys,true);
			vetosurfaces_phys.push_back(vetoSurface_phys);
	}
}

/* MRD z edges, from firing a geantino (set geantinogun in primary generator action)
volume column is the volume it is LEAVING
Step#      X         Y         Z        KineE    dEStep   StepLeng  TrakLeng    Volume     Process
    0      3 cm      3 cm      0 fm      1 GeV     0 eV      0 fm      0 fm    waterTank    initStep
    1      3 cm      3 cm   1.52 m       1 GeV     0 eV   1.52 m    1.52 m     waterTank  Transportation
    2      3 cm      3 cm   1.54 m       1 GeV     0 eV   2.03 cm   1.54 m          Hall  Transportation
    3      3 cm      3 cm   1.56 m       1 GeV     0 eV    2.1 cm   1.56 m    totMRDphys  Transportation
    4      3 cm      3 cm   1.57 m       1 GeV     0 eV      6 mm   1.57 m    paddle_phys  Transportation
    5      3 cm      3 cm   1.59 m       1 GeV     0 eV      2 cm   1.59 m    totMRDphys  Transportation
    6      3 cm      3 cm   1.64 m       1 GeV     0 eV      5 cm   1.64 m    steel_phys  Transportation
    7      3 cm      3 cm   1.64 m       1 GeV     0 eV      2 mm   1.64 m    totMRDphys  Transportation
    8      3 cm      3 cm   1.65 m       1 GeV     0 eV      6 mm   1.65 m    paddle_phys  Transportation
    9      3 cm      3 cm   1.69 m       1 GeV     0 eV    4.1 cm   1.69 m    totMRDphys  Transportation
   10      3 cm      3 cm   1.74 m       1 GeV     0 eV      5 cm   1.74 m    steel_phys  Transportation
   11      3 cm      3 cm   1.74 m       1 GeV     0 eV      2 mm   1.74 m    totMRDphys  Transportation
   12      3 cm      3 cm   1.75 m       1 GeV     0 eV      6 mm   1.75 m    paddle_phys  Transportation
   13      3 cm      3 cm   1.79 m       1 GeV     0 eV    4.1 cm   1.79 m    totMRDphys  Transportation
   14      3 cm      3 cm   1.84 m       1 GeV     0 eV      5 cm   1.84 m    steel_phys  Transportation
   15      3 cm      3 cm   1.84 m       1 GeV     0 eV      2 mm   1.84 m    totMRDphys  Transportation
   16      3 cm      3 cm   1.85 m       1 GeV     0 eV      6 mm   1.85 m    paddle_phys  Transportation
   17      3 cm      3 cm   1.89 m       1 GeV     0 eV    4.1 cm   1.89 m    totMRDphys  Transportation
   18      3 cm      3 cm   1.94 m       1 GeV     0 eV      5 cm   1.94 m    steel_phys  Transportation
   19      3 cm      3 cm   1.94 m       1 GeV     0 eV      2 mm   1.94 m    totMRDphys  Transportation
   20      3 cm      3 cm   1.95 m       1 GeV     0 eV      6 mm   1.95 m    paddle_phys  Transportation
   21      3 cm      3 cm   1.99 m       1 GeV     0 eV    4.1 cm   1.99 m    totMRDphys  Transportation
   22      3 cm      3 cm   2.04 m       1 GeV     0 eV      5 cm   2.04 m    steel_phys  Transportation
   23      3 cm      3 cm   2.04 m       1 GeV     0 eV      2 mm   2.04 m    totMRDphys  Transportation
   24      3 cm      3 cm   2.04 m       1 GeV     0 eV      6 mm   2.04 m    paddle_phys  Transportation
   25      3 cm      3 cm   2.09 m       1 GeV     0 eV    4.1 cm   2.09 m    totMRDphys  Transportation
   26      3 cm      3 cm   2.14 m       1 GeV     0 eV      5 cm   2.14 m    steel_phys  Transportation
   27      3 cm      3 cm   2.14 m       1 GeV     0 eV      2 mm   2.14 m    totMRDphys  Transportation
   28      3 cm      3 cm   2.14 m       1 GeV     0 eV      6 mm   2.14 m    paddle_phys  Transportation
   29      3 cm      3 cm   2.19 m       1 GeV     0 eV    4.1 cm   2.19 m    totMRDphys  Transportation
   30      3 cm      3 cm   2.23 m       1 GeV     0 eV      5 cm   2.23 m    steel_phys  Transportation
   31      3 cm      3 cm   2.24 m       1 GeV     0 eV      2 mm   2.24 m    totMRDphys  Transportation
   32      3 cm      3 cm   2.24 m       1 GeV     0 eV      6 mm   2.24 m    paddle_phys  Transportation
   33      3 cm      3 cm   2.28 m       1 GeV     0 eV    4.1 cm   2.28 m    totMRDphys  Transportation
   34      3 cm      3 cm   2.33 m       1 GeV     0 eV      5 cm   2.33 m    steel_phys  Transportation
   35      3 cm      3 cm   2.34 m       1 GeV     0 eV      2 mm   2.34 m    totMRDphys  Transportation
   36      3 cm      3 cm   2.34 m       1 GeV     0 eV      6 mm   2.34 m    paddle_phys  Transportation
   37      3 cm      3 cm   2.38 m       1 GeV     0 eV    4.1 cm   2.38 m    totMRDphys  Transportation
   38      3 cm      3 cm   2.43 m       1 GeV     0 eV      5 cm   2.43 m    steel_phys  Transportation
   39      3 cm      3 cm   2.44 m       1 GeV     0 eV      2 mm   2.44 m    totMRDphys  Transportation
   40      3 cm      3 cm   2.44 m       1 GeV     0 eV      6 mm   2.44 m    paddle_phys  Transportation
   41      3 cm      3 cm   2.48 m       1 GeV     0 eV    4.1 cm   2.48 m    totMRDphys  Transportation
   42      3 cm      3 cm   2.53 m       1 GeV     0 eV      5 cm   2.53 m    steel_phys  Transportation
   43      3 cm      3 cm   2.53 m       1 GeV     0 eV      2 mm   2.53 m    totMRDphys  Transportation
   44      3 cm      3 cm   2.54 m       1 GeV     0 eV      6 mm   2.54 m    paddle_phys  Transportation
   45      3 cm      3 cm   2.58 m       1 GeV     0 eV    4.1 cm   2.58 m    totMRDphys  Transportation
   46      3 cm      3 cm   2.63 m       1 GeV     0 eV      5 cm   2.63 m    steel_phys  Transportation
   47      3 cm      3 cm   2.63 m       1 GeV     0 eV      2 mm   2.63 m    totMRDphys  Transportation
   48      3 cm      3 cm   2.64 m       1 GeV     0 eV      6 mm   2.64 m    paddle_phys  Transportation
   49      3 cm      3 cm   2.68 m       1 GeV     0 eV    4.1 cm   2.68 m    totMRDphys  Transportation
   50      3 cm      3 cm   2.73 m       1 GeV     0 eV      5 cm   2.73 m    steel_phys  Transportation
   51      3 cm      3 cm   2.76 m       1 GeV     0 eV    2.9 cm   2.76 m    totMRDphys  Transportation
   52      3 cm      3 cm     25 m       1 GeV     0 eV   22.2 m      25 m    OutOfWorld  Transportation
*/
