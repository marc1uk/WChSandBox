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

G4double Xposition=0, Yposition=0, Zposition=0;		// used for positioning parameterisations.
G4int numpaddlesperpanel=30;								// paddles per scintillator panel
G4int numpanels=12;													// scintillator panels
G4int numrpcs=10;														// rpc panels
G4int numplates=12;													// steel plates
G4int numalustructs=13;											// number of supporting structs. We may be dropping one as we have fewer scintillators?

G4double steelfullxlen = 274*cm;
G4double steelfullylen = 305*cm;
G4double steelfullzlen = 5*cm;

G4double scintfullxlen = 20*cm;
G4double scintfullzlen= 0.6*cm;
G4double scintvfullylen = 155*cm;
G4double scinthfullylen= 138*cm;

G4double scinttapfullwidth = 17.1*cm; 			// width of the tapering part of scintillator paddles at the narrow end
G4double scinttapfullheight = 7.8*cm; 			// z length of tapering part of scint paddles.

G4double scintlgfullwidth = 10.16*cm; 			// tapered light guides
G4double scintlgfullheight = 16.67*cm; 			// mrdmodule file says 33.34 full length, but looks too long...??

G4double alufullxlen = steelfullxlen+5*cm;	// outer thicknesses - total frame dims are about those of the steel plate
G4double alufullylen = steelfullylen+5*cm;	// basically making this up
G4double alufullzlen = 1.9*cm;							// from eye and a skim of the mrdmodule.txt file, i'm guessing depth is ~0.75 inches
G4double alufullxthickness = 2.54*cm;				// as above, guessing frame to be 1 inch box cross-section
G4double alufullythickness = 2.54*cm;
G4double windowwidth = (steelfullxlen-(4*alufullxthickness))/3;	// (full length - 4 beams) / 3 windows
G4double windowheight= (steelfullylen-(4*alufullythickness))/3;
	
G4double mrdZlen = numplates*steelfullzlen + (numpanels+1)*scintfullzlen + numalustructs*alufullzlen; 
// add another panel to full length because last alu struct is placed back assuming it's there. Maybe need to change... 

G4double widths[] = {2*(scinthfullylen+scinttapfullheight+scintlgfullheight),((numpaddlesperpanel/2)*scintfullxlen)};	
// 2* and y dim because we stack 2 rotated h scint paddles end-to-end. 
//
G4double heights[] = {2*(scintvfullylen+scinttapfullheight+scintlgfullheight),((numpaddlesperpanel/2)*scintfullxlen)};
//
G4double maxwidth = *std::max_element(widths,widths+(sizeof(widths)/sizeof(widths[0])));
G4double maxheight = *std::max_element(heights,heights+(sizeof(heights)/sizeof(heights[0])));
        
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

G4Box* totMRD_box = new G4Box("totMRD",maxwidth/2,maxheight/2,mrdZlen/2);
	
// Define logical volumes 
//=======================
// G4LogicalVolume* variableName = new G4LogicalVolume(solidVariableName, Material, "logicalVolName");
// Can't do this outside a function - the materials aren't defined yet.

G4LogicalVolume* hpaddle_log;
			
G4LogicalVolume* vpaddle_log;

G4LogicalVolume* taper_log;

G4LogicalVolume* lg_log;

G4LogicalVolume* steel_log;

// Physical Volumes
// ================
// to place optical surface between tapers and LGs for optical detection we need pointers to the physical volumes
// declare containers for the pointers here, they'll be filled by the placement function.
std::vector<G4VPhysicalVolume*> tapers_phys;
std::vector<G4VPhysicalVolume*> lgs_phys;

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
			
		Zposition=(copyNo)*(steelfullzlen + scintfullzlen + alufullzlen);	// layer width offset is always constant
		Zposition=Zposition + (scintfullzlen + alufullzlen); 							// offset of first layer 
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
		G4int panelnum = floor(copyNo/numpaddlesperpanel);			// numbering from 0
		G4int paddlenum = copyNo%numpaddlesperpanel; 						// numering from 0 within a panel
		G4int pairnum = floor(paddlenum/2);					// Consider paddles 0&1 as a vertical pair; then X offset is the same for every pair		
		
		Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen);		// layer width offset is always constant (except first)
		if(panelnum==0){Zposition = Zposition + alufullzlen;}									// first layer intrudes into 'layer offset' 
		Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we are placing front face not centre
		Zposition=Zposition - (mrdZlen/2);																		// offset by half total length to shift to front
			
		// Y position is offset by half length, so that one end is at the origin. This needs to be the correct half-length
		// for the appropriate paddle (H or V type). Offsets are defined in the MOTHER ref frame. 
		// Then rotate paddles by 90deg in Y axis in every other panel. 
		if (panelnum%2==0){
			// horizontal panel
			if (paddlenum%2==0){
				Xposition=(scinthfullylen/2); 		// offset by +half length so one end is at x=0
			} else {
				Xposition=-(scinthfullylen/2); 		// offset by -half length so one end is at x=0
			}
			Yposition = pairnum*2*(scintfullxlen/2); 	// individual offset by pair number
			// shift whole set by 1/2 total X extent to shift center back to X=0: HalfLength cancels doubed num of paddles
			Yposition = Yposition - 0.5*((scintfullxlen/2)*numpaddlesperpanel)+(scintfullxlen/2); 
			
			rotmtx=rotatedmatx;
		} else {
			// vertical panel
			if (paddlenum%2==0){
				Yposition=(scintvfullylen/2); 
			} else {
				Yposition=-(scintvfullylen/2);
			}
			// for vertical panels need to shift Y for each pair
			Xposition = pairnum*2*(scintfullxlen/2); 	// individual offset by pair number
			// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
			Xposition = Xposition - 0.5*((scintfullxlen/2)*numpaddlesperpanel)+(scintfullxlen/2); 
			
			rotmtx=0;							// don't rotate vertical panels
		}
		
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(rotmtx);
		if (panelnum%2==0){
			physVol->GetLogicalVolume()->SetVisAttributes(scintvatts);	//can set visualisation attributes like this
		} else {
			physVol->GetLogicalVolume()->SetVisAttributes(scinthatts);	
		}
}

// Define Positioning of Trapezoidal Taper ends of MRD paddles
// ===========================================================
// Combined with trapezoidal light-guide tapers as the code is 90% the same
void ComputeTaperTransformation (const G4int copyNo, G4VPhysicalVolume* physVol, G4bool lgs) {

		Xposition=0, Yposition=0, Zposition=0;
		G4RotationMatrix* rotmtx;
		G4int panelnum = floor(copyNo/numpaddlesperpanel);			// same calculations as for paddles.
		G4int paddlenum = copyNo%numpaddlesperpanel; 						// numer from 0 within a panel
		G4int pairnum = floor(paddlenum/2);											// LGs 0,1 are a vertical pair; X offset is the same for every pair		
		
		// exact same z position calculations as paddles
		Zposition = panelnum*(steelfullzlen + alufullzlen + scintfullzlen);		// layer width offset is always constant (except first)
		if(panelnum==0){Zposition = Zposition + alufullzlen;}									// first layer intrudes into 'layer offset' 
		Zposition = Zposition + (scintfullzlen/2);														// offset by half depth so we are placing front face not centre
		Zposition=Zposition - (mrdZlen/2);									// offset by half total length to shift to front
	
		// Y offset is the full length of the paddle plus half length of LG. Paddle length needs to be the correct type
		// (H or V type). Same rotation as paddles. 
		if (panelnum%2==0){
			if(lgs){
			  Xposition=scinthfullylen+scinttapfullheight+(scintlgfullheight/2);
			} else {
				Xposition=scinthfullylen+(scinttapfullheight/2);
			}
			// horizontal panel
			if (paddlenum%2==0){
				rotmtx=rightmtx;
			} else {
				Xposition=-Xposition;
				rotmtx=leftmtx;
			}
			// Y offset exactly the same as paddles
			Yposition = pairnum*2*(scintfullxlen/2);
			Yposition = Yposition - 0.5*((scintfullxlen/2)*numpaddlesperpanel)+(scintfullxlen/2); 
		} else {
			if(lgs){
				Yposition=scintvfullylen+scinttapfullheight+(scintlgfullheight/2);
			} else {
				Yposition=scintvfullylen+(scinttapfullheight/2);
			}
			// vertical panel
			if (paddlenum%2==0){
				rotmtx=upmtx; 
			} else {
				Yposition=-Yposition;
				rotmtx=downmtx;
			}
			Xposition = pairnum*2*(scintfullxlen/2);
			Xposition = Xposition - 0.5*((scintfullxlen/2)*numpaddlesperpanel)+(scintfullxlen/2); 
		}

		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetRotation(rotmtx);
		physVol->SetTranslation(origin);
		if(lgs){physVol->GetLogicalVolume()->SetVisAttributes(scintlgatts);}
		else {physVol->GetLogicalVolume()->SetVisAttributes(scinttapatts);}
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
	G4LogicalVolume* paddle_log;
	G4VPhysicalVolume* paddle_phys;
	for(G4int i=0;i<(numpaddlesperpanel*numpanels);i++){
			G4int panelnum = floor(i/numpaddlesperpanel);			// numbering from 0
			if (panelnum%2==0){
				paddle_log = hpaddle_log;
			} else {
				paddle_log = vpaddle_log;
			}
			paddle_phys = new G4PVPlacement(noRot,G4ThreeVector(), paddle_log, "paddle_phys", totMRD_log, false, i);
			ComputePaddleTransformation(i, paddle_phys);
	}
}

// Code to do generation and placement of scintillator taper ends
// ==============================================================
void PlaceTapers(G4LogicalVolume* totMRD_log){

	taper_log = new G4LogicalVolume(mrdScintTap_box, G4Material::GetMaterial("Scinti"), "taper_log");			
	G4VPhysicalVolume* taper_phys;
	for(G4int i=0;i<(numpaddlesperpanel*numpanels);i++){
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
	for(G4int i=0;i<(numpaddlesperpanel*numpanels);i++){
			lg_phys = new G4PVPlacement(noRot,G4ThreeVector(), lg_log, "lg_phys", totMRD_log, false, i);
			ComputeTaperTransformation(i, lg_phys, true);
			lgs_phys.push_back(lg_phys);
	}
}

// Code to do generation and placement of steel plates
// ===================================================
void PlaceSteels(G4LogicalVolume* totMRD_log){

	steel_log = new G4LogicalVolume(steelMRDplate_box, G4Material::GetMaterial("Steel"), "steel_log");
	G4VPhysicalVolume* steel_phys;
	for(G4int i=0;i<(numplates);i++){
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
// using method with G4Transform3D means we can change the transform after as it is passed byval. If we used a passive transform we would need to maintain the G4Transform3D. 
			aluMRDstruc_sol2 = new G4SubtractionSolid("aluStruct",aluMRDstruc_sol,aluMRDwindow_box,transform);
			//delete aluMRDstruc_sol; -- can't do this! seems derived boolean solid requires it's originals are kept alive!
			aluMRDstruc_sol = aluMRDstruc_sol2;
		}
	}

	G4LogicalVolume* aluMRDstruc_log 
	= new G4LogicalVolume(aluMRDstruc_sol, G4Material::GetMaterial("Aluminum"), "aluStruct",0,0,0) ;

	Ra.set(0,0,0);
	Ta.set(0,0,0);
	for (G4int structnum=0;structnum<numalustructs;structnum++){
		G4double zpos = structnum*(steelfullzlen + alufullzlen + scintfullzlen);	// layer width offset is always constant (except first)
		if(structnum>0){zpos = zpos + scintfullzlen;}															// all layers > 1 have additional scint offset 
		zpos+= (alufullzlen/2); 																									// offset by half depth to place front face not centre
		zpos+= -(mrdZlen/2);																		// offset by half total length to shift to front
		Ta.set(0,0,zpos); 
		aluMRDassembly->AddPlacedVolume(aluMRDstruc_log,Ta,&Ra);
	}
	//aluMRDassembly->MakeImprint(expHall->GetLogicalVolume(), Tm,&Rm);  //placement <- in DetectorConstruction
}

/* scintillator panel z edges, in MRD ref frame (subtract 2000mm from global)
layer 0 at z=19 and 25
layer 1 at z=75 and 81
layer 2 at z=150 and 156
layer 3 at z=225 and 23
layer 4 at z=300 and 306
layer 5 at z=375 and 381
layer 6 at z=450 and 456
layer 7 at z=525 and 531
layer 8 at z=600 and 606
layer 9 at z=675 and 681
layer 10 at z=750 and 756
layer 11 at z=825 and 831 
*/
