// ====================================================================
//   MRDDetectorConstruction.cc
//
//   17/11/15 M. O'Flaherty
// ====================================================================
//===============================================================================================================
  //MRD DETECTOR DEFINITION
//===============================================================================================================
/* Notes on dimensions from papers
  From 'The Muon Range Detector at SciBooNE': 
  There are 362 channels in total; 182 horizontal, 180 vertical. 
  Each channel is 20 cm wide and either 155 cm (horizontal) or 138 cm (vertical) long. 
  Planes have active areas of 310×260 cm 2 (horizontal) and 300×276 cm 2 (vertical).
  
  From ANNIE proposal:
  The SciBooNE MRD consists of twelve iron plates, each 2 inches thick, sandwiched between thirteen layers of scintillator. 
  The scintillator planes were comprised of panels, each 20 cm wide and 0.6 cm thick, and were were arranged in alternating
  vertical and horizontal layers. The SciBooNE MRD was read out by 362 PMTs, each 2 inches in diameter. The iron plates 
  cover an area of 274 × 305 cm^2. ...
  The ANNIE MRD would keep the twelve steel plates intact and in their present configuration. 
  Instead of thirteen layers of scintillator, the new design calls for only three layers, with the remaining ten layers
  comprised of RPCs.
  Each RPC is layer is one square meter. The first (i.e., most upstream) layer of the detector is to be composed of scintillator
  and used as a trigger. This will be followed by five layers of RPCs. A second scintillator layer will serve as a coincidence
  trigger for noise reduction if necessary [determined during test phase I]. Following the second scintillator layer, there will be
  five more layers of RPCs, ending with a final layer of scintillator. The last layer will tag exiting events, for which only a lower
  bound on energy can be set.
  
  From ANNIE Whitepaper:
  This consists of 12 slabs of 2-inch-thick iron plates sandwiched between 13 alternating vertical and horizontal layers of 20-cm-wide,
  0.6-cm-thick plastic scintillation panels read out by 362 two-inch PMTs. The vertical layers are arranged in a 2x15 pattern of
  138-cm-long paddles, while the horizontal layers are in a 2x13 array of 155-cm-long paddles.
*/
//===============================================================================================================
/* File summary:
Define solids and logical volumes for steel plates, scintillator padels and rpc planes. 
We make use of the repetative nature of the construction, using parameterisation to define just one "physical"
volume to represent all 12 steel panels, all 10 rpcs, and all 3 scintillator padels. 
For steel panels and rpcs, we define a logical volume to represent one panel, then a "parameterisation class" 
that encapsulates functions for calculating where each element should be placed and its dimensions. The physical
class is then a 'G4VParameterised' class; a child of the G4VPhysicalVolume class, defined in terms of an 
element logical volume and the parameterisation class. 
For scintillator panels there is another layer of repetition: we have 3 identical panels, each comprised of 30 
identical paddles. We therefore define a logical volume to represent the panel, and use parameterisation to 
propagate this into a set as above. The logical volume representing the panel is then used as the mother volume 
for a parameterisation class representing a set of multiple paddles....
Unfortunately, however, the paddle dimensions are different for vertical and horizontal orientation, necessitating
selection of the solid based on the parent (panel) copyNo. At present, this is not supported (only the material
selection method may access the parent copyNo), so we need to use two parameterisations of panels. 
(ok, so this is a bit overkill since one panel will be a parameterised set of 1. Future versions of geant may
allow us to combine the two and only have one physical volume to represent all 3 panels) 
*/
 
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
#include "MRDDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"

G4double Xposition=0, Yposition=0, Zposition=0;		// used for positioning parameterisations.
G4int numpaddlesperpanel=30;	// paddles per scintillator panel
G4int numpanels=12;		// scintillator panels
G4int numrpcs=10;		// rpc panels
G4int numplates=12;		// steel plates
G4int numalustructs=13;		// number of supporting structs. We may be dropping one as we have fewer scintillators?

G4double steelfullxlen = 274*cm;
G4double steelfullylen = 305*cm;
G4double steelfullzlen = 5*cm;

G4double scintvfullxlen = 20*cm;
G4double scintvfullylen = 155*cm;
G4double scintvfullzlen= 0.6*cm;

G4double scinthfullxlen = 20*cm;
G4double scinthfullylen= 138*cm;
G4double scinthfullzlen= 0.6*cm;

G4double alufullxlen = steelfullxlen;	// outer thicknesses - total frame dims are about those of the steel plate
G4double alufullylen = steelfullylen;	// basically making this up
G4double alufullzlen = 1.9*cm;		// from eye and a skim of the mrdmodule.txt file, i'm guessing depth is ~0.75 inches
G4double alufullxthickness = 2.54*cm;	// as above, guessing frame to be 1 inch box cross-section
G4double alufullythickness = 2.54*cm;
G4double windowwidth = (steelfullxlen-(4*alufullxthickness))/3;	// (full length - 4 beams) / 3 windows
G4double windowheight= (steelfullylen-(4*alufullythickness))/3;
	
G4double mrdZlen = numplates*steelfullzlen + numpanels*scintvfullzlen + numalustructs*alufullzlen; 

//G4double maxwidth = 310*cm;
//G4double maxheight = 305*cm;
G4double widths[] = {steelfullxlen, 2*scinthfullylen, (numpaddlesperpanel/2)*scintvfullxlen, alufullxlen};	// 2* and y dim because we stack 2 rotated h scint paddles end-to-end. 
G4double heights[] = {steelfullylen, 2*scintvfullylen, (numpaddlesperpanel/2)*scinthfullxlen, alufullylen};
G4double maxwidth = *std::max_element(widths,widths+(sizeof(widths)/sizeof(widths[0])));
G4double maxheight = *std::max_element(heights,heights+(sizeof(heights)/sizeof(heights[0])));

//size method does not work if the array is declared with 'new' as sizeof(heights) returns the sizeof the pointer, not the number of bytes allocated for the array. For dynamically created either keep track, or use std::array/vector.

// N.B. Hall is 50*500*500m
        
// Define solids
//==============  
// G4Box* variableName = new G4Box("SolidName", x_length, y_length, z_length);

G4Box* steelMRDplate_box = new G4Box("steelPlate",steelfullxlen/2,steelfullylen/2,steelfullzlen/2);

G4Box* sciMRDhpaddle_box = new G4Box("scintHpaddle",scinthfullxlen/2,scinthfullylen/2,scinthfullzlen/2);

G4Box* sciMRDvpaddle_box = new G4Box("scintVpaddle",scintvfullxlen/2,scintvfullylen/2,scintvfullzlen/2);

//	G4Box* rpcMRDplate_box = new G4Box("rpcPlate",100*cm,100*cm,1*cm); 		// ASSUMED RPC DIMS  

//The alu support structure is roughly the external size of the steel plates...
G4Box* aluMRDstruc_box = new G4Box("outer_Box", alufullxlen/2, alufullylen/2, alufullzlen/2);

// ...with a 3x3 grid of ~even sized holes  
G4Box* aluMRDwindow_box = new G4Box("inner_Box", windowwidth/2, windowheight/2, alufullzlen/2);

G4Box* totMRD_box = new G4Box("totMRD",maxwidth/2,maxheight/2,mrdZlen/2);
	
// Define logical volumes (all done in parameterisation)
//=======================
// G4LogicalVolume* variableName = new G4LogicalVolume(solidVariableName, Material, "logicalVolName");
// Note:  We can only declare/define variables outside a function. This means we can't often define logical volumes outside of a function. Logical volumes require a material; materials are normally generated by the constructmaterials() function, but we can't call it outside a function. Nor can we define materials separately, as this requires defining them (OK), AND operating on them (not OK) - e.g. G4Material Air(blah); Air->AddElement(elN) << this can't be done. Unless we #include something that pulls our materials in, or maybe use 'extern'?

/*
G4LogicalVolume* steelMRDplate_log
= new G4LogicalVolume(steelMRDplate_box,G4Material::GetMaterial("Iron"),"steelPlate",0,0,0);		//material (Iron) not important because it's defined within parameterisation

G4LogicalVolume* sciMRDhpaddle_log
= new G4LogicalVolume(sciMRDhpaddle_box,G4Material::GetMaterial("Scinti"),"scintHpaddle",0,0,0);	//material (Scinti) not important because it's defined within parameterisation

G4LogicalVolume* sciMRDvpaddle_log
= new G4LogicalVolume(sciMRDvpaddle_box,G4Material::GetMaterial("Scinti"),"scintVpaddle",0,0,0);	//material (Scinti) not important because it's defined within parameterisation

G4LogicalVolume* totMRD_log
= new G4LogicalVolume(totMRD_box, G4Material::GetMaterial("Air"),"totMRD",0,0,0);
*/
	
// Define steel plates
//=======================	
	// 1. Define G4VPVParameterisation class - in header
	// This defines functions to retrieve the specifics of individual steel plates
/*	class steelPlateParameterisation : public G4VPVParameterisation{
		public:
		steelPlateParameterisation();							// constructor
		virtual ~steelPlateParameterisation();						// destructor
		//G4VPhysicalVolume here is the physical volume of the instance
		virtual void ComputeTransformation						// mandatory
		  (const G4int copyNo, G4VPhysicalVolume* physVol) const;
		virtual void ComputeDimensions 							// mandatory
		  (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
	};
*/	
/* Notes:

ComputeTransformation computes the transformation for this copy, and sets the physical output volume resulting from this transformation. 
Similarly the ComputeDimensions method is used to set the size of that copy. 
ComputeSolid can be defined to specify the type of solid. (optional)
ComputeMaterial can be used to specify the type of material. (optional - see G4VNestedParameterisation)

Currently daughter volumes are only allowed for parameterised volumes when all parameterised solids are the same (type and size). 
*/
	
	steelPlateParameterisation::steelPlateParameterisation(){};				// constructor
	steelPlateParameterisation::~steelPlateParameterisation(){};				// destructor
	
	// Compute translation of steel plate based on plate number
	void steelPlateParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const {
		Xposition=0, Yposition=0, Zposition=0;
		Zposition=10*(copyNo-1)*2*(steelMRDplate_box->GetZHalfLength() + sciMRDhpaddle_box->GetZHalfLength());
		//Zposition=Zposition + 5*(totMRD_box->GetZHalfLength()) + steelMRDplate_box->GetZHalfLength();		//offset by half total length to shift to front. why 5?
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(0);
	}
	
	// Compute dimensions. I'm not sure what &logVol i'm supposed to be giving here.. 
	void steelPlateParameterisation::ComputeDimensions (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const {
		logVol.SetXHalfLength(steelMRDplate_box->GetXHalfLength());
		logVol.SetYHalfLength(steelMRDplate_box->GetYHalfLength());
		logVol.SetZHalfLength(steelMRDplate_box->GetZHalfLength());
	}

	
// Define scintillator panels
//===========================	
	// Define functions to specify size and placement of individual scintillator panels - in header
/*	class sciPaddleParameterisation : public G4VPVParameterisation{
		public:
		sciPaddleParameterisation();						// constructor
		virtual ~sciPaddleParameterisation();					// destructor
		//G4VPhysicalVolume here is the physical volume of the instance
		virtual void ComputeTransformation					// mandatory
		  (const G4int copyNo, G4VPhysicalVolume* physVol) const;
		virtual void ComputeDimensions 						// mandatory
		  (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
		virtual G4VSolid* ComputeSolid 						// allow 2 different paddle types
		  (const G4int copyNo, G4VPhysicalVolume* physVol);			
	};
*/	

	sciPaddleParameterisation::sciPaddleParameterisation(){};			// constructor
	sciPaddleParameterisation::~sciPaddleParameterisation(){};			// destructor
	
	// Compute translation of scint panel based on panel number
	void sciPaddleParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const {
	
		Xposition=0, Yposition=0, Zposition=0;
		G4RotationMatrix* rotmtx  = new G4RotationMatrix();			// declare a rotation matrix
		G4int panelnum = floor(copyNo/numpaddlesperpanel);			// numbering from 0
/* (int) vs floor:
 (int)copyNo/30 will truncate toward zero. floor(copyNo/30) will truncate toward -infinity. This will give different values if copyNo were negative. */
		G4int paddlenum = copyNo%numpaddlesperpanel; 				// numering from 0 within a panel
		G4int pairnum = floor(paddlenum/2);					// Consider paddles 0&1 as a vertical pair; then X offset is the same for every pair
		
		// Z position is offset by (thickness of steel+scint pair) * number of preceeding pairs
		Zposition = panelnum*(steelMRDplate_box->GetZHalfLength() + sciMRDhpaddle_box->GetZHalfLength());
		Zposition=Zposition + 0.5*(totMRD_box->GetZHalfLength()) + sciMRDhpaddle_box->GetZHalfLength();		//offset by half total length to shift to front. why 0.5?
	
		// Y position is offset by half length, so that one end is at the origin. This needs to be the correct half-length
		// for the appropriate paddle (H or V type). Offsets are defined in the MOTHER ref frame. 
		// Then rotate paddles by 90deg in Y axis in every other panel. 
		// ASSUMED we'll start with a horizontal array first.
		if (panelnum%2==0){
			// horizontal panel
			if (paddlenum%2==0){
				Xposition=sciMRDhpaddle_box->GetYHalfLength(); 		// offset by +half length so one end is at x=0
			} else {
				Xposition=-sciMRDhpaddle_box->GetYHalfLength(); 	// offset by -half length so one end is at x=0
			}
			Yposition = pairnum*2*sciMRDvpaddle_box->GetXHalfLength(); 	// individual offset by pair number
			// shift whole set by 1/2 total X extent to shift center back to X=0: HalfLength cancels doubed num of paddles
			Yposition = Yposition - 0.5*(sciMRDvpaddle_box->GetXHalfLength()*numpaddlesperpanel)+sciMRDhpaddle_box->GetXHalfLength(); 
			
			
			rotmtx->rotateZ(90*deg);					// rotate horizontal panels
		} else {
			// vertical panel
			if (paddlenum%2==0){
				Yposition=sciMRDvpaddle_box->GetYHalfLength(); 
			} else {
				Yposition=-sciMRDvpaddle_box->GetYHalfLength();
			}
			// for vertical panels need to shift Y for each pair
			Xposition = pairnum*2*sciMRDvpaddle_box->GetXHalfLength(); 	// individual offset by pair number
			// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
			Xposition = Xposition - 0.5*(sciMRDvpaddle_box->GetXHalfLength()*numpaddlesperpanel)+sciMRDhpaddle_box->GetXHalfLength(); 
			
			rotmtx=0;							// don't rotate vertical panels
		}
		
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(rotmtx);

/* Rotations
The rotation associated to a physical volume represents the rotation of the reference system of this volume with respect to its mother.
A rotation matrix is constructed by instantiating the identity matrix and then applying a rotation to it. Example:

G4double phi = i*dPhi;									// define amount of rotation of element i
G4RotationMatrix rotm  = G4RotationMatrix();						// declare a rotation matrix
rotm.rotateY(90*deg); 									// rotate Y
rotm.rotateZ(phi);									// rotate Z

G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     		// declare a position 3-vector
G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;					// scale by some amount

G4Transform3D transform = G4Transform3D(rotm,position);					// combine position & rotation to form transformation

Note: transformation here includes both rotation of the box about it's origin, then translation by suitable amounts IN THE MOTHER FRAME. Therefore the translation is by quantitites calculated from the rotation. 
*/

	}
	
	// Compute dimensions based on copyNo
	void sciPaddleParameterisation::ComputeDimensions (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const {
		G4int panelnum = floor(copyNo/numpaddlesperpanel);			// numbering from 0
		// ASSUMED we'll start with a horizontal array first.
		if (panelnum%2==0){
			logVol.SetYHalfLength(sciMRDhpaddle_box->GetYHalfLength());
		} else {
			logVol.SetYHalfLength(sciMRDvpaddle_box->GetYHalfLength());
		}
		// Other dims the same
		logVol.SetXHalfLength(sciMRDvpaddle_box->GetXHalfLength());
		logVol.SetZHalfLength(sciMRDvpaddle_box->GetZHalfLength());
	}
	
	// Do i need a ComputeSolid function, as my horiz paddles and vert paddles are different solids? 
	// But they only differ by dimensions (both G4Box) and dimensions are specified in ComputeDimensions... 
	
// Define both combined
//=====================
// So you can't add more than one parameterisation volume to a single mother. Since the two parameterisations overlap, I can't use two mothers either. 
// Both must be combined into one unweildy parameterisation. 
	// Define functions to specify size and placement of combined scintillator paddles and steel plates - in header
/*	class mrdParamaterisedModel : public G4VPVParameterisation{
		public:
		mrdParamaterisedModel();						// constructor
		virtual ~mrdParamaterisedModel();					// destructor
		//G4VPhysicalVolume here is the physical volume of the instance
		virtual void ComputeTransformation					// mandatory
		  (const G4int copyNo, G4VPhysicalVolume* physVol) const;
		virtual void ComputeDimensions 						// mandatory
		  (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
		virtual G4VSolid* ComputeSolid 						// allow 2 different paddle types
		  (const G4int copyNo, G4VPhysicalVolume* physVol);			
	};
*/	

	mrdParamaterisedModel::mrdParamaterisedModel(){};				// constructor
	mrdParamaterisedModel::~mrdParamaterisedModel(){};				// destructor
	
	// Compute translation of scint panel based on panel number
	void mrdParamaterisedModel::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const {
	
		Xposition=0, Yposition=0, Zposition=0; 
		G4RotationMatrix* rotmtx  = new G4RotationMatrix();			// declare a rotation matrix
		if (copyNo<(numpaddlesperpanel*numpanels)) {	
			// COPY IS A PADDLE
			//=================
			G4int panelnum = floor(copyNo/numpaddlesperpanel);			// numbering from 0
			G4int paddlenum = copyNo%numpaddlesperpanel; 				// numering from 0 within a panel
			G4int pairnum = floor(paddlenum/2);					// Consider paddles 0&1 as a vertical pair; then X offset is the same for every pair
		
			// Z position is offset by (thickness of steel+scint pair) * number of preceeding pairs
			Zposition = panelnum*2*(steelMRDplate_box->GetZHalfLength() + sciMRDhpaddle_box->GetZHalfLength() + aluMRDstruc_box->GetZHalfLength());
			//offset by half total z depth to shift to front, then offset by depth of 1 alu struct.
			Zposition=Zposition + sciMRDhpaddle_box->GetZHalfLength() + 2*aluMRDstruc_box->GetZHalfLength();		
	
			// Y position is offset by half length, so that one end is at the origin. This needs to be the correct half-length
			// for the appropriate paddle (H or V type). Offsets are defined in the MOTHER ref frame. 
			// Then rotate paddles by 90deg in Y axis in every other panel. 
			// ASSUMED we'll start with a horizontal array first.
			if (panelnum%2==0){
				// horizontal panel
				if (paddlenum%2==0){
					Xposition=sciMRDhpaddle_box->GetYHalfLength(); 		// offset by +half length so one end is at x=0
				} else {
					Xposition=-sciMRDhpaddle_box->GetYHalfLength(); 	// offset by -half length so one end is at x=0
				}
				Yposition = pairnum*2*sciMRDvpaddle_box->GetXHalfLength(); 	// individual offset by pair number
				// shift whole set by 1/2 total X extent to shift center back to X=0: HalfLength cancels doubed num of paddles
				Yposition = Yposition - 0.5*(sciMRDvpaddle_box->GetXHalfLength()*numpaddlesperpanel)+sciMRDhpaddle_box->GetXHalfLength(); 
			
				rotmtx->rotateZ(90*deg);					// rotate horizontal panels
			} else {
				// vertical panel
				if (paddlenum%2==0){
					Yposition=sciMRDvpaddle_box->GetYHalfLength(); 
				} else {
					Yposition=-sciMRDvpaddle_box->GetYHalfLength();
				}
				// for vertical panels need to shift Y for each pair
				Xposition = pairnum*2*sciMRDvpaddle_box->GetXHalfLength(); 	// individual offset by pair number
				// shift whole set by 1/2 total Y extent to shift center back to Y=0: HalfLength cancels doubed num of paddles
				Xposition = Xposition - 0.5*(sciMRDvpaddle_box->GetXHalfLength()*numpaddlesperpanel)+sciMRDhpaddle_box->GetXHalfLength(); 
			
				rotmtx=0;							// don't rotate vertical panels
			}

			//=================
		} else {
			// COPY IS A STEEL PLATE
			//======================
			G4int plateNo=copyNo-(numpaddlesperpanel*numpanels);
			Zposition=(plateNo)*2*(steelMRDplate_box->GetZHalfLength() + sciMRDhpaddle_box->GetZHalfLength() + aluMRDstruc_box->GetZHalfLength());
			Zposition=Zposition + steelMRDplate_box->GetZHalfLength() + 2*(sciMRDhpaddle_box->GetZHalfLength()+aluMRDstruc_box->GetZHalfLength());
			//offset by half depth to shift to front. Also offset by paddle depth because we don't want both in the same place.
			Xposition=0, Yposition=0;
			rotmtx=0;
			//======================
		}
		Zposition=Zposition - (totMRD_box->GetZHalfLength());		//offset by half total length to shift to front.
		G4ThreeVector origin(Xposition,Yposition,Zposition);
		physVol->SetTranslation(origin);
		physVol->SetRotation(rotmtx);
	}
	
	// Compute dimensions based on copyNo
	void mrdParamaterisedModel::ComputeDimensions (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const {
		if (copyNo<(numpaddlesperpanel*numpanels)) {
			// COPY IS A PADDLE
			//=================
			G4int panelnum = floor(copyNo/numpaddlesperpanel);		// numbering from 0
			// ASSUMED we'll start with a horizontal array first.
			if (panelnum%2==0){
				logVol.SetYHalfLength(sciMRDhpaddle_box->GetYHalfLength());
			} else {
				logVol.SetYHalfLength(sciMRDvpaddle_box->GetYHalfLength());
			}
			// Other dims the same
			logVol.SetXHalfLength(sciMRDvpaddle_box->GetXHalfLength());
			logVol.SetZHalfLength(sciMRDvpaddle_box->GetZHalfLength());
			//=================
		} else {
			// COPY IS A STEEL PLATE
			//======================
			logVol.SetXHalfLength(steelMRDplate_box->GetXHalfLength());
			logVol.SetYHalfLength(steelMRDplate_box->GetYHalfLength());
			logVol.SetZHalfLength(steelMRDplate_box->GetZHalfLength());
			//======================
		}
	}

	// Compute material, sensitivity, visAttributes based on copyNo
        G4Material* mrdParamaterisedModel::ComputeMaterial (const G4int copyNo, G4VPhysicalVolume* physVol, const G4VTouchable * parentTouch) {
		G4Material* material;
		G4VisAttributes* myVisAtt = new G4VisAttributes; 
  		myVisAtt->SetVisibility(true);
		if (copyNo<(numpaddlesperpanel*numpanels)){
			material = G4Material::GetMaterial("Scinti");
			myVisAtt->SetColour(1.0,0.0,1.0);
			
		} else {
			material = G4Material::GetMaterial("Steel");
			myVisAtt->SetColour(1.0,1.0,0.0);
		}
		physVol->GetLogicalVolume()->SetVisAttributes(myVisAtt);	//can set visualisation attributes like this
		return material;
	}
// Define Alu Support Structure
//=============================
// Couldn't do this with a parameterisation because it requires subtraction solid (parameterisation only works with CGS solids)
// Instead combine them into an assembly, then make imprints for each layer.

void makeAlu(G4AssemblyVolume* aluMRDassembly){
	// this code must be in a function as it uses for loops. :| 

	//N.B. subtraction solids can only use CGS solids or the output of boolean operations, not assemblies or other stuff
	G4RotationMatrix  Ra (0,0,0);
	G4ThreeVector  Ta (0,0,0);
	Ta.set(windowwidth+alufullxthickness,windowheight+alufullythickness,0);
	G4Transform3D transform = G4Transform3D(Ra,Ta);
	G4SubtractionSolid* aluMRDstruc_sol = new G4SubtractionSolid("aluStruct",aluMRDstruc_box,aluMRDwindow_box,transform);
	G4SubtractionSolid* aluMRDstruc_sol2 = aluMRDstruc_sol;
	aluMRDstruc_sol = aluMRDstruc_sol2;
	
	for (G4int row=-1;row<2;row++){
		for (G4int col=-1;col<2;col++){
			G4double xoffset=col*(alufullxthickness+windowwidth);
			G4double yoffset=row*(alufullythickness+windowheight);
			Ta.set(xoffset,yoffset,0);                        
			transform = G4Transform3D(Ra,Ta);
// using method with G4Transform3D means we can change the transform after as it is passed byval. If we used a passive transform we would need to maintain the G4Transform3D. 
			aluMRDstruc_sol2=0;
			delete aluMRDstruc_sol2;
			G4SubtractionSolid* aluMRDstruc_sol2
			 = new G4SubtractionSolid("aluStruct",aluMRDstruc_sol,aluMRDwindow_box,transform);
			aluMRDstruc_sol = aluMRDstruc_sol2;
		}
	}

	G4LogicalVolume* aluMRDstruc_log 
	= new G4LogicalVolume(aluMRDstruc_sol, G4Material::GetMaterial("Aluminum"), "aluStruct",0,0,0) ;

	Ra.set(0,0,0);
	Ta.set(0,0,0);
	G4double layeroffset = 2*(aluMRDstruc_box->GetZHalfLength()+sciMRDhpaddle_box->GetZHalfLength()+steelMRDplate_box->GetZHalfLength());
	for (G4int structnum=0;structnum<numalustructs;structnum++){
		G4double zpos = structnum*layeroffset+aluMRDstruc_box->GetZHalfLength()- (totMRD_box->GetZHalfLength());
		Ta.set(0,0,zpos); //offset by half total length to shift to front.
		aluMRDassembly->AddPlacedVolume(aluMRDstruc_log,Ta,&Ra);
	}
	//aluMRDassembly->MakeImprint(expHall->GetLogicalVolume(), Tm,&Rm);  //placement <- in DetectorConstruction
}
