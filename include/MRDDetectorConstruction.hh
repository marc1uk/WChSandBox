

#ifndef MRDDetectorConstruction_h
#define MRDDetectorConstruction_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4AssemblyVolume.hh"
	
class steelPlateParameterisation : public G4VPVParameterisation { 
  public:
    steelPlateParameterisation();
    virtual ~steelPlateParameterisation();
    void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;
    void ComputeDimensions (G4Box & trackerLayer, const G4int copyNo,
                            const G4VPhysicalVolume* physVol) const;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class sciPaddleParameterisation : public G4VPVParameterisation { 
  public:
    sciPaddleParameterisation();
    virtual ~sciPaddleParameterisation();
    void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;
    void ComputeDimensions (G4Box & trackerLayer, const G4int copyNo,
                            const G4VPhysicalVolume* physVol) const;
};

class mrdParamaterisedModel : public G4VPVParameterisation{

	public:
	mrdParamaterisedModel();						// constructor
	virtual ~mrdParamaterisedModel();					// destructor
	//G4VPhysicalVolume here is the physical volume of the instance
	virtual void ComputeTransformation					// mandatory
	  (const G4int copyNo, G4VPhysicalVolume* physVol) const;
	virtual void ComputeDimensions 						// mandatory
	  (G4Box& logVol, const G4int copyNo, const G4VPhysicalVolume* physVol) const;
	//virtual G4VSolid* ComputeSolid 					// allow 2 different paddle types
	//  (const G4int copyNo, G4VPhysicalVolume* physVol);
        virtual G4Material* ComputeMaterial (const G4int copyNo , G4VPhysicalVolume* physVol , const G4VTouchable * parentTouch =0);
        									// don't use touchable unless you're using nested parameterisation (just leave it null)
};
#endif

void makeAlu(G4AssemblyVolume* aluMRDassembly);
