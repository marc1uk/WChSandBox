#ifndef WCSIMVERTEXGEOMETRY_HH
#define WCSIMVERTEXGEOMETRY_HH

#include "TObject.h"

#include <vector>

class WCSimVertexGeometry : public TObject {

 public:
  static WCSimVertexGeometry* Instance();

  //  void LoadEvent(WCSimRecoEvent* evt);

  void CalcResiduals(Double_t vx, Double_t vy, Double_t vz, Double_t vtxTime,
		      Double_t px, Double_t py, Double_t pz);

  void CalcTResid(Double_t vx, Double_t vy, Double_t vz, Double_t vtxTime,
		      Double_t px, Double_t py, Double_t pz, Double_t hx, Double_t hy, Double_t hz, Double_t ht, Double_t &fTResid, Double_t &Lphoton, Double_t &Ltrack);

  void CalcTResidColor(Double_t vx, Double_t vy, Double_t vz, Double_t vtxTime,
		      Double_t px, Double_t py, Double_t pz, Double_t hx, Double_t hy, Double_t hz, Double_t ht, Double_t flmda, Double_t fvmuon, Double_t &fTResid, Double_t &Lphoton, Double_t &Ltrack);

  void CalcPointResiduals(Double_t vx, Double_t vy, Double_t vz, Double_t vtxTime,
		           Double_t px, Double_t py, Double_t pz);

  void CalcExtendedResiduals(Double_t vx, Double_t vy, Double_t vz, Double_t vtxTime,
		              Double_t px, Double_t py, Double_t pz);

  //  void CalcVertexSeeds(WCSimRecoEvent* evt, Int_t NSeeds = 1);

  void CalcVertexSeeds(Int_t NSeeds = 1);

//  WCSimRecoVertex* CalcSimpleVertex(WCSimRecoEvent* evt);

//  WCSimRecoVertex* CalcSimpleVertex();

//  WCSimRecoVertex* CalcSimpleDirection(WCSimRecoEvent* evt, WCSimRecoVertex* vtx);

//  WCSimRecoVertex* CalcSimpleDirection(WCSimRecoVertex* vtx);

  Int_t GetNDigits() { return fNDigits; }

  Int_t GetNFilterDigits() { return fNFilterDigits; }

  Bool_t IsFiltered(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fIsFiltered[idigit]; 
    }
    else return 0;
  }

  Double_t GetDigitX(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDigitX[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDigitY(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDigitY[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDigitZ(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDigitZ[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDigitT(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDigitT[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDigitQ(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDigitQ[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetConeAngle(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fConeAngle[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetAngle(Int_t idigit) { 
    return GetZenith(idigit);
  }

  Double_t GetZenith(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fZenith[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetAzimuth(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fAzimuth[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetSolidAngle(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fSolidAngle[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDistPoint(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDistPoint[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDistTrack(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDistTrack[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDistPhoton(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDistPhoton[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDistScatter(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDistScatter[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaTime(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaTime[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaSigma(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaSigma[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaAngle(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaAngle[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaPoint(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaPoint[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaTrack(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaTrack[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaPhoton(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaPhoton[idigit]; 
    }
    else return 0.0;
  }  

  Double_t GetDeltaScatter(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDeltaScatter[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetPointPath(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fPointPath[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetExtendedPath(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fExtendedPath[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetPointResidual(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fPointResidual[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetExtendedResidual(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fExtendedResidual[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDelta(Int_t idigit) { 
    if( idigit>=0 && idigit<fNDigits ) {
      return fDelta[idigit]; 
    }
    else return 0.0;
  }

  Double_t GetDeltaCorrection(Int_t idigit, Double_t length);

  UInt_t GetNSeeds() {
    return vSeedVtxX.size();
  }

  Double_t GetSeedVtxX(UInt_t ivtx) {
    if( ivtx>=0 && ivtx<GetNSeeds() ){
      return vSeedVtxX.at(ivtx);
    }
    else return 0.0;
  }

  Double_t GetSeedVtxY(UInt_t ivtx) {
    if( ivtx>=0 && ivtx<GetNSeeds() ){
      return vSeedVtxY.at(ivtx);
    }
    else return 0.0;
  }

  Double_t GetSeedVtxZ(UInt_t ivtx) {
    if( ivtx>=0 && ivtx<GetNSeeds() ){
      return vSeedVtxZ.at(ivtx);
    }
    else return 0.0;
  } 

  Double_t GetSeedVtxTime(UInt_t ivtx) {
    if( ivtx>=0 && ivtx<GetNSeeds() ){
      return vSeedVtxTime.at(ivtx);
    }
    else return 0.0;
  }

 private:
  WCSimVertexGeometry();
  ~WCSimVertexGeometry();

  void CalcSimpleVertex(Double_t& vtxX, Double_t& vtxY, Double_t& vtxZ, 
                        Double_t& vtxTime);

  void ChooseNextQuadruple(Double_t& x0, Double_t& y0, Double_t& z0, Double_t& t0,
                           Double_t& x1, Double_t& y1, Double_t& z1, Double_t& t1,
                           Double_t& x2, Double_t& y2, Double_t& z2, Double_t& t2,
                           Double_t& x3, Double_t& y3, Double_t& z3, Double_t& t3);

  void ChooseNextDigit(Double_t& x, Double_t& y, Double_t& z, Double_t& t);

  Int_t fRunNum;
  Int_t fEventNum;
  Int_t fTriggerNum;

  Int_t fPMTs;
  Int_t fNDigits;
  Int_t fNFilterDigits;

  Int_t fThisDigit;
  Int_t fLastEntry;
  Int_t fCounter;

  Double_t fVtxX1;
  Double_t fVtxY1;
  Double_t fVtxZ1;
  Double_t fVtxTime1;

  Double_t fVtxX2;
  Double_t fVtxY2;
  Double_t fVtxZ2;
  Double_t fVtxTime2;

  Double_t fMeanQ;
  Double_t fTotalQ; 

  Double_t fMeanFilteredQ;
  Double_t fTotalFilteredQ;

  Double_t fMinTime;
  Double_t fMaxTime;
  
  Double_t fQmin;

  Bool_t* fIsFiltered; 

  Double_t* fDigitX;           // Digit X (cm)
  Double_t* fDigitY;           // Digit Y (cm)
  Double_t* fDigitZ;           // Digit Z (cm)
  Double_t* fDigitT;           // Digit T (ns)
  Double_t* fDigitQ;           // Digit Q (PEs)
  Int_t* fDigitPE;

  Double_t* fZenith;           // Zenith (degrees)
  Double_t* fAzimuth;          // Azimuth (degrees)
  Double_t* fSolidAngle;       // SolidAngle = sin(angle)/sin(42)
  Double_t* fConeAngle;        // Cone Angle (degrees)

  Double_t* fDistPoint;        // Distance from Vertex (S)
  Double_t* fDistTrack;        // Distance along Track (S1)
  Double_t* fDistPhoton;       // Distance along Photon (S2)
  Double_t* fDistScatter;      // Distance of Scattering 

  Double_t* fDeltaTime;        // DeltaT = T - VtxT
  Double_t* fDeltaSigma;       // SigmaT(Q)

  Double_t* fDeltaAngle;       // DeltaAngle  = ConeAngle - Angle
  Double_t* fDeltaPoint;       // DeltaPoint  = DistPoint/(fC/fN)
  Double_t* fDeltaTrack;       // DeltaTrack  = DistTrack/(fC)
  Double_t* fDeltaPhoton;      // DeltaPhoton = DistPhoton/(fC/fN)
  Double_t* fDeltaScatter;     // DeltaScatter = DeltaScatter/(fC/fN)

  Double_t* fPointPath;        // PointPath = fN*DistPoint
  Double_t* fExtendedPath;     // ExtendedPath = DistTrack + fN*DistPhoton

  Double_t* fPointResidual;    // PointResidual = DeltaTime - PointPath
  Double_t* fExtendedResidual; // ExtendedResidual = DeltaTime - ExtendedPath
  Double_t* fDelta;            // Chosen Residual [Point/Extended/Null]
  

  std::vector<Double_t> vSeedVtxX;         
  std::vector<Double_t> vSeedVtxY;
  std::vector<Double_t> vSeedVtxZ;
  std::vector<Double_t> vSeedVtxTime;
  std::vector<Int_t> vSeedDigitList;

  //  std::vector<WCSimRecoVertex*> vVertexList;

  ClassDef(WCSimVertexGeometry,0)

};

#endif
