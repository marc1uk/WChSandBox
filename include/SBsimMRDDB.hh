// ====================================================================
//   SBsimMRDDB.hh
//
//   2006/09/04 K. Hiraide
// ====================================================================
#ifndef SBSIM_MRD_DB
#define SBSIM_MRD_DB

struct MRDModule{
  int NLayer;
  double IronSizeXY[2];
  double IronSizeZ[12];
  double VScintiSize[3], HScintiSize[3], TScintiSize[5];
  double VetoVSize[3], VetoHSize[3], VetoESize[3];
  double ScintiGap;
  double LayerGap;
  double IronScintiGap;
  double AvgEdep;
  double Attlength[2];
  double InitIntensity[2];
  double LGSize[5];
  double AlSizeV1[4];
  double AlSizeV2[4];
  double AlSizeV3[4];
  double AlSizeV4[4];
  double AlSizeV5[4];
  double AlSizeH1[4];
  double AlSizeH2[4];
  double AlSizeH3[4];
};

struct MRDPosition {
  double GlobalPosition[3];
  double PlanePosition[13][3];
  double PlaneRotation[13][3];
  double IronPosition[13][3];
  double IronRotation[13][3]; 
};

struct MRDNoiseHit {
  int ModuleNum[512];
  double NoiseRate[512];
};

class SBsimMRDDB {
public:
  SBsimMRDDB();
  ~SBsimMRDDB();

  void SetModuleInfo(char* fname);
  void SetPositionInfo(char* fname);
  void SetAlignmentInfo(char* fname);
  void SetNoiseHitInfo(char* fname);

  void PrintAll();
  void PrintModuleInfo();
  void PrintPositionInfo();

  inline MRDModule*   GetMRDModuleInfo(){return &mod;};
  inline MRDPosition* GetMRDPositionInfo(){return &pos;};
  inline MRDNoiseHit* GetMRDNoiseHitInfo(){return &noise;};

private:
  MRDModule   mod;
  MRDPosition pos;
  MRDNoiseHit noise;

};

#endif
