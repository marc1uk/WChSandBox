// ====================================================================
//   SBsimMRDDB.cc
//   
//   2006/09/04 K. Hiraide
// ====================================================================
#include "SBsimMRDDB.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

////////////////////////
SBsimMRDDB::SBsimMRDDB()
////////////////////////
{
  //SetModuleInfo((char*)"database/mrdgeom/mrdmodule.txt");
  SetModuleInfo((char*)"../WChSandBox_v1/src/mrdmodule.txt");
  SetPositionInfo((char*)"../WChSandBox_v1/src/mrdposition.txt");
  SetAlignmentInfo("../WChSandBox_v1/src/mrd-align.txt");
  SetNoiseHitInfo((char*)"../WChSandBox_v1/src/mrdnoisehit.txt");
}

/////////////////////////
SBsimMRDDB::~SBsimMRDDB()
/////////////////////////
{

}

///////////////////////////////////////////
void SBsimMRDDB::SetModuleInfo(char* fname)
///////////////////////////////////////////
{
  cout << fname << endl;
  ifstream ifs(fname);
  if ( ifs.fail() ) {
    cerr << "SBsimMRDDB::SetModuleInfo: file open error!" << endl;
    exit(-1);
  }

  char buf[256];
  string name;

  while ( !ifs.eof() ) {
    ifs.getline(buf,256);
    istringstream iss(buf);
    iss >> name;

    if      (name=="NLayer")      iss >> mod.NLayer;
    else if (name=="IronSizeXY")    iss >> mod.IronSizeXY[0]
                                      >> mod.IronSizeXY[1];
    else if (name=="IronSizeZ") {
      if (mod.NLayer > 12){
        cout << "Error: too many NLayer specified. Force set to 12" << endl;
        mod.NLayer = 12;
      }
      for (int iLayer = 0; iLayer < mod.NLayer ; iLayer ++){
        iss >> mod.IronSizeZ[iLayer];
        mod.IronSizeZ[iLayer] *= 0.5; // change to half size
      }
    }
    else if (name=="VScintiSize") iss >> mod.VScintiSize[0]
				      >> mod.VScintiSize[1]
				      >> mod.VScintiSize[2];
    else if (name=="HScintiSize") iss >> mod.HScintiSize[0]
				      >> mod.HScintiSize[1]
				      >> mod.HScintiSize[2];
    else if (name=="TScintiSize") iss >> mod.TScintiSize[0]
				      >> mod.TScintiSize[1]
				      >> mod.TScintiSize[2]
				      >> mod.TScintiSize[3]
				      >> mod.TScintiSize[4];
    else if (name=="VetoVSize")   iss >> mod.VetoVSize[0]
				      >> mod.VetoVSize[1]
				      >> mod.VetoVSize[2];
    else if (name=="VetoHSize")   iss >> mod.VetoHSize[0]
				      >> mod.VetoHSize[1]
				      >> mod.VetoHSize[2];
    else if (name=="VetoESize")   iss >> mod.VetoESize[0]
				      >> mod.VetoESize[1]
				      >> mod.VetoESize[2];
    else if (name=="ScintiGap")   iss >> mod.ScintiGap;
    else if (name=="LayerGap")    iss >> mod.LayerGap;
    else if (name=="IronScintiGap")    iss >> mod.IronScintiGap;
    else if (name=="AvgEdep")    iss >> mod.AvgEdep;
    else if (name=="Attlength")   iss >> mod.Attlength[0]
				      >> mod.Attlength[1];
    else if (name=="InitIntensity") iss >> mod.InitIntensity[0]
				        >> mod.InitIntensity[1];
    else if (name=="LGSize") iss >> mod.LGSize[0]
                                 >> mod.LGSize[1]
                                 >> mod.LGSize[2]
                                 >> mod.LGSize[3]
                                 >> mod.LGSize[4];
    else if (name=="AlSizeV1") iss >> mod.AlSizeV1[0]
                                   >> mod.AlSizeV1[1]
                                   >> mod.AlSizeV1[2]
                                   >> mod.AlSizeV1[3];
    else if (name=="AlSizeV2") iss >> mod.AlSizeV2[0]
                                   >> mod.AlSizeV2[1]
                                   >> mod.AlSizeV2[2]
                                   >> mod.AlSizeV2[3];
    else if (name=="AlSizeV3") iss >> mod.AlSizeV3[0]
                                   >> mod.AlSizeV3[1]
                                   >> mod.AlSizeV3[2]
                                   >> mod.AlSizeV3[3];
    else if (name=="AlSizeV4") iss >> mod.AlSizeV4[0]
                                   >> mod.AlSizeV4[1]
                                   >> mod.AlSizeV4[2]
                                   >> mod.AlSizeV4[3];
    else if (name=="AlSizeV5") iss >> mod.AlSizeV5[0]
                                   >> mod.AlSizeV5[1]
                                   >> mod.AlSizeV5[2]
                                   >> mod.AlSizeV5[3];
    else if (name=="AlSizeH1") iss >> mod.AlSizeH1[0]
                                   >> mod.AlSizeH1[1]
                                   >> mod.AlSizeH1[2]
                                   >> mod.AlSizeH1[3];
    else if (name=="AlSizeH2") iss >> mod.AlSizeH2[0]
                                   >> mod.AlSizeH2[1]
                                   >> mod.AlSizeH2[2]
                                   >> mod.AlSizeH2[3];
    else if (name=="AlSizeH3") iss >> mod.AlSizeH3[0]
                                   >> mod.AlSizeH3[1]
                                   >> mod.AlSizeH3[2]
                                   >> mod.AlSizeH3[3];
  }
}

////////////////////////////////////////////
void SBsimMRDDB::SetPositionInfo(char* fname)
////////////////////////////////////////////
{
  cout << fname << endl;
  ifstream ifs(fname);
  if ( ifs.fail() ) {
    cerr << "SBsimMRDDB::SetPositionInfo: file open error!" << endl;
    exit(-1);
  }

  char buf[256];
  string name;

  while ( !ifs.eof() ) {
    ifs.getline(buf,256);
    istringstream iss(buf);
    iss >> name;

    if (name=="GlobalPos") iss >> pos.GlobalPosition[0]
			       >> pos.GlobalPosition[1]
			       >> pos.GlobalPosition[2];
  }
}

//////////////////////////////////////////////
void SBsimMRDDB::SetAlignmentInfo(char* fname)
//////////////////////////////////////////////
{
  cout << fname << endl;
  ifstream ifs(fname);
  if ( ifs.fail() ) {
    cerr << "SBsimMRDDB::SetAlignmentInfo: file open error!" << endl;
    exit(-1);
  }

  char buf[256];
  int layer=0;

  while ( !ifs.eof() ) {
    ifs.getline(buf,256);
    if (buf[0]!='#') {
      istringstream iss(buf);
      iss >> layer;
      iss >> pos.PlanePosition[layer][0]
          >> pos.PlanePosition[layer][1]
          >> pos.PlanePosition[layer][2]
          >> pos.PlaneRotation[layer][0]
          >> pos.PlaneRotation[layer][1]
          >> pos.PlaneRotation[layer][2];
    }
  }
}

////////////////////////////////////////////
void SBsimMRDDB::SetNoiseHitInfo(char* fname)
////////////////////////////////////////////
{
  cout << fname << endl;
  ifstream ifs(fname);
  if ( ifs.fail() ) {
    cerr << "SBsimMRDDB::SetNoiseHitInfo: file open error!" << endl;
    exit(-1);
  }

  char buf[256];
  string name;
  int pmtnum = 0;
  while ( !ifs.eof() ) {
    ifs.getline(buf,256);
    istringstream iss(buf);

    iss >> noise.ModuleNum[pmtnum] >> noise.NoiseRate[pmtnum] ;    
    pmtnum++;    
  }
}

//////////////////////////
void SBsimMRDDB::PrintAll()
//////////////////////////
{
  PrintModuleInfo();
  PrintPositionInfo();

  for (int i=0;i<13;i++) {
    cout << pos.PlanePosition[i][0] << " "
	 << pos.PlanePosition[i][1] << " "
	 << pos.PlanePosition[i][2] << endl;
  }
}

/////////////////////////////////
void SBsimMRDDB::PrintModuleInfo()
/////////////////////////////////
{
  cout << "##### MRD Module Variables #####" << endl;
  cout << "NLayer      : " << mod.NLayer     << endl;
  cout << "IronSizeXY    : ("
       << mod.IronSizeXY[0] << ", "
       << mod.IronSizeXY[1] << ")" 
       << endl;
  cout << "IronSizeZ    : (";
  for (int iLayer = 0; iLayer < mod.NLayer ; iLayer ++){
    cout << mod.IronSizeZ[iLayer] << ", ";
  }
  cout << ")" << endl;
  cout << "VScintiSize : ("
       << mod.VScintiSize[0] << ", "
       << mod.VScintiSize[1] << ", "
       << mod.VScintiSize[2] << ")" 
       << endl;
  cout << "HScintiSize : ("
       << mod.HScintiSize[0] << ", "
       << mod.HScintiSize[1] << ", "
       << mod.HScintiSize[2] << ")" 
       << endl;
  cout << "HScintiSize : ("
       << mod.TScintiSize[0] << ", "
       << mod.TScintiSize[1] << ", "
       << mod.TScintiSize[2] << ", " 
       << mod.TScintiSize[3] << ", "
       << mod.TScintiSize[4] << ")" 
       << endl;
  cout << "VetoVSize   : ("
       << mod.VetoVSize[0] << ", "
       << mod.VetoVSize[1] << ", "
       << mod.VetoVSize[2] << ")" 
       << endl;
  cout << "VetoHSize   : ("
       << mod.VetoHSize[0] << ", "
       << mod.VetoHSize[1] << ", "
       << mod.VetoHSize[2] << ")" 
       << endl;
  cout << "VetoESize   : ("
       << mod.VetoESize[0] << ", "
       << mod.VetoESize[1] << ", "
       << mod.VetoESize[2] << ")" 
       << endl;
  cout << "ScintiGap   : " << mod.ScintiGap  << endl;
  cout << "LayerGap    : " << mod.LayerGap   << endl;
  cout << "IronSintiGap    : " << mod.IronScintiGap   << endl;
  cout << "AvgEdep   : " << mod.AvgEdep   << endl;
  cout << "Attlength   : ("
       << mod.Attlength[0] << ", "
       << mod.Attlength[1] << ") "
       << endl;
  cout << "InitIntensity   : ("
       << mod.InitIntensity[0] << ", "
       << mod.InitIntensity[1] << ") "
       << endl;
  cout << "LGSize   : ("
       << mod.LGSize[0] << ", "
       << mod.LGSize[1] << ", "
       << mod.LGSize[2] << ", " 
       << mod.LGSize[3] << ", "
       << mod.LGSize[4] << ")" 
       << endl;
  cout << endl;
}

///////////////////////////////////
void SBsimMRDDB::PrintPositionInfo()
///////////////////////////////////
{
  cout << "##### MRD Position Variables #####" << endl;
  cout << "GlobalPos : ("
       << pos.GlobalPosition[0] << ", "
       << pos.GlobalPosition[1] << ", "
       << pos.GlobalPosition[2] << ")" 
       << endl;
  cout << endl;
}
