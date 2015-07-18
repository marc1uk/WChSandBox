#ifndef EVENTDISPLAY3D_HH
#define EVENTDISPLAY3D_HH

#include "TObject.h"
#include "TH1.h"

#include <iostream>
#include <cmath>

#include "TGeometry.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"

#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveViewer.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TEveStraightLineSet.h"
#include "TEveArrow.h"
#include "TEveText.h"



class EventDisplay3D : public TObject {

 public: 

  EventDisplay3D();

  ~EventDisplay3D();

  void AddHits(std::vector<std::vector<double> > hits) {_hits = hits;}

  void DoGen(std::vector<std::vector<double> > gcoor) {_genpositions = gcoor; _drawgen=true;}

  void DrawDetector();

  void PlotEvent();

  void SetTimeRanges(int ri, double tt) {_tb[ri]=tt;}


 private: 

  bool _drawgen;

  std::vector<std::vector<double> > _hits;
  std::vector<std::vector<double> > _genpositions;
 
  double _tb[8];
 
  TGeoManager *geom;	
  TGeoMaterial *matVacuum;
  TGeoMaterial *matWater;
  TGeoMedium* Vacuum;
  TGeoMedium* Water;     

  TGeoVolume *top; 


  
  ClassDef(EventDisplay3D,0)

};

#endif
