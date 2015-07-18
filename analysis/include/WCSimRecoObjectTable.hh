#ifndef WCSIMRECOOBJECTTABLE_HH
#define WCSIMRECOOBJECTTABLE_HH

#include "TObject.h"

class WCSimRecoObjectTable : public TObject {

 public:
  static WCSimRecoObjectTable* Instance();

  void NewDigit(){ numDigits++; }
  void DeleteDigit(){ numDigits--; }
  Int_t NumberOfDigits(){ return numDigits; }

  void NewCluster() { numClusters++; }
  void DeleteCluster(){ numClusters--; }
  Int_t NumberOfClusters(){ return numClusters; }

  void NewClusterDigit(){ numClusterDigits++; }
  void DeleteClusterDigit(){ numClusterDigits--; }
  Int_t NumberOfClusterDigits(){ return numClusterDigits; }

  void NewVertex(){ numVertices++; }
  void DeleteVertex(){ numVertices--; }
  Int_t NumberOfVertices(){ return numVertices; }

  void NewRing(){ numRings++; }
  void DeleteRing(){ numRings--; }
  Int_t NumberOfRings(){ return numRings; }

  void NewEvent(){ numEvents++; }
  void DeleteEvent(){ numEvents--; }
  Int_t NumberOfEvents(){ return numEvents; }

  void Reset();
  void Print();

 private:
  WCSimRecoObjectTable();
  ~WCSimRecoObjectTable();

  Int_t numDigits;
  Int_t numClusters;
  Int_t numClusterDigits;
  Int_t numVertices;
  Int_t numRings;
  Int_t numEvents;

  ClassDef(WCSimRecoObjectTable,0)

};

#endif







