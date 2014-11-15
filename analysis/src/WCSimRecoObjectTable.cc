
#include "WCSimRecoObjectTable.hh"

#include <iostream>
#include <cassert>

ClassImp(WCSimRecoObjectTable)

static WCSimRecoObjectTable* fgRecoObjectTable = 0;

WCSimRecoObjectTable* WCSimRecoObjectTable::Instance()
{
  if( !fgRecoObjectTable ){
    fgRecoObjectTable = new WCSimRecoObjectTable();
  }

  if( !fgRecoObjectTable ){
    assert(fgRecoObjectTable);
  }

  if( fgRecoObjectTable ){

  }

  return fgRecoObjectTable;
}

WCSimRecoObjectTable::WCSimRecoObjectTable()
{
  this->Reset();
}

WCSimRecoObjectTable::~WCSimRecoObjectTable()
{
  
}

void WCSimRecoObjectTable::Reset()
{
  numDigits = 0;
  numClusters = 0;
  numClusterDigits = 0;
  numVertices = 0;
  numRings = 0;
  numEvents = 0;
}

void WCSimRecoObjectTable::Print()
{
  std::cout << " *** WCSimRecoObjectTable::Print() *** " << std::endl;
  std::cout << numDigits << "\t Digits " << std::endl;
  std::cout << numClusterDigits << "\t ClusterDigits " << std::endl;
  std::cout << numClusters << "\t Clusters " << std::endl;
  std::cout << numVertices << "\t Vertices " << std::endl;
  std::cout << numRings << "\t Rings " << std::endl;
  std::cout << numEvents << "\t Events " << std::endl;
}
