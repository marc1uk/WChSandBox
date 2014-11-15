#include "WCSimRecoVertex.hh"

#include "WCSimRecoObjectTable.hh"
#include "WCSimParameters.hh"

ClassImp(WCSimRecoVertex)

WCSimRecoVertex::WCSimRecoVertex()
{
  this->Reset();

  WCSimRecoObjectTable::Instance()->NewVertex();
}

WCSimRecoVertex::WCSimRecoVertex( Double_t x, Double_t y, Double_t z )
{
  this->Reset();
  
  this->SetVertex(x,y,z);
  this->SetFOM(0.0,1,1);

  WCSimRecoObjectTable::Instance()->NewVertex();
}
  
WCSimRecoVertex::WCSimRecoVertex( Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz )
{
  this->Reset();

  this->SetVertex(x,y,z);
  this->SetDirection(px,py,pz);
  this->SetFOM(0.0,1,1);
  
  WCSimRecoObjectTable::Instance()->NewVertex();
}


WCSimRecoVertex::WCSimRecoVertex(Double_t x, Double_t y, Double_t z, Double_t t, Double_t px, Double_t py, Double_t pz, Double_t fom, Int_t nsteps, Bool_t pass, Int_t status )
{
  this->Reset();

  this->SetVertex(x,y,z);
  this->SetDirection(px,py,pz);
  this->SetFOM(fom,nsteps,pass);
  this->SetStatus(status);

  WCSimRecoObjectTable::Instance()->NewVertex();
}

WCSimRecoVertex::WCSimRecoVertex(Double_t x, Double_t y, Double_t z, Double_t t, Double_t px, Double_t py, Double_t pz, Double_t angle, Double_t length, Double_t fom, Int_t nsteps, Bool_t pass, Int_t status )
{
  this->Reset();

  this->SetVertex(x,y,z,t);
  this->SetDirection(px,py,pz);
  this->SetConeAngle(angle);
  this->SetTrackLength(length);
  this->SetFOM(fom,nsteps,pass);
  this->SetStatus(status);

  WCSimRecoObjectTable::Instance()->NewVertex();
}

WCSimRecoVertex::~WCSimRecoVertex()
{
  WCSimRecoObjectTable::Instance()->DeleteVertex();
}

void WCSimRecoVertex::SetVertex( Double_t x, Double_t y, Double_t z )
{
  this->SetVertex(x,y,z,950.0);
}

void WCSimRecoVertex::SetVertex( Double_t x, Double_t y, Double_t z, Double_t t )
{
  fX = x;
  fY = y;
  fZ = z;
  fTime  = t;
  fFoundVertex = 1;
}

void WCSimRecoVertex::SetDirection( Double_t px, Double_t py, Double_t pz )
{
  fDirX = px;
  fDirY = py;
  fDirZ = pz;
  fFoundDirection = 1;
}

void WCSimRecoVertex::SetConeAngle( Double_t angle )
{
  fConeAngle = angle;
}

void WCSimRecoVertex::SetTrackLength( Double_t length )
{
  fTrackLength = length;
}

void WCSimRecoVertex::SetFOM( Double_t fom, Int_t nsteps, Bool_t pass )
{
  fFOM = fom;
  fIterations = nsteps;
  fPass = pass;
}

void WCSimRecoVertex::SetStatus( Int_t status )
{
  fStatus = status;
}

void WCSimRecoVertex::Reset()
{ 
  fX = 0.0;
  fY = 0.0;
  fZ = 0.0;
  fTime = 950.0;
  fFoundVertex = 0;

  fDirX = 0.0;
  fDirY = 0.0;
  fDirZ = 0.0;
  fFoundDirection = 0;

  fConeAngle = WCSimParameters::CherenkovAngle();
  fTrackLength = 0.0;

  fFOM = 0.0;
  fIterations = 0;
  fPass = 0;
  
  fStatus = WCSimRecoVertex::kOK;
}

