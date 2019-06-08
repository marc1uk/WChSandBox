//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4TrajectoryPoint.cc,v 1.19 2006/06/29 21:16:15 gunter Exp $
// GEANT4 tag $Name: geant4-09-01-patch-01 $
//
// ---------------------------------------------------------------
//
// G4TrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4TrajectoryPoint.hh"
#include "WCLiteTrajectoryPoint.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UnitsTable.hh"

//#define G4ATTDEBUG
#ifdef G4ATTDEBUG
#include "G4AttCheck.hh"
#endif

G4Allocator<WCLiteTrajectoryPoint> WCLiteTrajectoryPointAllocator;

WCLiteTrajectoryPoint::WCLiteTrajectoryPoint()
{
  fPosition = G4ThreeVector(0.,0.,0.);
  fDir = G4ThreeVector(0.,0.,0.);
  _globaltime=0;
  _TotalEnergy = -666;
  _KinEnergy = -666;
  _ProcessName="Empty";
  _DetectorLocation="Empty";

 G4cout<<"NOOOOOOOOOOOOOOOOOOOOOOOOO!123"<<G4endl;
}

WCLiteTrajectoryPoint::WCLiteTrajectoryPoint(G4ThreeVector pos, G4ThreeVector dir, G4double theGlobaltime, G4double thetotE, G4double theKE, G4String theproc, G4String theDetLoc)
{
  //  G4cout<<"CONSTRUCTING TRAJPOINT"<<G4endl;

  fPosition = pos;
  fDir = dir;
  _globaltime = theGlobaltime;
  _TotalEnergy = thetotE;
  _KinEnergy = theKE;
  _ProcessName=theproc;
  _DetectorLocation=theDetLoc;

  //  G4cout<<"&*&*&*&*&*&*&* Traj Point: "<<pos.x()<<" "<<pos.y()<<" "<<pos.z()<<" "<<_KinEnergy<<" "<<_DetectorLocation<<G4endl;
}

WCLiteTrajectoryPoint::WCLiteTrajectoryPoint(const WCLiteTrajectoryPoint &right)
  : G4TrajectoryPoint(),fPosition(right.fPosition),fDir(right.fDir),_globaltime(right._globaltime),_TotalEnergy(right._TotalEnergy),_KinEnergy(right._KinEnergy),_ProcessName(right._ProcessName),_DetectorLocation(right._DetectorLocation)
{
 G4cout<<"NOOOOOOOOOOOOOOOOOOOOOOOOO!"<<G4endl;
}

WCLiteTrajectoryPoint::~WCLiteTrajectoryPoint()
{
}

G4double WCLiteTrajectoryPoint::GetGlobalTime()
{
  return _globaltime;
}

G4double WCLiteTrajectoryPoint::GetTotalEnergy()
{
  return _TotalEnergy;
}

G4double WCLiteTrajectoryPoint::GetKineticEnergy()
{
  return _KinEnergy;
}

G4String WCLiteTrajectoryPoint::GetProcessName()
{
  return _ProcessName;
}

G4String WCLiteTrajectoryPoint::GetDetectorLocation()
{
  return _DetectorLocation;
}
