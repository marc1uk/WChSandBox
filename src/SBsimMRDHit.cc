// ====================================================================
//   SBsimMRDHit.cc
//
//   2006/03/03 K. Hiraide
// ====================================================================
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "SBsimMRDHit.hh"

#include <iostream>
#include "G4SystemOfUnits.hh"

// allocator
G4Allocator<SBsimMRDHit> SBsimMRDHitAllocator;

//////////////////////////
SBsimMRDHit::SBsimMRDHit()
//////////////////////////
{
}

///////////////////////////
SBsimMRDHit::~SBsimMRDHit()
///////////////////////////
{
}


////////////////////////
void SBsimMRDHit::Draw()
////////////////////////
{
  /*
  G4VVisManager* pVVisManager= G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(hitpos);
    circle.SetScreenSize(5.0);
    circle.SetFillStyle(G4Circle::filled);

    G4Color color, goodColor(1.,0.,0.), badColor(0.,0.,1.);

    color=goodColor;

    G4VisAttributes attribs(color);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
  */
  /*
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Transform3D trans(rot.inverse(),hitpos);
      G4VisAttributes attribs;
      const G4VisAttributes* pVA = pLogV->GetVisAttributes();
      if(pVA) attribs = *pVA;
      G4Colour colour(1.,0.,0.);
      attribs.SetColour(colour);
      attribs.SetForceWireframe(false);
      attribs.SetForceSolid(true);
      pVVisManager->Draw(*pLogV,attribs,trans);
    }
  */
}

/////////////////////////
void SBsimMRDHit::Print()
/////////////////////////
{
  //  G4cout.setf( ios::fixed );
  G4cout.precision(2);
  G4cout << "id= "    << id << ", "
	 << "track= " << itrack   << "("
	 << parname   << "), "
         << "t= "    << hittime/ns     << " ns, ";
  G4cout.precision(0);
  G4cout << "pos = "  << hitpos*(1./cm) << " cm,"
	 << G4endl;
}
