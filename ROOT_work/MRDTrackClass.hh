#ifndef _MRDTrack_Class_
#define _MRDTrack_Class_ 1

#include <TObject.h>
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TVector3.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include <exception>	// for stdexcept
#include <vector>

//#include "MRDStrikeClass.hh"
#include "MRDspecs.hh"
extern Double_t mrdzstart;	// note; "extern Double_t mrdzstart, mrdzlen;" will NOT work. Must be separate lines. 
extern Double_t mrdzlen;
extern Int_t mrdnumlayers;
extern Int_t mrdpaddlesperpanel;
extern Double_t scintzedges[24];

class cMRDTrack : public TObject {
	private:	
	Int_t particleID;	// particle type code
//	Double_t xStart;	// don't need to store these separately
//	Double_t yStart;	
//	Doublt_t zStart;
//	Double_t xEnd;		
//	Double_t yEnd;
//	Double_t zEnd;
	Double_t KEStart;
	Double_t KEEnd;
	Int_t valid;	// track validity; has it been checked that all hit panels are consistent
	// can we/do we need to do p as well as m?
	std::vector<ROOT::Math::XYZTVector> points;	// estimated times and positions of interaction points, or 'strikes'.
	// Due to layered structure, estimate of each strike's position is calculated using info from multiple strikes.
	// Having best-estimate positions in LorentzVector format may be helpful for reco plotting and analysis.
	std::vector<cMRDStrike> strikes;			//vector of cMRDStrike class objects. 
	// strike class objects contain all information derivable from just one strike object - layer number, number and vector of 
	// PMT hits, number of separate PMTs hit, energy deposition, time, positional error from struck PMTs... 
	
	public:		
	cMRDTrack() : particleID(0), KEStart(0), KEEnd(0), valid(0){}
	~cMRDTrack(){}	// vector container automatically handles resource deallocation: destruction of points and strikes
	cMRDTrack(std::vector<cMRDStrike> strikesin) : particleID(0), KEStart(0), KEEnd(0), valid(0), strikes(strikesin){DoReconstruction();}
	
	void SetParticleID(Int_t particleIDin){particleID = particleIDin;}
	Int_t GetParticleID(){return particleID;}

//	void SetxStart(Double_t xStartin){xStart = xStartin;}
	Double_t GetxStart(){return points.front().X();}
//	void SetyStart(Double_t yStartin){yStart = yStartin;}
	Double_t GetyStart(){return points.front().Y();}
//	void SetzStart(Double_t zStartin){zStart = zStartin;}
	Double_t GetzStart(){return points.front().Z();}

//	void SetxEnd(Double_t xEndin){xEnd = xEndin;}
	Double_t GetxEnd(){return points.back().X();}
//	void SetyEnd(Double_t yEndin){yEnd = yEndin;}
	Double_t GetyEnd(){return points.back().Y();}
//	void SetzEnd(Double_t zEndin){zEnd = zEndin;}
	Double_t GetzEnd(){return points.back().Z();}	

	void SetKEStart(Double_t KEStartin){KEStart = KEStartin;}
	Double_t GetKEStart(){return KEStart;}
	void SetKEEnd(Double_t KEEndin){KEEnd = KEEndin;}
	Double_t GetKEEnd(){return KEEnd;}

	Int_t GetNumStrikes(){return strikes.size();}
	Double_t GetTotalEdep(){return KEStart-KEEnd;}
	Double_t GetPenetration(){return points.back().Z()-points.front().Z();}
	
	std::vector<ROOT::Math::XYZTVector> GetPoints(){return points;}
	ROOT::Math::XYZTVector* GetFirstPoint(){return &(points.front());}		// should check we have at least one.
	ROOT::Math::XYZTVector* GetLastPoint(){return &(points.back());}		// should check we have at least one.
	ROOT::Math::XYZTVector* GetPoint(Int_t i){
		try{return &(points.at(i));}
		catch(const std::out_of_range& oor){return 0;}
	}	
	
	class RecoRegion : public TObject {
	
	//we also need to store EDGES around a known pass-through region so we can project back at known angles from
	//a given point. centre vertex is not really suitable - we need x, y planes at a particular z
		private:
		TVector3 hvertex;
		TVector3 vvertex;
		Double_t hanglemin;
		Double_t hanglemax;
		Double_t vanglemin;
		Double_t vanglemax;
		public:
		TVector3* GimmeHvertex(){return hvertex;}
		TVector3* GimmeVvertex(){return vvertex;}
		void SetHangleMin(Double_t anglein){hanglemin=anglein;}
		void SetHangleMax(Double_t anglein){hanglemax=anglein;}
		void SetVangleMin(Double_t anglein){vanglemin=anglein;}
		void SetVangleMax(Double_t anglein){vanglemax=anglein;}
		Bool_t AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin);
		Bool_t CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin)
		const std::vector<Double_t> GetAngles(){
			std::vector<Double_t> holder;
			holder.push_back(hanglemin);
			holder.push_back(hanglemax);
			holder.push_back(vanglemin);
			holder.push_back(vanglemax);
			return holder;
		}
	}
	
	std::vector<cMRDStrike> GetStrikes(){return strikes;}	
	cMRDStrike* GetFirstStrike(){return &(strikes.front());}	// should check we have at least one.
	cMRDStrike* GetLastStrike(){return &(strikes.back());}		// should check we have at least one.
	cMRDStrike* GetStrike(Int_t i){
		try{return &(strikes.at(i));}
		catch(const std::out_of_range& oor){return 0;}
	}	
	void AppendStrike(cMRDStrike strikein){strikes.push_back(strikein);}

	ClassDef(cMRDTrack,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

// ********* trajectory reconstruction: main *****************
	void DoRegionReconstruction(){
		std::vector<TVector3> bottompoints;
		std::vector<TVector3> toppoints;
		for(Int_t i=0;i<strikes.size();i++){
			cMRDStrike thisstrike=strikes.at(i);
			bottompoints.push_back(TVector3(thisstrike.GetXrange().first,thisstrike.GetYrange().first,thisstrike.GetZrange().first));
			toppoints.push_back(TVector3(thisstrike.GetXrange().second,thisstrike.GetYrange().second,thisstrike.GetZrange().second));
		}
		Double_t angmax = (0.0/0.0);	// set to NaN angmax represents angle of uppermost limit of valid trajectories
		Double_t angmin = (0.0/0.0);	// set to NaN. angmin represents angle of lowermost limit of valid trajectories 
		Int_t layermax = -1;
		Int_t layermin = -1;
		
		// scan through from all bottom points to all top points to find the greatest upward angle a particle may
		// have had and still hit all recorded layers (note: this may be negative, indicating a particle was
		// going downward - in this case angmin is the angle at which it went downward _least_ steeply)
		// angle is measured from the positive z axis
		for(Int_t startlayer=0;startlayer<strikes.size();startlayer++){
			for(Int_t nextlayer=startlayer+1;nextlayer<strikes.size(); nextlayer++){
				Double_t opp = toppoints.at(nextlayer).X()-bottompoints.at(startlayer).X();
				Double_t adj = toppoints.at(nextlayer).Z()-bottompoints.at(startlayer).Z();
				Double_t ang = TMath::Atan(opp/adj);
				if(TMath::IsNaN(angmin)||(ang<angmin)){
					if(AngValid(startlayer,ang,1)){
						angmin=ang;	
						layermin=startlayer;
					}					
				}
			}
			for(Int_t nextlayer=hitlayer+1;nextlayer<strikes.size(); nextlayer++){
				Double_t opp = bottompoints.at(nextlayer).X()-toppoints.at(hitlayer).X();	// negative if downgoing
				Double_t adj = toppoints.at(nextlayer).Z()-bottompoints.at(hitlayer).Z();
				Double_t ang = TMath::Atan(opp/adj);										// negative if downgoing
				if(TMath::IsNaN(angmax)||(ang>angmax)){
					if(AngValid(nextlayer,ang,0)){     // project forward both forward and back from bottom vertex (0) of
						angmax=ang;						// nextlayerf to all other layers and check within limits
						layermax=nextlayer;
					}										
				}

		}
		if(TMath::IsNaN(angmax)||Tmath::IsNaN(angmin)){
			valid=1;
			if(layermax>layermin){
				// use cosine rule from layermin
			} else {
				// use cosine rule from layermax
			} else {valid=0;}
		}
			
}

// ********* trajectory reconstruction: scan projections to all layers *****************
Bool_t cMRDStrike::AngValid(Int_t layerstart, Double_t angle, Int_t MaxMin){
	Bool_t angvalid=true;
	for(Int_t layerend=0;layerend<strikes.size();layerend++){
		if(layerstart==layerend){continue;}	// don't check from a layer to itself
		angvalid = CheckIntersection(layerstart, layerend, angle, MaxMin);
		if(!angvalid){return false;}
	}
}

// ********* trajectory reconstruction: check single projection *****************
Bool_t cMRDStrike::CheckIntersection(Int_t layerstart, Int_t layerend, Double_t angle, Int_t MaxMin){
	Double_t z1 = strikes.at(layerstart).GetZrange().first;
	Double_t z2 = strikes.at(layerend).GetZrange().first;
	Double_t zstart;
	Double_t zend;
	cMRDStrike LHlayer;
	cMRDStrike RHlayer;
	if(z1<z2){
		LHlayer = strikes.at(layerstart);
		RHlayer = strikes.at(layerend);
		zstart=z1; 
		zend=z2;
	} else {
		LHlayer = strikes.at(layerend);
		RHlayer = strikes.at(layerstart);
		zstart=z2; 
		zend=z1;
	}
	Double_t ydiff = (zend-zstart)*TMath::Tan(angle);
	if(MaxMin==1){	
	// layerstart defines a maximum angle - check projection from it's top to layerend at 'angle'
	// downwards (or at 'angle' upwards if projecting backward) strikes layerend
		if(z1==zstart){	// layerstart is on LHS - project forwards (-'ve y displacement)
			yproj = strikes.at(layerstart).GetYrange().second - ydiff;
		} else {		// layerstart is on RHS - project backwards (+'ve y displacement)
			yproj = strikes.at(layerstart).GetYrange().second + ydiff;
		}
	} else {
	// layerstart defines a min angle - check projection from it's bottom at angle upwards 
	// (if projecting forwards, or downwards if backwards) to layerend, strikes within 
	// layerend's y boundaries
		if(z1==zstart){	// layerstart is on LHS - project forwards (+'ve y displacement)
			yproj = strikes.at(layerstart).GetYrange().first + ydiff;
		} else {		// layerstart is on RHS - project backwards (-'ve y displacement)
			yproj = strikes.at(layerstart).GetYrange().first - ydiff;
		}
	}
	if( (yproj >= strikes.at(layerend).GetYrange().first)&&
		(yproj <= strikes.at(layerend).GetYrange().second) ){ return true; }
	else { return false; }
}

#endif

#ifdef __CINT__
#pragma link C++ class cMRDTrack+;
//#pragma link C++ class ROOT::Math::XYZTVector+;
//#pragma link C++ class std::vector<ROOT::Math::XYZTVector>+;
//#pragma link C++ class cMRDStrike+;
#pragma link C++ class std::vector<cMRDStrike>+;
#endif
