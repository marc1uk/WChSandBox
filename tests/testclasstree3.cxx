#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TNtuple.h"
#include <TObject.h>

Double_t mrdzstart=1990;
Double_t mrdzlen=919;
Int_t mrdnumlayers=12;		 
Int_t mrdpaddlesperpanel=30;	 
Double_t scintfullxlen = 20;
Double_t scintfullzlen= 0.6;
Double_t scintvfullylen = 155;
Double_t scinthfullylen= 138;
Double_t scintzedges[24] = {19, 25, 75, 81, 150, 156, 225, 231, 300, 306, 375, 381, 450, 456, 525, 531, 600, 606, 675, 681, 750, 756, 825, 831};


class cMRDStrike : public TObject {
// not quite 'hits' because one strike may comprise multiple 'hits', but multiple hits cannot be distinguished from PMT response.
	private: 		// even though these are private they can be histogrammed directly! Even from TBrowser!!
	//Double_t xEstimate;			// this information is stored in the points vector in MRDTrackClass
	//Double_t yEstimate;
	//Double_t zEstimate;
	//Double_t tEstimate;
	Int_t MRDlayer;
	Int_t MRDpaddle;
	Int_t PMTnum;				// assume photons are confined to a paddle.
	std::pair<Double_t, Double_t> xrange;	// from width of panel and whether multiple PMTs were hit.
	std::pair<Double_t, Double_t> yrange;
	std::pair<Double_t, Double_t> zrange;	// from depth of panel
	std::pair<Double_t, Double_t> trange; 	// from uncertainty in PMT timing resoluton
	//Double_t eDeposited;			// can only use a scaled count of PMThits 
	std::vector<Double_t> PMThitTimes;	// vector of hit times.
	
	public:
	cMRDStrike() : MRDlayer(-1), MRDpaddle(-1), PMTnum(-1){}	// default constructor.
	~cMRDStrike(){}	// nothing allocated with new, nothing to delete.
	
	cMRDStrike(std::vector<Double_t> PMThitTimesin, Int_t PMTnumin){	// actual constructor. 
		PMThitTimes = PMThitTimesin;
		PMTnum = PMTnumin;
		
	//std::vector<Int_t> hitpmts = GetNumPMTsHit();	
	//for(std::vector<Int_t>::iterator it=hitpmts.begin(),it!=hitpmts.end(),it++){
		MRDlayer=-1;
		MRDpaddle=-1;
		if(PMThitTimesin.size()==0){MRDlayer=-1; MRDpaddle=-1; PMTnum=-1;}
		else {			
		MRDlayer = TMath::Floor(PMTnum/mrdpaddlesperpanel);
		zrange.first = scintzedges[(Int_t)TMath::Floor(MRDlayer/2)];
		zrange.second = scintzedges[(Int_t)TMath::Floor(MRDlayer/2)+1];
		MRDpaddle = PMTnum-(mrdpaddlesperpanel*MRDlayer);
		if(MRDlayer%2==0){	//horizontal panels, use PMT number to determine y position
			if(MRDpaddle<mrdpaddlesperpanel/2){	// row 1 or 2?
				xrange.first = 0;
				xrange.second = scintvfullylen;
				yrange.first = (scintfullxlen*(MRDpaddle))-(mrdpaddlesperpanel*scintfullxlen/4);
				yrange.second = (scintfullxlen*(MRDpaddle+1))-(mrdpaddlesperpanel*scintfullxlen/4);

			} else {
				xrange.first = -scintvfullylen;
				xrange.second = 0;
				yrange.first = (scintfullxlen*(MRDpaddle-(mrdpaddlesperpanel/2)))-(mrdpaddlesperpanel*scintfullxlen/4);
				yrange.second = (scintfullxlen*(MRDpaddle-(mrdpaddlesperpanel/2)+1))-(mrdpaddlesperpanel*scintfullxlen/4);
			}			
		} else {
			if(MRDpaddle<mrdpaddlesperpanel/2){	// row 1 or 2?
				yrange.first = 0;
				yrange.second = scinthfullylen;
				xrange.first = (scintfullxlen*(MRDpaddle))-(mrdpaddlesperpanel*scintfullxlen/4);
				xrange.second = (scintfullxlen*(MRDpaddle+1))-(mrdpaddlesperpanel*scintfullxlen/4);

			} else {
				yrange.first = -scinthfullylen;
				yrange.second = 0;
				xrange.first = (scintfullxlen*(MRDpaddle-(mrdpaddlesperpanel/2)))-(mrdpaddlesperpanel*scintfullxlen/4);
				xrange.second = (scintfullxlen*(MRDpaddle-(mrdpaddlesperpanel/2)+1))-(mrdpaddlesperpanel*scintfullxlen/4);
			}
		}
		trange.first = PMThitTimes.front();
		trange.second = PMThitTimes.back();
		}
	}

	Int_t GetLayerNum(){return MRDlayer;}
	Int_t GetPaddleNum(){return MRDpaddle;}	// number of paddle within this panel
	Int_t GetPMTNumber(){return PMTnum;}
	Int_t GetNumPMTHits(){return PMThitTimes.size();}
	Double_t GetEdeposited(){ return PMThitTimes.size()*3;}	// units eV
	std::pair<Double_t, Double_t> GetXrange(){return xrange;}
	std::pair<Double_t, Double_t> GetYrange(){return yrange;}
	std::pair<Double_t, Double_t> GetZrange(){return zrange;}
	std::pair<Double_t, Double_t> GetTrange(){return trange;}
	std::vector<Double_t> GetPMThitTimes(){return PMThitTimes;}
	
	ClassDef(cMRDStrike,1);	// INCREMENT VERSION NUM EVERY TIME CLASS MEMBERS CHANGE
};

#ifdef __CINT__
#pragma link C++ class cMRDStrike+;
#endif

void addmrdplots(){
// using custom classes in ROOT test
//**********************************

// Create a ROOT file
TFile* rootfileout = new TFile("TestFile.root", "RECREATE");

TTree* testtree = new TTree("testtree","Tree for testing custom class use");
cMRDStrike* aStrike=0;
testtree->Branch("cMRDStrike",&aStrike);

TRandom R; 


  
for(Int_t e=0;e<10;e++){
	std::vector<Double_t> myhittimes;
	Int_t numhits = TMath::Floor(R.Gaus(10,10));
	cout<<numhits<<endl;
	for(Int_t i=0;i<numhits;i++){myhittimes.push_back(R.Gaus(0,10));}
	Int_t pmtnum = TMath::Floor(R.Gaus(0,300));
	for(Int_t j=0;j<numhits;j++){cout<<"hit time at "<<myhittimes.at(j)<<endl;}
	aStrike = new cMRDStrike(myhittimes, pmtnum);

	testtree->Fill();
}
testtree->Write();
rootfileout->Close();
delete rootfileout;

}


// try to write then read branch of vector of ints
// then vector of class objects
// then class containing just vector of ints
// then class containing vector of class objects


