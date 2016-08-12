#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "strikeclasstest.hh"

//#ifdef __CINT__
//#pragma link C++ class cMRDStrike+;
//#endif

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






