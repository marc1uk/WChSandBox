#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCut.h"
#include "TFormula.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TPie.h"
#include "TPieSlice.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include <new> // for operator new
#include <stdlib.h>     /* abs */

Double_t mrdzstart=1990;
Double_t mrdzlen=919;
Int_t mrdnumlayers=12;		// 12 layers of scintillator, steel. 13 alu structs, but they're not really important. 
Int_t mrdpaddlesperpanel=30;	// remember TWO ROWS of 15 paddles!!
Double_t maxtrackduration=30;	// used to split tracks in MRD 

//------------------------------------------------------------------
//------------------------------------------------------------------


void addmrdplots(const Char_t* fullinfile="FullEvent.root"){
// will add standard set of graphs to MRDEvents.root file
//*******************************************************
/*// Open mrd hit file and check valid
TFile mrdfile = TFile(mrdinfile);
if (mrdfile.IsZombie()) {
	std::cout << "Error opening MRD file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* mrdtree = (TTree*)mrdfile.Get("MRDTree"); // not strictly necessary - can access MRDTree directly, at least in CINT.
*/
// Open full event file and check valid
TFile fullfile = TFile(fullinfile);
if (fullfile.IsZombie()) {
	std::cout << "Error opening FullEvent file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* fulltree = (TTree*)fullfile.Get("EventTree");

// PMT hit check: Let's check the maximum x and y positions of photons in the MRD.
//Steps:
//1) get number of hits this entry
//2) create dynamic array of appropriate size for all z positions of hits
//3) get z position branch into that array
//4) scan through the array and increment zmax found so far.

Int_t numentries = fulltree->GetEntries();
cout<<"numentries:"<<numentries<<endl;
TBranch* numhitsbranch = fulltree->GetBranch("nphot");
Int_t arraysize = 0;
numhitsbranch->SetAddress(&arraysize);
TBranch* xEndBranch = fulltree->GetBranch("phot_xEnd");
TBranch* yEndBranch = fulltree->GetBranch("phot_yEnd");
TBranch* zEndBranch = fulltree->GetBranch("phot_zEnd");
Double_t* xEndVals;
Double_t xEndVal = 0;
Double_t xEndValMax=0;
Double_t* yEndVals;
Double_t yEndVal = 0;
Double_t yEndValMax=0;
Double_t* zEndVals;
Double_t zEndVal = 0;
for (Int_t i=0;i<numentries;i++){
	numhitsbranch->GetEntry(i);
	cout<<arraysize<<endl;
	if(arraysize!=0){ 
		xEndVals = new Double_t[arraysize];
		xEndBranch->SetAddress(xEndVals);
		xEndBranch->GetEntry(i);
		yEndVals = new Double_t[arraysize];
		yEndBranch->SetAddress(yEndVals);
		yEndBranch->GetEntry(i);
		zEndVals = new Double_t[arraysize];
		zEndBranch->SetAddress(zEndVals);
		zEndBranch->GetEntry(i);
		for(size_t j=0;j<arraysize;j++){
			zEndVal=zEndVals[j];
			xEndVal=xEndVals[j];
			if((abs(xEndVal)>abs(xEndValMax))&&zEndVal>2000){xEndValMax=xEndVal;}
			yEndVal=yEndVals[j];
			if((abs(yEndVal)>abs(yEndValMax))&&zEndVal>2000){yEndValMax=yEndVal;}
		}
		delete[] xEndVals;
		xEndVals=0;
		delete[] yEndVals;
		yEndVals=0;
	}
}
cout<<"maximum x of a photon end: "<<xEndValMax<<"maximum y: "<<yEndValMax<<endl;
// 276cm end to end horizontally (300cm full width of 15x20cm vertical strips)
// 310cm end to end vertically. 
fulltree->ResetBranchAddresses();
fullfile.Close();
cout<<"marco"<<endl;
// for 6 events; max xmin = -93cm, ymin = -101cm. nowhere near. 
}

