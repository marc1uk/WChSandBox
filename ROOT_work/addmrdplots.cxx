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
#include <exception>	// for stdexcept
#include "TClonesArray.h"

#include "MRDStrikeClass.hh"
//#include "MRDTrackClass.hh"

/*
&&&&&&&&&&&&&&&&&&&&&&&&&&
NOTE: ALL LOCAL DECLARATIONS TCANV THIS = TCANV("NAME","TITLE"), TFILE F = TFILE("FILENAME.ROOT","READ") ETC
WILL NEED TO BE REPLACED BY POINTER VERSIONS???
&&&&&&&&&&&&&&&&&&&&&&&&&&

#include "MRDspecs.hh"
extern Double_t mrdzstart;	// note; "extern Double_t mrdzstart, mrdzlen;" will NOT work. Must be separate lines. 
extern Double_t mrdzlen;
extern Int_t mrdnumlayers;
extern Int_t mrdpaddlesperpanel;
extern Double_t scintzedges[24];

std::map<TString,Int_t> ParticleNamesToCodes;
Int_t ConvertParticleNameToCode(TString particleName);

std::map<TString,Int_t> ProcessNamesToCodes;
Int_t ConvertProcessNameToCode(TString processName);	

std::map<Int_t,TString> ParticleCodesToNames;
TString ConvertParticleCodeToName(Int_t particleCode);

std::map<Int_t,TString> ProcessCodesToNames;
TString ConvertProcessCodeToName(Int_t processCode);


typedef std::vector< std::pair<TString,Int_t> > stringintpairvec;
void GenerateMaps();
void GenerateVectors(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void GeneratePieCharts(stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void MakeEnergyDepPlots(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void PlotPenetration(TTree* mrdtree);
void EnergyDepositioninLayers(TTree* mrdtree);
void EnergyDepositionInPaddles(TTree* mrdtree);
void countmrdhits(TTree* fulltree);
void splitindividualsources(TTree* mrdtree, TTree* recotree);
void splitPMTtracks(TTree* pmttree, TTree* recotree);
Int_t findmode(Int_t* arrayin, Int_t arraysize);
void checkpmtclass(TTree* recotree);

Double_t maxtrackduration=30;	  // used to split tracks in MRD 
Double_t maxpmttrackduration=300; // used to split tracks in PMT hits

TFile* rootfileout;

//------------------------------------------------------------------
//------------------------------------------------------------------


void addmrdplots(const Char_t* mrdinfile="MRDEvents2.root",const Char_t* fullinfile="FullEvent2.root"){
// will add standard set of graphs to MRDEvents.root file
//*******************************************************
//gSystem->Load("libPhysics");
//gSystem->Load("libGenVector");

// Open mrd hit file and check valid
TFile mrdfile = TFile(mrdinfile, "READ");
if (mrdfile.IsZombie()) {
	std::cout << "Error opening MRD file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* mrdtree = (TTree*)mrdfile.Get("MRDTree"); // not strictly necessary - can access MRDTree directly, at least in CINT.
// Open PMT events tree
TTree* pmttree = (TTree*)mrdfile.Get("PMTTree"); // not strictly necessary - can access MRDTree directly, at least in CINT.

// Open full event file and check valid
TFile fullfile = TFile(fullinfile, "READ");
if (fullfile.IsZombie()) {
	std::cout << "Error opening FullEvent file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* fulltree = (TTree*)fullfile.Get("EventTree");

// Open output file
rootfileout = new TFile("MRDEventsOut.root", "RECREATE");

GenerateMaps(); // fills maps to convert particle and process names to codes and vice versa.

/*stringintpairvec mrdparticleids;
stringintpairvec mrdprocessids;

GenerateVectors(mrdtree,mrdparticleids,mrdprocessids);

GeneratePieCharts(mrdparticleids, mrdprocessids);

MakeEnergyDepPlots(mrdtree, mrdparticleids, mrdprocessids);

PlotPenetration(mrdtree);


EnergyDepositionInPaddles(mrdtree);
EnergyDepositioninLayers(mrdtree);
*/
TTree* recotree = new TTree("recotree","Tree for reconstruction data");
recotree->SetEntries(mrdtree->GetEntries());	
// number of entries (events) must be set manually if using TBranch->Fill() rather than TTree->Fill(). 
// we do this so we can fill MRD hit reco info and PMT hit reco info separately.

splitindividualsources(mrdtree, recotree);

splitPMTtracks(pmttree,recotree);

checkpmtclass(recotree);

//countmrdhits(fulltree);

mrdfile.Close();
fullfile.Close();
//rootfileout->Write(); //Write all objects to disk. 
recotree->AddFriend("MRDTree",mrdinfile);	// only do this after the mrdtree file is closed.
//recotree->Write();
rootfileout->Close();
delete rootfileout;


//delete mrdfile;	//close the file (only for pointers)
cout<<"marco"<<endl;
}


void GenerateVectors(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids){
TCanvas vectCanv = TCanvas("vectCanv","Title");
vectCanv.cd();	// ensure selected

// Find the particle IDs of particles that have interacted (probably e's, mu's and gammas):

// histogram needs a LARGE number of bins - large range of particle codes, but must have resolution <1.0
mrdtree->Draw("mrdhit_particleID>>histpid(3000)");	
TH1F* histpid = (TH1F*)vectCanv.GetPrimitive("histpid");
Int_t numbins = histpid->GetNbinsX();
// scan through the histogram for non-empty bins - these are the particle IDs that appear in our tree
Int_t binconts=0;
//std::vector< std::pair<TString,Int_t> > mrdparticleids;
std::pair <TString,Int_t> typeandnumhits;
for(Int_t i=0;i<numbins;i++){
	binconts = histpid->GetBinContent(i);
	if(binconts!=0){
		TString pname = ConvertParticleCodeToName(histpid->GetXaxis()->GetBinLowEdge(i));
		typeandnumhits = std::make_pair(pname,binconts);
		mrdparticleids.push_back(typeandnumhits);
		cout<<binconts<<" hits of type "<<pname<<endl;
	}
} 

// Find the process IDs of processes that have interacted (probably scint, phot, brem):
mrdtree->Draw("mrdhit_process>>histproc");
TH1F* histproc = (TH1F*)vectCanv.GetPrimitive("histproc");
numbins = histproc->GetNbinsX();
// scan through the histogram for non-empty bins - these are the particle IDs that appear in our tree
binconts=0;
//std::vector< std::pair<TString,Int_t> > mrdprocessids;
for(Int_t i=0;i<numbins;i++){
	binconts = histproc->GetBinContent(i);
	if(binconts!=0){
		TString pname = ConvertProcessCodeToName(histproc->GetXaxis()->GetBinLowEdge(i));
		typeandnumhits = std::make_pair(pname,binconts);
		mrdprocessids.push_back(typeandnumhits);
		cout<<binconts<<" hits of type "<<pname<<endl;
	}
}

}

void GeneratePieCharts(stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids){
// this function generates pie charts of the particles and processes that have interacted/occurred in the MRD. 

// Produce pie chart of particles that interacted in the MRD, and their frequency:
// ============================================================================
TCanvas pieCanv = TCanvas("pieCanv","Pie Charts",385,110,700,867);
pieCanv.Divide(1,2);
TPie mrdparticlepie = TPie("mrdparticlepie", "Particle hits on the MRD",mrdparticleids.size());
for(size_t i=0;i<mrdparticleids.size();i++){
	//TPieSlice pieslice = mrdparticlepie.GetSlice(i);
	//pieslice->SetValue(mrdparticleids.at(i).second);
	//pieslice->SetName((mrdparticleids.at(i)).first); // doesn't seem to work.
	TString pname = (mrdparticleids.at(i)).first;
	Double_t pcount = (mrdparticleids.at(i)).second;
	mrdparticlepie.SetEntryVal(i,pcount);
	mrdparticlepie.SetEntryLabel(i,pname.Data());
}
mrdparticlepie.SetAngularOffset(333); 
mrdparticlepie.SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
mrdparticlepie.SetValueFormat("%4.0f");
mrdparticlepie.SetPercentFormat("%3.0f");
mrdparticlepie.SetCircle(0.5, 0.4702026, 0.3302274);
mrdparticlepie.SetTextSize(0.03455766);
//mrdparticlepie.SetLabelsOffset(0.05);	// shift labels from pie centre so they don't mash together.
mrdparticlepie.Write();
pieCanv.cd(1);
mrdparticlepie.Draw("3d");
//** local print for check
/* // divide works fine but this lets set position. just 'cos its a mess.
TCanvas* canvloc1 = new TCanvas("canvloc1","Pie Charts",385,110,700,867); //for pie charts
//canvloc1->Divide(1,2);
canvloc1->cd();
TPad *canvloc1_1 = new TPad("canvloc1_1", "canvloc1_1",0.01,0.51,0.99,0.99);
canvloc1_1->Draw();
canvloc1->cd();
TPad *canvloc1_2 = new TPad("canvloc1_2", "canvloc1_2",0.01,0.01,0.99,0.49);
canvloc1_2->Draw();

TPie* mrdparticlepielocal = new TPie(mrdparticlepie);
mrdparticlepielocal->SetAngularOffset(333.);
mrdparticlepielocal->SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
mrdparticlepielocal->SetValueFormat("%4.0f");
mrdparticlepielocal->SetPercentFormat("%3.0f");
mrdparticlepielocal->SetCircle(0.5, 0.4702026, 0.3302274);
mrdparticlepielocal->SetTextSize(0.03455766);
//mrdparticlepielocal->SetLabelsOffset(0.05);	// shift labels from pie centre so they don't mash together.
//canvloc1->cd(1);
canvloc1_1->cd();
mrdparticlepielocal->Draw("3d");
*/

// Produce pie chart of processes that occurred in the MRD, and their frequency
// ============================================================================
TPie *mrdprocesspie = new TPie("mrdprocesspie", "Interaction processes in the MRD",mrdprocessids.size());
for(size_t i=0;i<mrdprocessids.size();i++){
	TString pname = (mrdprocessids.at(i)).first;
	Double_t pcount = (mrdprocessids.at(i)).second;
	mrdprocesspie->SetEntryVal(i,pcount);
	mrdprocesspie->SetEntryLabel(i,pname.Data());
}
mrdprocesspie->SetAngularOffset(30.);
mrdprocesspie->SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
mrdprocesspie->SetValueFormat("%4.0f");
mrdprocesspie->SetPercentFormat("%2.0f");
mrdprocesspie->SetCircle(0.5, 0.4702026, 0.3302274);	// reduce size so labels actually fit on the canvas
mrdprocesspie->SetTextSize(0.03455766);
//mrdprocesspie->SetLabelsOffset(0.05);
rootfileout->cd();
mrdprocesspie->Write();
pieCanv.cd(2);
mrdprocesspie->Draw("3d");
pieCanv.SaveAs("piecharts.png");
//** local print for check
/*TPie* mrdprocesspielocal = new TPie(*mrdprocesspie);
mrdprocesspielocal->SetLabelFormat("#splitline{%txt}{#splitline{%val}{(%perc)}}");
mrdprocesspielocal->SetValueFormat("%4.0f");
mrdprocesspielocal->SetPercentFormat("%1.0f");
mrdprocesspielocal->SetCircle(0.5, 0.4702026, 0.3302274);	// reduce size so labels actually fit on the canvas
mrdprocesspielocal->SetTextSize(0.03455766);
//mrdprocesspielocal->SetLabelsOffset(0.05);
//canvloc1->cd(2);
canvloc1_2->cd();
mrdprocesspielocal->Draw("3d");
//canvloc1->SaveAs("piecharts.png");
*/

delete mrdprocesspie;
}

void MakeEnergyDepPlots(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids){
// Produce histogram of energy deposition, with curves for independent processes and particles
// ===========================================================================================
EColor mycolours[11] = {kBlack, kBlue, (EColor)TColor::GetColorDark(kGreen), kRed, kViolet, kOrange, kMagenta,(EColor)(kAzure+2),(EColor)(kOrange+4),(EColor)(kViolet-6),(EColor)(kTeal-6)};

TCanvas stackCanv = TCanvas("stackCanv","Title");
std::vector<TH1F*> allhistoslocal;
mrdtree->Draw("mrdhit_edep>>histedep(100)");
TH1F* histedep = (TH1F*)stackCanv.GetPrimitive("histedep");
stackCanv.cd();
THStack* alledeposits = new THStack("alledeposits","Energy Desposition in MRD");
Char_t partcut[80];
TCut thiscut;
TLegend* alegend = new TLegend(0.6460565,0.5833333,0.8862275,0.8830409,"Sources");
Int_t colournum=2;
TString currentparticle;
Int_t particleidnum;
TString currentprocess;
Int_t processnum;
Char_t legendentry[80];
for (size_t i=0;i<mrdparticleids.size();i++){
	for (size_t j=0;j<mrdprocessids.size();j++){
		histedep->SetLineColor(mycolours[colournum]); //set line colour 
		histedep->SetLineWidth(2);
		currentparticle = (mrdparticleids.at(i)).first;
		//cout<< "currentparticle: "<<currentparticle<<", ";
		particleidnum = ConvertParticleNameToCode(currentparticle);
		sprintf(partcut,"%s %d","mrdhit_particleID ==",particleidnum);
		currentprocess = (mrdprocessids.at(j)).first;
		//cout<< "currentprocess: "<<currentprocess<<endl;
		processnum = ConvertProcessNameToCode(currentprocess);
		sprintf(partcut,"%s %s %d",partcut,"&& mrdhit_process ==",processnum);
		//cout << "partcut: " << partcut << endl;
		thiscut = partcut;
		mrdtree->Draw("mrdhit_edep>>histedep",thiscut);
		//cout<<"cut passes "<<(histedep->GetEntries())<<" entries"<<endl;
		if(histedep->GetEntries()!=0){
			TH1F* newhisto = new TH1F(*histedep);
			if(colournum<((sizeof(mycolours)/sizeof(EColor))-1)){colournum++;}
			stackCanv.cd();
			newhisto->Draw();
			stackCanv.SaveAs(Form("%s-%d%s","edephist",colournum-2,".png"));
			newhisto->SetStats(0);
			alledeposits->Add(newhisto);
			//cout<<"added "<<newhisto<<" to the stack - now has: "<<(alledeposits->GetStack()->GetEntries())<<" histos"<<endl;
			const char* buffer1=currentparticle.Data();
			const char* buffer2=currentprocess.Data();
			sprintf(legendentry,"%s %s",buffer1,buffer2);
			alegend->AddEntry(newhisto,legendentry); //add a line with function and message //histoelocal
			//** local version for check
			/*TH1F* histoelocal = new TH1F(*newhisto);
			histoelocal->SetDirectory(0);
			allhistoslocal.push_back(histoelocal);*/
		}
	}
}

//cout<<"finished adding: the stack has: "<<(alledeposits->GetStack()->GetEntries())<<" histos"<<endl;
TH1F* lastentry = (TH1F*)alledeposits->GetStack()->Last();
TH1F* sumedep = new TH1F(*lastentry);
sumedep->SetStats(0);
sumedep->SetLineColor(mycolours[1]);
sumedep->SetLineWidth(3);
alledeposits->Add(sumedep);
stackCanv.cd();
alledeposits->Draw("nostack");
//stackCanv.BuildLegend(); - this doesn't work.
alegend->AddEntry(sumedep,"Total"); //add a line with function and message
alegend->SetTextSize(0.03);
alegend->SetFillColor(0);
alegend->Draw(); //draw legend.
stackCanv.SaveAs("plotstack.png");
rootfileout->cd();
alledeposits->Write();

//** local print for check - also need to comment out local version code in for loop.
//Note: legend is now generated with file resident histograms - uncommenting to show local may not display the legend as it takes pointers to objects that will go out of scope when the macro ends. If necessary change legend line in loop to use histoelocal instead of newhisto and use the below code for the sum entry. This legend will work for file too, but sum will need to be generated before save -> i.e. move this 'local' code before the paragraph above. 
/*TCanvas* canvloc2 = new TCanvas("canvloc2","Energy Deposition Chart",79,59,1173,712); // for energy deposition vs length
canvloc2->SetLogy();
TH1F* sumhistolocal = (TH1F*)allhistoslocal.at(0)->Clone("sumhistolocal");
sumhistolocal->SetDirectory(0);
// Note when plotting using 'same' only the titles of the first histogram are kept.
allhistoslocal.at(0)->SetTitle("Energy Deposition in the MRD;Energy Deposited (keV); Frequency"); 
allhistoslocal.at(0)->GetXaxis()->CenterTitle(true);
allhistoslocal.at(0)->GetYaxis()->CenterTitle(true);
canvloc2->cd();
allhistoslocal.at(0)->Draw();
for(size_t i=1;i<allhistoslocal.size()-1;i++){
	sumhistolocal->Add(allhistoslocal.at(i));
	allhistoslocal.at(i)->Draw("same");
}
sumhistolocal->SetLineColor(mycolours[1]);	//http://root.cern.ch/root/html/TAttLine.html (kBlue)
sumhistolocal->SetLineWidth(3);
sumhistolocal->Draw("same");
alegend->AddEntry(sumhistolocal,"Total"); //add a line with function and message
alegend->SetTextSize(0.03);
alegend->SetFillColor(0);
alegend->Draw();
*/

delete alledeposits;
delete alegend;
delete sumedep;
}

void PlotPenetration(TTree* mrdtree){
// Produce histogram of penetration depth with bins corresponding to MRD layers
// =============================================================================
TCanvas penetrCanv = TCanvas("penetrCanv","Title");
penetrCanv.cd();
TH1F mrdzlayer = TH1F("mrdzlayer","Layer Penetration in MRD;Layer Number;Count",mrdnumlayers,1,mrdnumlayers+1);
// TODO: next geant4 run, need to update 1990->2000, 2909->2919 to account for change in geant4 code. 
//water tank radius is 1.98m, geant4 sim adds 2cm offset so MRD is at 2000mm. MRD total length is 919mm calculated from 12x5cm steel plates, 12x0.6cm scintillator panels and 13x1.9cm aluminium support structures. Alu frame thicknesses may not be quite right. No spacing is accounted for between the layers. 
mrdtree->SetAlias("mrdhit_zlayer",Form("(((mrdhit_z-%f)/%f)*%d)+1",mrdzstart,mrdzlen,mrdnumlayers));
//cout<<"alias def: "<<(Form("(((mrdhit_z-%f)/%f)*%d)+1",mrdzstart,mrdzlen,mrdnumlayers))<<endl;
penetrCanv.cd();
mrdtree->Draw("mrdhit_zlayer>>mrdzlayer");
mrdzlayer.GetXaxis()->SetNdivisions(mrdnumlayers+1);
mrdzlayer.GetXaxis()->CenterLabels();	//undocumented!
mrdzlayer.SetStats(0);
penetrCanv.SaveAs("zlayer.png");
mrdzlayer.Write();
//** local print for check
/*TCanvas* canvloc3 = new TCanvas("canvloc3","Frequency of Layer Penetration");
TH1F* mrdzlayerloc = new TH1F(mrdzlayer);
mrdzlayerloc->SetDirectory(0);
canvloc3->cd();
mrdzlayerloc->SetStats(0);
mrdzlayerloc->Draw();
*/

// histogram binning check
//cout<<"xmin: "<<(mrdzlayer.GetXaxis()->GetXmin())<<endl;
//cout<<"xmax: "<<(mrdzlayer.GetXaxis()->GetXmax())<<endl;
Int_t mostbins = mrdzlayer.GetNbinsX();
/*cout<<"numbins: "<<mostbins<<endl;
for(Int_t i=1;i<(mostbins+1);i++){	// NOTE BIN NUMBERING CONVENTION: BIN 0 = UNDERFLOW BIN, BIN NBINS+1 = OVERFLOW BIN
	Double_t binlowedge = mrdzlayer.GetXaxis()->GetBinLowEdge(i);
	Double_t binupedge = mrdzlayer.GetXaxis()->GetBinUpEdge(i);
	cout<<"bin low edge: "<<binlowedge<<" up edge: "<<binupedge<<endl;
}*/
Double_t underconts = mrdzlayer.GetBinContent(0);
Double_t overconts = mrdzlayer.GetBinContent(mostbins+1);
cout<<"hits in underflow bin: "<<underconts<<" and overflow bin: "<<overconts<<endl;
// How can we plot this as an extra bin on the histogram?

// histogram filling check: Let's check the maximum z.
//Steps:
//1) get number of hits this entry
//2) create dynamic array of appropriate size for all z positions of hits
//3) get z position branch into that array
//4) scan through the array and increment zmax found so far.
/*
Int_t numentries = mrdtree->GetEntries();
cout<<"numentries:"<<numentries<<endl;
TBranch* numhitsbranch = mrdtree->GetBranch("hitnum");
numhitsbranch->SetAutoDelete(kTRUE);
Int_t arraysize = 0;
numhitsbranch->SetAddress(&arraysize);
TBranch* zbranch = mrdtree->GetBranch("mrdhit_z");
zbranch->SetAutoDelete(kTRUE); //can we omit this if we defined event as static?
Double_t* hitzvals;
Double_t zval = 0;
Double_t zvalmax=0;
for (Int_t i=0;i<numentries;i++){
	numhitsbranch->GetEntry(i);
	if(arraysize!=0){ 
		hitzvals = new Double_t[arraysize];
		zbranch->SetAddress(hitzvals);
		zbranch->GetEntry(i);
		for(size_t j=0;j<arraysize;j++){
			zval=hitzvals[j];
			if(zval>zvalmax){zvalmax=zval;}
		}
		delete[] hitzvals;
		hitzvals=0;
	}
}
cout<<"maximum z of a hit in the MRD: "<<zvalmax<<endl;
//result: max z is 2837.52mm = 847mm into the MRD = 0.922 of the total depth = 11.9/13 th of the depth = i.e. in plate 12 

mrdtree->ResetBranchAddresses();
*/

}


void EnergyDepositioninLayers(TTree* mrdtree){
// Produce plot of energy deposition vs depth, with bins corresponding to MRD layers
// ==================================================================================
//Steps:
//1) get number of hits this entry
//2) create dynamic array of appropriate size for all z positions of hits
//3) get z position branch into that array
//4) scan through the array and increment zmax found so far.
TCanvas edeposCanv = TCanvas("edeposCanv","Title");
Int_t numentries = mrdtree->GetEntries();
TBranch* numhitsbranch = mrdtree->GetBranch("hitnum");
Int_t arraysize = 0;
numhitsbranch->SetAddress(&arraysize);
TBranch* zbranch = mrdtree->GetBranch("mrdhit_z");
//zbranch->SetAutoDelete(kTRUE);	
// this will ensure when each GetEntry is called, the previous array used for the branch will be deleted first before being reused
Double_t* hitzvals;
Double_t zval = 0;
TBranch* edepbranch = mrdtree->GetBranch("mrdhit_edep");
//edepbranch->SetAutoDelete(kTRUE);
Double_t* hitedeps;				// read from tree
Double_t edep = 0;
Double_t* edeps = new Double_t[mrdnumlayers];	// summed for graph
Double_t* layercentres = new Double_t[mrdnumlayers];
for(Int_t layer=1;layer<(mrdnumlayers+1);layer++){layercentres[layer-1]=(Double_t)layer;edeps[layer-1]=0.;}//cout<<"layer: "<<layer<<endl;
Int_t overflowcount=0;
Int_t underflowcount=0;
Int_t mrdhitcount=0;
for (Int_t i=0;i<numentries-1;i++){
	numhitsbranch->GetEntry(i);
	if(arraysize!=0){
		hitzvals = new Double_t[arraysize];
		zbranch->SetAddress(hitzvals);
		zbranch->GetEntry(i);
		hitedeps = new Double_t[arraysize];
		edepbranch->SetAddress(hitedeps);
		edepbranch->GetEntry(i);
		for(size_t j=0;j<arraysize;j++){
			mrdhitcount++;
			//cout<<"j="<<j<<endl;
			zval=hitzvals[j];
			edep=hitedeps[j];
			for(Int_t k=0;k<mrdnumlayers;k++){
				//cout<<"k="<<k<<endl;
				Double_t layerstart = mrdzstart+((mrdzlen/mrdnumlayers))*k;
				Double_t layerend = mrdzstart+((mrdzlen/mrdnumlayers))*(k+1);
				//cout<<"checking if z "<<zval<<" lies in "<<layerstart<<" to "<<layerend<<endl;
				if(zval>=layerstart&&zval<layerend){edeps[k]+=edep;break;}
			}
			// skipped by break statement if binned 
			if (zval>=(mrdzstart+mrdzlen)){overflowcount++;}	
			else if (zval<mrdzstart){underflowcount++;}
		}
		//cout<<"deleting"<<endl;
		delete[] hitzvals; hitzvals=0; 
		delete[] hitedeps; hitedeps=0;
	}
	//cout<<"i="<<i<<endl;
}

//for(Int_t i=0;i<mrdnumlayers;i++){cout<<"layer: "<<layercentres[i]<<" edep: "<<edeps[i]<<endl;}
//cout<<overflowcount<<" overflows and "<<underflowcount<<" underflows in manual Edep binning"<<endl;
cout<<"total of "<<mrdhitcount<<" hits, vs "<<mrdtree->GetEntries()<<" entries"<<endl;
edeposCanv.cd();
edeposCanv.Clear();
TGraph edepvsz = TGraph(mrdnumlayers,layercentres,edeps);
//https://root.cern.ch/root/html534/guides/users-guide/Graphs.html#graph-draw-options
edepvsz.SetFillColor(40); 	//light blue
edepvsz.Draw("AB");		// axes and bar graphs - Axes are not drawn unless requested!?
edepvsz.SetTitle("Total Energy Deposited in each MRD Layer");
edepvsz.GetXaxis()->SetTitle("MRD Layer");
edepvsz.GetXaxis()->SetRangeUser(0,mrdnumlayers+0.5);
edepvsz.GetXaxis()->CenterTitle(true);
edepvsz.GetXaxis()->SetNdivisions(15);
edepvsz.GetYaxis()->SetTitle("Total Energy Deposit (keV)");
edepvsz.GetYaxis()->CenterTitle(true);
edeposCanv.Modified();
// Note: You can only set axes titles AFTER drawing - until a graph is drawn it does not HAVE axes!
//https://root.cern.ch/root/html534/guides/users-guide/Graphs.html#setting-the-graphs-axis-title
//edepvsz.Write(); method 2 is more efficient. 
edeposCanv.SaveAs("edepvsz.png");
//** local print for check
/*TCanvas* canvloc4 = new TCanvas("canvloc4","Energy Deposited in MRD Layers");
TGraph* edepvszlocal = new TGraph(edepvsz);//mrdnumlayers,layercentres,edeps);
canvloc4->cd();
canvloc4->Clear();
edepvszlocal->Draw("AB");
edepvszlocal->GetXaxis()->SetTitle("MRD Layer");
edepvszlocal->GetXaxis()->CenterTitle(true);
edepvszlocal->GetXaxis()->SetNdivisions(15);
edepvszlocal->GetYaxis()->SetTitle("Total Energy Deposit (keV)");
edepvszlocal->GetYaxis()->CenterTitle(true);
canvloc4->Modified();*/

mrdtree->ResetBranchAddresses();

//delete[] hitzvals; hitzvals=0; // do for last loop
//delete[] hitedeps; hitedeps=0; // do for last loop
delete[] edeps;
delete[] layercentres;

}

void EnergyDepositionInPaddles(TTree* mrdtree){
// Produce Distributions of Energy Deposition in each Panel, with bins corresponding to Scintillator Paddles
// =========================================================================================================
/* Steps
setalias of panelnum as floor(copynum%15).
setalias of paddlenum as copynum - (panelnum*30).
plot paddlenum on x weighting by edep with cut panelnum = each. 
*/
TCanvas edeposCanv2 = TCanvas("edeposCanv2","Title",8,59,1266,737);
//TCanvas* canvloc5 = new TCanvas("canvloc5","Energy deposited in MRD panels",8,59,1266,737);
//canvloc5->Divide(4,3);	// automatically creates canvases canvloc5_1 to canvloc5_12. 
edeposCanv2.Clear();
edeposCanv2.Divide(4,3);	// hard coded reliant on mrdnumlayers. 
Double_t* totedepinpanels = new Double_t[mrdnumlayers];
mrdtree->SetAlias("mrdhit_zlayer",TString::Format("TMath::Floor(((mrdhit_z-%f)/%f)*%d)+1",mrdzstart,mrdzlen,mrdnumlayers));
mrdtree->SetAlias("mrdhit_paddlenum",TString::Format("(mrdhit_copynum+1)-((mrdhit_zlayer*%d)-1)",mrdpaddlesperpanel));
for(Int_t panelnum=1;panelnum<(mrdnumlayers+1);panelnum++){
	edeposCanv2.cd(panelnum);
	mrdtree->Draw(TString::Format("mrdhit_paddlenum>>edepinpanel%d(%d,1,%d)",panelnum,mrdpaddlesperpanel,mrdpaddlesperpanel+1),TString::Format("mrdhit_edep*(mrdhit_zlayer==%d)",panelnum));
	// note that 'cut' is actually a weighting. For boolean expressions this becomes a cut, with weighting 1 or 0. 
	// Here we specify both a weighting and a cut by multiplying the weight by the boolean.
	TCanvas* currentcanvas = (TCanvas*)edeposCanv2.GetPrimitive(TString::Format("edeposCanv2_%d",panelnum));
	TH1F* edepinpanel = (TH1F*)currentcanvas->GetPrimitive(TString::Format("edepinpanel%d",panelnum));
	edepinpanel->SetTitle(TString::Format("Energy Deposition in MRD Layer %d;Paddle Number;Total Energy Deposited",panelnum));
	edepinpanel->GetXaxis()->CenterLabels();
	rootfileout->cd();
	edepinpanel->Write();
	Double_t totedepinpanel = edepinpanel->Integral();
	totedepinpanel+=(edepinpanel->GetBinContent(0))+(edepinpanel->GetBinContent(edepinpanel->GetNbinsX()+1));
	// this is still needed because copynum does not seem to properly reflect paddle/panel num. Copynums > 40 appear for hits with z in layer 1 as shown by Scan below.
	totedepinpanels[panelnum-1]=totedepinpanel;
	/*if(panelnum==1){
		mrdtree->Scan("mrdhit_z:mrdhit_zlayer:mrdhit_copynum:mrdhit_paddlenum:mrdhit_edep","mrdhit_z<=2066.58","",500,1);
	}*/
	//** local print for check
	/*TH1F* oneofthehistos = new TH1F(*edepinpanel);
	oneofthehistos->SetDirectory(0);
	canvloc5->cd(panelnum+1);
	oneofthehistos->Draw();*/
}
edeposCanv2.cd();
edeposCanv2.Modified();
edeposCanv2.SaveAs("deposition in panels.png");
edeposCanv2.Clear();
/*mrdtree->Draw("mrdhit_paddlenum>>paddlenumoverflow",Form("mrdhit_paddlenum>%d",mrdpaddlesperpanel));
TH1F* paddlenumoverflow = (TH1F*)edeposCanv2.GetPrimitive("paddlenumoverflow");
Int_t paddnumoverflowcount = paddlenumoverflow->GetEntries();
mrdtree->Draw("mrdhit_zlayer>>panelnumoverflow",Form("mrdhit_zlayer>%d",mrdnumlayers));
TH1F* panelnumoverflow = (TH1F*)edeposCanv2.GetPrimitive("panelnumoverflow");
Int_t panelnumoverflowcount = panelnumoverflow->GetEntries();
cout<<paddnumoverflowcount<<" overflows in alias paddle numbering"<<endl;
cout<<panelnumoverflowcount<<" overflows in alias panel numbering"<<endl;*/


TCanvas edeposCanv22 = TCanvas("edeposCanv22","Title");
Double_t* layercentres2 = new Double_t[mrdnumlayers];
for(Int_t layer=1;layer<(mrdnumlayers+1);layer++){layercentres2[layer-1]=(Double_t)layer;}
TGraph edepvsz2 = TGraph(mrdnumlayers,layercentres2,totedepinpanels);
edepvsz2.SetFillColor(46); 	//light red
edepvsz2.Draw("AB");		// axes and bar graphs
edepvsz2.SetTitle("Total Energy Deposited in each MRD Layer 2");
edepvsz2.GetXaxis()->SetTitle("MRD Layer");
edepvsz2.GetXaxis()->SetRangeUser(0,mrdnumlayers+0.5);
edepvsz2.GetXaxis()->CenterTitle(true);
edepvsz2.GetXaxis()->SetNdivisions(15);
edepvsz2.GetYaxis()->SetTitle("Total Energy Deposit (keV)");
edepvsz2.GetYaxis()->CenterTitle(true);
edeposCanv22.Modified();
edepvsz2.Write();
edeposCanv22.SaveAs("edepvsz2.png");
//Double_t manualintegral = 0;
//for(Int_t i=0;i<(edepvsz2.GetN());i++){manualintegral=manualintegral+edepvsz2.GetY()[i];}

delete[] totedepinpanels;
delete[] layercentres2;

}

void splitindividualsources(TTree* mrdtree, TTree* recotree){
// Split hits from each event into groups of hits corresponding to a single 'track' - in case we have more than one track per event
// ================================================================================================================================
TCanvas splittracksCanv = TCanvas("splittracksCanv","Title");
Int_t numentries = mrdtree->GetEntries();
Int_t numhitsthisevent=0;
Int_t numhitsthiseventtemp=0;
TBranch* numhitsbranch = mrdtree->GetBranch("hitnum");
numhitsbranch->SetAddress(&numhitsthiseventtemp);
TBranch* hittimebranch = mrdtree->GetBranch("mrdhit_t");
Double_t* hittimes; 
TBranch* trackidsbranch = mrdtree->GetBranch("mrdhit_trackID");
Int_t* trackids;

TH1F* splittrackshist=0;

Int_t numtracksthisevent;
mrdtree->Draw("Entry$>>hitnumhist","hitnum");	// 'hitnum' branch is number of hits in each event.
TH1F* hitnumhist = (TH1F*)splittracksCanv.GetPrimitive("hitnumhist");
Int_t maxtracksperevent = hitnumhist->GetMaximum();
Int_t* originalhitnum = new Int_t[maxtracksperevent];	// index for retrieval of all properties from original MRD tree.
Int_t* tracknumthisevent = new Int_t[maxtracksperevent];
TBranch* tracksthiseventb = recotree->Branch("numtracksthisevent",&numtracksthisevent);
TBranch* hitsthiseventb = recotree->Branch("numhitsthisevent",&numhitsthisevent);
TBranch* originalhitnumb = recotree->Branch("originalhitnum",originalhitnum,"originalhitnum[numhitsthisevent]/I");
TBranch* tracknumthiseventb = recotree->Branch("tracknumthisevent",tracknumthisevent,"tracknumthisevent[numhitsthisevent]/I");

// to be made a friend tree of mrdtree
cout<<"scanning over "<<numentries<<" entries."<<endl;
for(Int_t entry=0;entry<numentries;entry++){
	numhitsbranch->GetEntry(entry); // first check if we have any hits this entry. 
	numhitsthisevent = numhitsthiseventtemp;	// numhitsthiseventtemp WILL GET OVERWRITTEN BY THE DRAW() CALL BELOW!
	if(numhitsthisevent==0){
		numtracksthisevent=0;
		originalhitnum[0]=0;
		tracknumthisevent[0]=0;
		//recotree->Fill();	// fill the tree anyway so entries align
		tracksthiseventb->Fill();
		hitsthiseventb->Fill();
		originalhitnumb->Fill();
		tracknumthiseventb->Fill();
		continue;
	}	// skip remainder
	//cout<<"entry "<<entry<<" has "<<numhitsthisevent<<" hits."<<endl;
	// if continued to here, we have at least some hits. 
	//mrdtree->ResetBranchAddresses();	// does not like this (unless all branches have addresses?)
	hittimebranch->ResetAddress(); 	// so we don't have to ensure sufficient memory allocation for 'Draw()' call.
	trackidsbranch->ResetAddress();
	// isolate each track by hit time - multiple hits from the same track generally all occur within 1ns span (1 unit in G4 native)
	// need to plot histogram with sufficient resolution that there will be at least one empty bin between
	// hits from adjacent tracks
	if(splittrackshist!=0){delete splittrackshist; splittrackshist=0;}
	mrdtree->Draw("mrdhit_t>>splittrackshist",TString::Format("Entry$==%d",entry));
	splittrackshist = (TH1F*)splittracksCanv.GetPrimitive("splittrackshist");
	// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
	Double_t eventendtime = splittrackshist->GetXaxis()->GetXmax();
	Double_t eventstarttime = splittrackshist->GetXaxis()->GetXmin();
	Double_t eventduration = (eventendtime - eventstarttime);
	//cout<<"event start: "<<eventstarttime<<" end : "<<eventendtime<<endl;
	//cout<<"event duration is "<<eventduration<<endl;
	if(eventduration<maxtrackduration){
		//cout<<"all hits this event within one track."<<endl;
		// all hits are within a short period of time - assume they are from the same track.
		numtracksthisevent=1;
		for(Int_t j=0;j<numhitsthisevent;j++){originalhitnum[j]=j;tracknumthisevent[j]=0;}
		//recotree->Fill();
		tracksthiseventb->Fill();
		hitsthiseventb->Fill();
		originalhitnumb->Fill();
		tracknumthiseventb->Fill();
	} else {
		// this event has multiple tracks. Need to split hits into which track they belong to.
		// first, to keep track of what's been assigned, initialize all to -1.		
		for(Int_t j=0;j<numhitsthisevent;j++){tracknumthisevent[j]=-1;}
		
		// Now ensure binning is fine enough to have at least one empty bin between hits.
		Int_t numbins = (TMath::Floor(eventduration))/5;
		if(numbins<100){numbins=100;}	// a minimum of 100 bins - too few bins and ROOT clips the range, omitting points.
		if(splittrackshist!=0){delete splittrackshist;splittrackshist=0;}
		mrdtree->Draw(TString::Format("mrdhit_t>>splittrackshist(%d,%f,%f)",numbins,eventstarttime,eventendtime),TString::Format("Entry$==%d",entry));
		splittrackshist = (TH1F*)splittracksCanv.GetPrimitive("splittrackshist");
		//cout<<"using "<<numbins<<", "<<splittrackshist->GetNbinsX()<<" bins of width "<<splittrackshist->GetXaxis()->GetBinWidth(1)<<endl;

		/* couldn't get this to work reliably.
		TSpectrum *sp1 = new TSpectrum();
		sp1->Search(splittrackshist,1
		.0.,"nobackground",0.01);	// Search(TH1*,width, option, threshold)
		numtracksthisevent = sp1->GetNPeaks();
      		cout<<numtracksthisevent<<" tracks this event"<<endl;
    		Float_t* trackhittimes = sp1->GetPositionX();
    		// do we need to convert to vector? return is probably already ordered...
  		std::vector<Float_t> trackhittimesv(trackhittimes, trackhittimes+numtracksthisevent); // convert to vector so we can...
  		std::sort(trackhittimesv.begin(), trackhittimesv.end());           // sort ascending
  		*/
  		bool inpeak=false;
  		std::vector<Float_t> trackhittimesv;
  		for(Int_t i=0;i<numbins;i++){
  			if(splittrackshist->GetBinContent(i)>0){inpeak=true;}
  			else {if(inpeak==true){
  				trackhittimesv.push_back(splittrackshist->GetXaxis()->GetBinLowEdge(i));
  				//cout<<"Setting track time threshold at "<<trackhittimesv.back()<<endl;
  				inpeak=false;}
  			}
  		}
  		numtracksthisevent = trackhittimesv.size();
      		//cout<<numtracksthisevent<<" tracks this event"<<endl;
  		
		hittimes = new Double_t[numhitsthisevent];
		hittimebranch->SetAddress(hittimes);
		hittimebranch->GetEntry(entry);
		trackids = new Int_t[numhitsthisevent];
		trackidsbranch->SetAddress(trackids);
		trackidsbranch->GetEntry(entry);
		
		for(Int_t j=0; j<numtracksthisevent; j++){
			//cout<<"Track struck MRD at = "<<trackhittimesv.at(j)<<"ns in event "<<entry<<endl;
			// Search will find peak in hit time, but we need to encapsulate all hits, so convert to bin number, then 
			// scan forward from the peak bin until we have an empty bin. 
			// don't need to worry about lower bound as we start from lowest t peak and exclude already allocated hits
			/* associated with Spectrum method - not needed with new manual method
			Int_t binnumber = splittrackshist->FindBin(trackhittimesv.at(j));
			while(splittrackshist->GetBinContent(binnumber+1)!=0){binnumber++;}	
			Double_t trackmaxtime = splittrackshist->GetXaxis()->GetBinUpEdge(binnumber);
			cout<<"Setting track time threshold at "<<trackmaxtime<<endl;
			*/
			
			// hit times are not ordered, so scan through them all 
			for(Int_t thishit=0;thishit<numhitsthisevent;thishit++){
				if(tracknumthisevent[thishit]<0&&hittimes[thishit]<trackhittimesv.at(j)){
					originalhitnum[thishit]=thishit;
					tracknumthisevent[thishit]=j;
					//if(entry==5){cout<<"put hit "<<thishit<<" into track "<<j<<" of this event";
					//cout<<" original trackID = "<<trackids[thishit]<<" hit time: "<<hittimes[thishit]<<endl;}
				}
			}
		}
		//cout<<"next event"<<endl;
		// quick scan to check for any unallocated hits (i.e. in a peak not detected by TSpectrum search)
		for(Int_t k=0;k<numhitsthisevent;k++){if(tracknumthisevent[k]==-1){cout<<"*****unbinned hit!"<<k<<" "<<hittimes[k]<<endl;}}
		//cout<<"filling tree event "<<entry<<" with "<<numhitsthisevent<<" hits"<<endl;
		//recotree->Fill();
		tracksthiseventb->Fill();
		hitsthiseventb->Fill();
		originalhitnumb->Fill();
		tracknumthiseventb->Fill();
		//delete sp1;
		delete[] hittimes;
		delete[] trackids;
	}	// next hit in entry.
}	// next entry in mrdtree

rootfileout->cd();
recotree->Write();

delete[] originalhitnum;
delete[] tracknumthisevent;
mrdtree->ResetBranchAddresses();

}		


void splitPMTtracks(TTree* pmttree, TTree* recotree){
// Split hits from each event into groups of hits corresponding to a single 'track' - in case we have more than one track per event
// ================================================================================================================================
TCanvas splittracksCanv = TCanvas("splittracksCanv","Title");
Int_t numentries = pmttree->GetEntries();
Int_t numhitsthisevent=0;
Int_t numhitsthiseventtemp=0;
TBranch* numhitsbranch = pmttree->GetBranch("pmt_hitnum");
numhitsbranch->SetAddress(&numhitsthiseventtemp);
TBranch* hittimebranch = pmttree->GetBranch("pmthit_t");
Double_t* hittimes; 
TBranch* trackidsbranch = pmttree->GetBranch("pmthit_trackID");	// only used for checking
Int_t* trackids;
TBranch* copynumbranch = pmttree->GetBranch("pmthit_copynum");
Int_t* copynums;

TH1F* splittrackshist=0;

Int_t numtracksthisevent;
pmttree->Draw("Entry$>>hitnumhist","pmt_hitnum");	// 'hitnum' branch is number of hits in each event.
TH1F* hitnumhist = (TH1F*)splittracksCanv.GetPrimitive("hitnumhist");
Int_t maxtracksperevent = hitnumhist->GetMaximum();

//cMRDTrack* aTrack=0;
Bool_t* hitunallocated = new Bool_t[maxtracksperevent];
//cMRDTrack* aTrack = new cMRDTrack[maxtracksperevent];

TClonesArray* aTrack = new TClonesArray("cMRDTrack");	// string is class name
TClonesArray &aTracka = *aTrack;
std::vector<Double_t> hittimesthisstrike;
std::vector<cMRDStrike> strikesthistrack;

TBranch* pmtnumtracksthiseventb = recotree->Branch("pmtnumtracksthisevent",&numtracksthisevent);
TBranch* pmtnumhitsthiseventb = recotree->Branch("pmtnumhitsthisevent",&numhitsthisevent);
TBranch* pmtTracksBranch = recotree->Branch("PMTrecotracks",&aTrack, maxtracksperevent); 
//TBranch* pmtTracksBranch  = recotree->Branch("PMTrecotracks[pmtnumtracksthisevent]/cMRDTrack",&aTrack);

std::vector<Double_t>* hittimesbystrike = new std::vector<Double_t>[mrdnumlayers];
std::vector<Int_t>* copynumsbystrike = new std::vector<Int_t>[mrdnumlayers];

cout<<"scanning over "<<numentries<<" entries."<<endl;
for(Int_t entry=0;entry<numentries;entry++){
	aTracka.Clear();	
	// if you class contains pointers, use aTrack.Clear("C"). You MUST then provide a Clear() method in your class that properly
	// performs clearing and memory freeing. (or "implements the reset procedure for pointer objects") 
	// https://root.cern.ch/doc/master/classTClonesArray.html#a025645e1e80ea79b43a08536c763cae2
	numhitsbranch->GetEntry(entry); // first check if we have any hits this entry. 
	numhitsthisevent = numhitsthiseventtemp;	// numhitsthiseventtemp WILL GET OVERWRITTEN BY THE DRAW() CALL BELOW!
	if(numhitsthisevent==0){
		numtracksthisevent=0;
		new(aTracka[0]) cMRDTrack();
		//recotree->Fill();	// fill the tree anyway so entries align
		pmtnumtracksthiseventb->Fill();
		pmtnumhitsthiseventb->Fill();
		pmtTracksBranch->Fill();
		continue;
	}	// skip remainder
	//cout<<"entry "<<entry<<" has "<<numhitsthisevent<<" hits."<<endl;
	// if continued to here, we have at least some hits. 
	hittimebranch->ResetAddress(); 	// so we don't have to ensure sufficient memory allocation for 'Draw()' call.
	trackidsbranch->ResetAddress();
	copynumbranch->ResetAddress();
	// isolate each track by hit time - multiple hits from the same track generally all occur within 1ns span (1 unit in G4 native)
	// need to plot histogram with sufficient resolution that there will be at least one empty bin between
	// hits from adjacent tracks
	if(splittrackshist!=0){delete splittrackshist; splittrackshist=0;}
	pmttree->Draw("pmthit_t>>splittrackshist",TString::Format("Entry$==%d",entry));
	splittrackshist = (TH1F*)splittracksCanv.GetPrimitive("splittrackshist");
	// first check: are all hits within a 30ns window (maxtrackduration) If so, just one track. 
	Double_t eventendtime = splittrackshist->GetXaxis()->GetXmax();
	Double_t eventstarttime = splittrackshist->GetXaxis()->GetXmin();
	Double_t eventduration = (eventendtime - eventstarttime);
	//cout<<"event start: "<<eventstarttime<<" end : "<<eventendtime<<endl;
	//cout<<"event duration is "<<eventduration<<endl;
	if(eventduration<maxpmttrackduration){
		//cout<<"all hits this event within one track."<<endl;
		// all hits are within a short period of time - assume they are from the same track.
		numtracksthisevent=1;
		hittimesthisstrike.clear();
		hittimes = new Double_t[numhitsthisevent];
		hittimebranch->SetAddress(hittimes);
		hittimebranch->GetEntry(entry);
		hittimesthisstrike.assign(hittimes,hittimes+numhitsthisevent);
		copynums = new Int_t[numhitsthisevent];	
		copynumbranch->SetAddress(copynums);
		copynumbranch->GetEntry(entry);
		Int_t PMTmosthit = findmode(copynums, numhitsthisevent);
		cMRDStrike aStrike = cMRDStrike(hittimesthisstrike,PMTmosthit);
		strikesthistrack.clear();
		strikesthistrack.push_back(aStrike);
		new(aTracka[0]) cMRDTrack(strikesthistrack);	
		// can also use 'cMRDTrack* = (cMRDTrack*)aTrack.ConstructedAt(0);' followed by a bunch of 'Set' calls to 
		// set all relevant fields. This bypasses the constructor, calling it only when necessary, saving time.
		// in this case we do not need to call aTracka.Clear(); 
		//recotree->Fill();
		pmtnumtracksthiseventb->Fill();
		pmtnumhitsthiseventb->Fill();
		pmtTracksBranch->Fill();

	} else {
		// this event has multiple tracks. Need to split hits into which track they belong to.
		// first, to keep track of what's been assigned, initialize all to -1.		
		for(Int_t j=0;j<numhitsthisevent;j++){hitunallocated[j]=true;}
		
		// Now ensure binning is fine enough to have at least one empty bin between tracks.
		Int_t numbins = TMath::Floor(eventduration)/50;
		if(splittrackshist!=0){delete splittrackshist;splittrackshist=0;}
		pmttree->Draw(TString::Format("pmthit_t>>splittrackshist(%d,%f,%f)",numbins,eventstarttime,eventendtime),TString::Format("Entry$==%d",entry));
		// unless we explicitly give the start/end times ROOT decides to arbitrarily clip the range when specifying # bins
		splittrackshist = (TH1F*)splittracksCanv.GetPrimitive("splittrackshist");
		//cout<<"using "<<numbins<<", "<<splittrackshist->GetNbinsX()<<" bins of width "<<splittrackshist->GetXaxis()->GetBinWidth(1)<<endl;

  		bool inpeak=false;
  		std::vector<Float_t> trackhittimesv;
  		for(Int_t i=0;i<numbins;i++){
  			if(splittrackshist->GetBinContent(i)>0){inpeak=true;}
  			else {if(inpeak==true){
  				trackhittimesv.push_back(splittrackshist->GetXaxis()->GetBinLowEdge(i));
  				//cout<<"Setting track time threshold at "<<trackhittimesv.back()<<endl;
  				inpeak=false;}
  			}
  		}
  		numtracksthisevent = trackhittimesv.size();
      		cout<<numtracksthisevent<<" tracks and "<<numhitsthisevent<<" hits in event "<<entry<<endl;
  		
		hittimes = new Double_t[numhitsthisevent];
		hittimebranch->SetAddress(hittimes);
		hittimebranch->GetEntry(entry);
		trackids = new Int_t[numhitsthisevent];	// only used for checking
		trackidsbranch->SetAddress(trackids);
		trackidsbranch->GetEntry(entry);
		copynums = new Int_t[numhitsthisevent];
		copynumbranch->SetAddress(copynums);
		copynumbranch->GetEntry(entry);

		// for each track time threshold, scan through all hits and collect those in this track
		// collect hits in a given track into an array of vectors; each array element respresents an MRD layer/potential strike
		for(Int_t layer=0;layer<mrdnumlayers;layer++){
			hittimesbystrike[layer].clear();
			copynumsbystrike[layer].clear();
		}

		for(Int_t j=0; j<numtracksthisevent; j++){
			strikesthistrack.clear();
			//cout<<"processing track "<<j<<endl;
			// hit times are not ordered, so scan through them all 
			for(Int_t thishit=0;thishit<numhitsthisevent;thishit++){
				if(hitunallocated[thishit]&&hittimes[thishit]<trackhittimesv.at(j)){
					//cout<<"allocating hit "<<thishit<<endl;
					hitunallocated[thishit]=false;
					// timing gives the track - now check z value to determine appropriate strike
					Int_t panelnum = TMath::Floor(copynums[thishit]/mrdpaddlesperpanel);
					hittimesbystrike[panelnum].push_back(hittimes[thishit]);
					copynumsbystrike[panelnum].push_back(copynums[thishit]);
					//if(entry==94){cout<<"put hit "<<thishit<<" into track "<<j<<" of this event;";
					//cout<<" hit time: "<<hittimes[thishit]<<endl;}
				}
			}

			// we have an array of 12 potential strikes; check which have hits and for those with >0 hits,
			// create a strike object.
			for(Int_t layer=0;layer<mrdnumlayers;layer++){
				//cout<<"forming strike for layer "<<layer<<endl;
				// each layer with hits in is a strike
				if(hittimesbystrike[layer].size()>0){
					cout<<hittimesbystrike[layer].size()<<" hits in strike in layer "<<layer;
					// create a strike for this layer. 
					hittimesthisstrike.clear();
					hittimesthisstrike=hittimesbystrike[layer];
					std::vector<Int_t> copynumsthisstrike = copynumsbystrike[layer];
					Int_t PMTmosthit = findmode(&copynumsthisstrike.front(),hittimesthisstrike.size());
					cout<<", PMT most hit: "<<PMTmosthit<<endl;
					cMRDStrike aStrike = cMRDStrike(hittimesthisstrike,PMTmosthit);
					strikesthistrack.push_back(aStrike);
				} //else {cout<<"no strikes this layer"<<endl;}
			}
			new(aTracka[j]) cMRDTrack(strikesthistrack);
			cout<<"creating track of "<<strikesthistrack.size()<<" strikes."<<endl;
		}
					
		// quick scan to check for any unallocated hits (i.e. in a peak not detected by TSpectrum search)
		for(Int_t k=0;k<numhitsthisevent;k++){if(hitunallocated[k]){cout<<"*****unbinned hit!"<<k<<" "<<hittimes[k]<<endl;}}
		//cout<<"filling tree event "<<entry<<" with "<<numhitsthisevent<<" hits"<<endl;
		//recotree->Fill();
		pmtnumtracksthiseventb->Fill();
		pmtnumhitsthiseventb->Fill();
		pmtTracksBranch->Fill();
		delete[] hittimes;
		delete[] trackids;
		delete[] copynums;
	}	// more than one track
	//cout<<"moving to next entry"<<endl; // doesn't show for every loop because of 'continue' statement above
}	// next entry in pmttree
rootfileout->cd();
recotree->Write("",TObject::kOverwrite);
delete[] hittimesbystrike;	// using delete rather than delete[] results in segfault!
delete[] copynumsbystrike;
pmttree->ResetBranchAddresses();
// NOTE: WHEN STORING AN OBJECT IN A BRANCH: The object must not be destroyed (i.e. be deleted) until the TTree is deleted or TTree::ResetBranchAddress is called. 
cout<<"resetting recotree"<<endl;
recotree->ResetBranchAddresses(); // so this MUST be called BEFORE aTrack gets deleted!
//cout<<"deleting aTrack"<<endl;	// doesn't seem necessary for TClonesArray? 
//delete[] aTrack;		// using delete rather than delete[] results in segfault!
cout<<"exiting"<<endl;

}		
	
void checkpmtclass(TTree* recotree){
cout<<"********** STARTING TREE CHECKS **********"<<endl;
//cMRDTrack* recotracks;
TClonesArray* recotracks = new TClonesArray("cMRDTrack");
TBranch* recotracksbranch = recotree->GetBranch("PMTrecotracks");
recotracksbranch->SetAutoDelete(kFALSE);
recotracksbranch->SetAddress(&recotracks);
Int_t numrecotracks;
TBranch* numrecotracksbranch = recotree->GetBranch("pmtnumtracksthisevent");
numrecotracksbranch->SetAddress(&numrecotracks);

for(Int_t event=0;event<recotree->GetEntries();event++){
	cout<<"getting event "<<event<<endl;
	numrecotracksbranch->GetEntry(event);
	cout<<numrecotracks<<" reconstructed tracks this event"<<endl;
	if(numrecotracks>0){
		//recotracks = new cMRDTrack[numrecotracks];
		recotracks->Clear();
		cout<<"getting entry"<<endl;
		recotracksbranch->GetEntry(event);
		cout<<"retrieved reconstructed tracks branch."<<endl;
		cout<<"num tracks in clonesarray: "<<recotracks->GetEntriesFast()<<endl;
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cout<<"track "<<tracknum<<" has particle code: "<<((cMRDTrack*)recotracks->At(tracknum))->GetParticleID()<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cout<<"track "<<tracknum<<" has "<<((cMRDTrack*)recotracks->At(tracknum))->GetNumStrikes()<<" strikes"<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cout<<"track "<<tracknum<<" has first strike at "<<((cMRDTrack*)recotracks->At(tracknum))->GetFirstStrike()<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cout<<"track "<<tracknum<<" has last strike at "<<((cMRDTrack*)recotracks->At(tracknum))->GetLastStrike()<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cMRDStrike* firststrike = ((cMRDTrack*)recotracks->At(tracknum))->GetFirstStrike();
			cout<<"track "<<tracknum<<" first strike struck PMT "<<firststrike->GetPMTNumber()<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cMRDStrike* firststrike = ((cMRDTrack*)recotracks->At(tracknum))->GetFirstStrike();
			cout<<"track "<<tracknum<<" first strike struck MRD layer "<<firststrike->GetLayerNum()<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cMRDStrike* firststrike = ((cMRDTrack*)recotracks->At(tracknum))->GetFirstStrike();
			std::vector<Double_t> hittimes = firststrike->GetPMThitTimes();
			cout<<"track "<<tracknum<<" first strike has "<<hittimes.size()<<" hits"<<endl;
		}
		for(Int_t tracknum=0;tracknum<numrecotracks;tracknum++){
			cMRDStrike* firststrike = ((cMRDTrack*)recotracks->At(tracknum))->GetFirstStrike();
			std::vector<Double_t> hittimes = firststrike->GetPMThitTimes();
			for(Int_t hit=0;hit<5;hit++){
				cout<<"track "<<tracknum<<" first strike hit "<<hit<<" is at "<<hittimes.at(hit)<<"ns"<<endl;
			}
			for(Int_t hit=(hittimes.size()-5);hit<hittimes.size();hit++){
				cout<<"track "<<tracknum<<" first strike hit "<<hit<<" is at "<<hittimes.at(hit)<<"ns"<<endl;
			}
		}
		cout<<"deleting recotracks allocated array"<<endl;
		//delete[] recotracks;
	}
}
cout<<"resetting branch addresses"<<endl;
recotree->ResetBranchAddresses();
cout<<"done"<<endl;
}


void countmrdhits(TTree* fulltree){
// Calculate MRD muon stopping efficiency
// ======================================
TCanvas countCanv = TCanvas("countCanv","Muon Counts");
countCanv.cd();
fulltree->Draw("part_trackid>>muons","part_pid==13||part_pid==-13");
TH1F* muons = (TH1F*)countCanv.GetPrimitive("muons");
Int_t allmuons = muons->GetEntries();
fulltree->Draw("part_trackid>>muons","part_mrdhit==1&&(part_pid==13||part_pid==-13)");
muons = (TH1F*)countCanv.GetPrimitive("muons");
Int_t hitmuons = muons->GetEntries();
fulltree->Draw("part_trackid>>muons","part_mrdhit==1&&part_processEnd==0&&(part_pid==13||part_pid==-13)");
muons = (TH1F*)countCanv.GetPrimitive("muons");
Int_t lostmuons = muons->GetEntries();
fulltree->Draw("part_trackid>>muons","part_mrdhit==1&&part_processEnd!=0&&(part_pid==13||part_pid==-13)");
muons = (TH1F*)countCanv.GetPrimitive("muons");
Int_t stoppedmuons = muons->GetEntries();

cout<<"total muons: "<<allmuons<<" of which "<<hitmuons<<" hit the MRD, and "<<stoppedmuons<<" were stopped, while "<<lostmuons<<" were not."<<endl;

}

//------------------------------------------------------------------
//------------------------------------------------------------------
TString ConvertParticleCodeToName(Int_t particleCode){

	// only do this once to fill map
	/*static bool makeparticlecodestonames = true;
	if(makeparticlecodestonames){
	    makeparticlecodestonames=false;
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(100,"opticalphoton"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(-11,"e+"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(11,"e-"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(22,"gamma"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(-13,"mu+"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(13,"mu-"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(111,"pi0"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(211,"pi+"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(-211,"pi-"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(2112,"neutron"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(2212,"proton"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(14,"nu_mu"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(12,"nu_e"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(-12,"anti_nu_e"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(-14,"anti_nu_mu"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3328,"alpha"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3329,"deuteron"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3330,"triton"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3331,"C10"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3332,"C12"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3333,"N15"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3334,"N16"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3335,"O16"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3336,"Gd158"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3337,"Gd156"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3338,"Gd157"));
	    ParticleCodesToNames.insert(std::pair<Int_t,TString>(3339,"Gd155"));
	}
	*/

	std::map<Int_t,TString>::iterator it;
	it = ParticleCodesToNames.find(particleCode);
	if (it!=ParticleCodesToNames.end()){
		return it->second;
	} else {
		Char_t buffer[80];
		sprintf(buffer,"%d",particleCode);
		TString particleName=buffer;
		return particleName;
	}
}

//------------------------------------------------------------------
//------------------------------------------------------------------

TString ConvertProcessCodeToName(Int_t processCode){

	// only do this once to fill map
	/*static bool makeprocesscodestonames = true;
	if(makeprocesscodestonames){
	    makeprocesscodestonames=false;
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(0,"Transportation"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(1,"Scintillation"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(2,"Cerenkov"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(3,"phot"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(4,"OpAbsorption"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(5,"eIoni"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(6,"muIoni"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(7,"hIoni"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(8,"eBrem"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(9,"muBrem"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(10,"HadronElastic"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(11,"nCapture"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(12,"compt"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(13,"Decay"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(14,"muMinusCaptureAtRest"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(15,"NeutronInelastic"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(16,"ProtonInelastic"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(17,"conv"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(18,"annihil"));
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(19,"PiMinusAbsorptionAtRest")); 
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(20,"PionMinusInelastic")); 
	    ProcessCodesToNames.insert(std::pair<Int_t,TString>(21,"PionPlusInelastic")); 
	}
	*/

	std::map<Int_t,TString>::iterator it;
	it = ProcessCodesToNames.find(processCode);
	if (it!=ProcessCodesToNames.end()){
		return it->second;
	} else {
		Char_t buffer[80];
		sprintf(buffer,"%d",processCode);
		TString processName=buffer;
		return processName;
	}
}

//------------------------------------------------------------------
//------------------------------------------------------------------

Int_t ConvertParticleNameToCode(TString particleName){
	std::map<TString,Int_t>::iterator it;
	it = ParticleNamesToCodes.find(particleName);
	if (it!=ParticleNamesToCodes.end()){ return it->second;} else { return -1;}
}

Int_t ConvertProcessNameToCode(TString processName){
	std::map<TString,Int_t>::iterator it;
	it = ProcessNamesToCodes.find(processName);
	if (it!=ProcessNamesToCodes.end()){ return it->second;} else { return -1;}
}

void GenerateMaps(){

	// only do this once to fill map
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("opticalphoton",100));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("e+",-11));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("e-",11));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("gamma",22));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("mu+",-13));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("mu-",13));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("pi0",111));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("pi+",211));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("pi-",-211));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("neutron",2112));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("proton",2212));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("nu_mu",14));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("nu_e",12));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("anti_nu_e",-12));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("anti_nu_mu",-14));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("alpha",3328));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("deuteron",3329));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("triton",3330));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Li7",3351));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("C10",3331));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("B11",3345));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("C12",3332));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("C13",3350));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("N13",3349));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("N14",3340));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("N15",3333));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("N16",3334));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("O16",3335));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Al27",3346));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Fe54",3341));
    	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Mn54",3348));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Mn55",3342));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Mn56",3352));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Fe56",3343));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Fe57",3344));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Fe58",3347));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Eu154",3353));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Gd158",3336));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Gd156",3337));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Gd157",3338));
	    ParticleNamesToCodes.insert(std::pair<TString,Int_t>("Gd155",3339));
	
	    for (std::map<TString,Int_t>::iterator it=ParticleNamesToCodes.begin(); it!=ParticleNamesToCodes.end(); ++it){
	    	ParticleCodesToNames.insert(std::pair<Int_t,TString>(it->second,it->first));
	    }
	
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("Transportation",0));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("Scintillation",1));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("Cerenkov",2));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("phot",3));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("OpAbsorption",4));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("eIoni",5));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("muIoni",6));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("hIoni",7));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("eBrem",8));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("muBrem",9));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("HadronElastic",10));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("nCapture",11));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("compt",12));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("Decay",13));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("muMinusCaptureAtRest",14));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("NeutronInelastic",15));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("ProtonInelastic",16));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("conv",17));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("annihil",18));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("PiMinusAbsorptionAtRest",19));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("PionMinusInelastic",20));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("PionPlusInelastic",21));
    	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("Electromagnetic",22));
	    ProcessNamesToCodes.insert(std::pair<TString,Int_t>("MultipleScattering",23));
	    
	    for (std::map<TString,Int_t>::iterator it=ProcessNamesToCodes.begin(); it!=ProcessNamesToCodes.end(); ++it){
	    	ProcessCodesToNames.insert(std::pair<Int_t,TString>(it->second,it->first));
	    }
}
/*
std::vector<ROOT::Math::XYZTVector> * pTracks = &tracks; 
// note ROOT::Math::XYZTVector is typedef of ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
// or in namespace ROOT::Math, XYZTVector is a LorentzVector<PxPyPzE4D<double> >
// Here PxPyPzE4D<double> is a 4-vector CLASS that stores all sorts of useful vector functions such as Mag, Perp, Phi, Et etc. 
// The same class is used for both position and time or momentum and energy. (There is no XYZT<double>.)
// see https://root.cern.ch/root/html/ROOT__Math__PxPyPzE4D_double_.html
// LorentzVector appears to be effectively the same, representing one position/momentum 4-vector, but provides additional
// members and methods such as Beta, Gamma, IsTimelike and BoostToCM. Maybe not necessary, but could be useful in the future.
t1.Branch("tracks","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pTracks);
*/

// for finding the PMT num for a given Strike. Should filter out noise by finding mode. 
// note this implementation scans over the array range (not array elements) so will be more efficient for arrays with many
// elements within short range (6,6,11,5,8,9) but inefficient for sparse large range (1000, 2000, 3000)
Int_t findmode(Int_t* arrayin, Int_t arraysize){
    Int_t largestNum = 0, largestCount = 0;
    Int_t* upperLimit = std::max_element(&arrayin[0], &arrayin[arraysize]);
    Int_t* lowerLimit = std::min_element(&arrayin[0], &arrayin[arraysize]);
    for (Int_t i=*lowerLimit; i<=*upperLimit; ++i)    {
        Int_t currentCount = std::count(&arrayin[0], &arrayin[arraysize], i);
        if (currentCount > largestCount)        {
            largestCount = currentCount;
            largestNum = i;
        }
    }
    return largestNum;
}
