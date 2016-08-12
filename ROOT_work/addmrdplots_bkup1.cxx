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

TString ConvertProcessCodeToName(Int_t processCode);
TString ConvertParticleCodeToName(Int_t particleCode);
Int_t ConvertParticleNameToCode(TString particleName);
Int_t ConvertProcessNameToCode(TString particleName);	
typedef std::vector< std::pair<TString,Int_t> > stringintpairvec;
void GenerateVectors(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void GeneratePieCharts(stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void MakeEnergyDepPlots(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids);
void PlotPenetration(TTree* mrdtree);
void EnergyDepositioninLayers(TTree* mrdtree);
void EnergyDepositionInPaddles(TTree* mrdtree);
void countmrdhits(TTree* fulltree);
void splitindividualsources(TTree* mrdtree);

Double_t mrdzstart=1990;
Double_t mrdzlen=919;
Int_t mrdnumlayers=12;		// 12 layers of scintillator, steel. 13 alu structs, but they're not really important. 
Int_t mrdpaddlesperpanel=30;	// remember TWO ROWS of 15 paddles!!
Double_t maxtrackduration=30;	// used to split tracks in MRD 

//------------------------------------------------------------------
//------------------------------------------------------------------


void addmrdplots(const Char_t* mrdinfile="MRDEvents2.root",const Char_t* fullinfile="FullEvent2.root"){
// will add standard set of graphs to MRDEvents.root file
//*******************************************************
// Open mrd hit file and check valid
TFile mrdfile = TFile(mrdinfile);
if (mrdfile.IsZombie()) {
	std::cout << "Error opening MRD file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* mrdtree = (TTree*)mrdfile.Get("MRDTree"); // not strictly necessary - can access MRDTree directly, at least in CINT.

// Open full event file and check valid
TFile fullfile = TFile(fullinfile);
if (mrdfile.IsZombie()) {
	std::cout << "Error opening FullEvent file" << std::endl;
	return;//exit(-1);
}
// Open MRD events tree
TTree* fulltree = (TTree*)fullfile.Get("EventTree");

// Open output file
TFile* rootfileout = new TFile("MRDEventsOut.root", "RECREATE");

stringintpairvec mrdparticleids;
stringintpairvec mrdprocessids;

GenerateVectors(mrdtree,mrdparticleids,mrdprocessids);

GeneratePieCharts(mrdparticleids, mrdprocessids);

MakeEnergyDepPlots(mrdtree, mrdparticleids, mrdprocessids);

PlotPenetration(mrdtree);


EnergyDepositionInPaddles(mrdtree);
EnergyDepositioninLayers(mrdtree);

splitindividualsources(mrdtree);

countmrdhits(fulltree);

//rootfileout->Write(); //Write all objects to disk. 
rootfileout->Close();
delete rootfileout;
mrdfile.Close();
fullfile.Close();
//delete mrdfile;	//close the file (only for pointers)
cout<<"marco"<<endl;
}


void GenerateVectors(TTree* mrdtree, stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids){
TCanvas vectCanv = TCanvas("vectCanv","Title");
vectCanv.cd();	// ensure selected
// Find the particle IDs of particles that have interacted (probably e's, mu's and gammas):
mrdtree->Draw("mrdhit_particleID>>histpid");	
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
		typeandnumhits = std::make_pair(ConvertProcessCodeToName(histproc->GetXaxis()->GetBinLowEdge(i)),binconts);
		mrdprocessids.push_back(typeandnumhits);
	}
}

}

void GeneratePieCharts(stringintpairvec &mrdparticleids, stringintpairvec &mrdprocessids){
// this function generates pie charts of the particles and processes that have interacted/occurred in the MRD. 

// Produce pie chart of particles that interacted in the MRD, and their frequency:
// ============================================================================
//TCanvas* canvloc1 = new TCanvas("canvloc1","Pie Charts",385,110,700,867); //for pie charts
//canvloc1->Divide(1,2);
//TPie* mrdparticlepielocal;
//TPie* mrdprocesspielocal;
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
mrdparticlepie.SetAngularOffset(30.);
mrdparticlepie.SetLabelFormat("#splitline{%txt}{%val (%perc)}");
mrdparticlepie.SetValueFormat("%4.0f");
mrdparticlepie.SetPercentFormat("%3.0f");
mrdparticlepie.SetCircle(0.5, 0.4702026, 0.3302274);
mrdparticlepie.Write();
pieCanv.cd(1);
mrdparticlepie.Draw("3d");
//** local print for check
/*mrdparticlepielocal = new TPie(mrdparticlepie);
mrdparticlepielocal->SetAngularOffset(30.);
mrdparticlepielocal->SetLabelFormat("#splitline{%txt}{%val (%perc)}");
mrdparticlepielocal->SetValueFormat("%4.0f");
mrdparticlepielocal->SetPercentFormat("%3.0f");
mrdparticlepielocal->SetCircle(0.5, 0.4702026, 0.3302274);
canvloc1->cd(1);
//mrdparticlepielocal->Draw("3d");
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
mrdprocesspie->SetLabelFormat("#splitline{%txt}{%val (%perc)}");
mrdprocesspie->SetValueFormat("%4.0f");
mrdprocesspie->SetPercentFormat("%3.0f");
mrdprocesspie->SetCircle(0.5, 0.4702026, 0.3302274);	// reduce size so labels actually fit on the canvas ¬_¬
mrdprocesspie->Write();
pieCanv.cd(2);
mrdprocesspie->Draw("3d");
pieCanv.SaveAs("piecharts.png");
//** local print for check
/*mrdprocesspielocal = new TPie(*mrdprocesspie);
mrdprocesspielocal->SetLabelFormat("#splitline{%txt}{%val (%perc)}");
mrdprocesspielocal->SetPercentFormat("%3.0f");
mrdprocesspielocal->SetValueFormat("%4.0f");
canvloc1->cd(2);
//mrdprocesspielocal->Draw("3d");
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
alledeposits->Write();

//** local print for check
//Note: legend is now generated with file resident histograms - uncommenting to show local may not display the legend as it takes pointers to objects that will go out of scope when the macro ends. If necessary change legend line in loop to use histoelocal instead of newhisto and use the below code for the sum entry. This legend will work for file too, but sum will need to be generated before save -> i.e. move this 'local' code before the paragraph above. 
/*TCanvas* canvloc2 = new TCanvas("canvloc2","Energy Deposition Chart",79,59,1173,712); // for energy deposition vs length
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
alegend->Draw();*/

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

void splitindividualsources(TTree* mrdtree){
// Split hits from each event into groups of hits corresponding to a single 'track' - in case we have more than one track per event
// ================================================================================================================================
TCanvas splittracksCanv = TCanvas("splittracksCanv","Title");
Int_t numentries = mrdtree->GetEntries();
Int_t numhitsthisevent=0;
Int_t numhitsthiseventbranch=0;
TBranch* numhitsbranch = mrdtree->GetBranch("hitnum");
numhitsbranch->SetAddress(&numhitsthiseventbranch);
TBranch* hittimebranch = mrdtree->GetBranch("mrdhit_t");
Double_t* hittimes; 

TH1F* splittrackshist;
TTree* tracksplittree = new TTree("tracksplittree","Tree for grouping hits in an event into tracks");

Int_t numtracksthisevent;
mrdtree->Draw("Entry$>>hitnumhist","hitnum");
TH1F* hitnumhist = (TH1F*)splittracksCanv.GetPrimitive("hitnumhist");
Int_t maxtracksperevent = hitnumhist->GetMaximum();
Int_t* tracknumthisevent = new Int_t[maxtracksperevent];
tracksplittree->Branch("hitnum",&numtracksthisevent);
tracksplittree->Branch("tracknumthisevent",tracknumthisevent,"tracknumthisevent[hitnum]/D");

// to be made a friend tree of mrdtree
cout<<"scanning over "<<numentries<<" entries."<<endl;
for(Int_t entry=0;entry<numentries;entry++){
	numhitsbranch->GetEntry(entry); // first check if we have any hits this entry. 
	numhitsthisevent = numhitsthiseventbranch;	// numhitsthiseventbranch WILL GET OVERWRITTEN BY THE DRAW() CALL BELOW!
	if(numhitsthisevent==0){
		numtracksthisevent=0;
		tracknumthisevent[0]=0;
		tracksplittree->Fill();	// fill the tree anyway so entries align
		continue;
	}	// skip remainder
	cout<<"entry "<<entry<<" has "<<numhitsthisevent<<" hits."<<endl;
	// if continued to here, we have at least some hits. 
//	hittimes = new Double_t[numhitsthisevent];	// necessary to allocate memory now for draw.
//	hittimebranch->SetAddress(hittimes);
	hittimebranch->ResetAddress();
	mrdtree->Draw("mrdhit_t>>splittrackshist(100)",TString::Format("Entry$==%d",entry));
	splittrackshist = (TH1F*)splittracksCanv.GetPrimitive("splittrackshist");
	// first check: are all hits within a 30ms window? If so, just one track. 
	if((splittrackshist->GetXaxis()->GetXmax() - splittrackshist->GetXaxis()->GetXmin())<maxtrackduration){
		cout<<"all hits this event within one track."<<endl;
		// all hits are within a short period of time - assume they are from the same track.
		numtracksthisevent=1;
		for(Int_t j=0;j<numhitsthisevent;j++){tracknumthisevent[j]=0;}
		tracksplittree->Fill();
	} else {
		cout<<"multiple tracks this event"<<endl;
		// this event has multiple tracks. Need to split hits into which track they belong to.
		// first, to keep track of what's been assigned, initialize all to -1. 
		for(Int_t j=0;j<numhitsthisevent;j++){tracknumthisevent[j]=-1;}
		// split tracks by plotting times as a histogram and finding peaks.
		/* simple way, but cannot find the number of peaks: pm is pointer to an array, but cannot be used to retrieve size 
		splittrackshist->ShowPeaks(.0005,"nobackground new nodraw",0.1)	
		// nobackground seems to cure error "Too large clipping window" and seems to improve actual peak finding,
		// including locating peak at 0. 'new' not sure about - it will replace histogram unless 'nodraw' option is also used. 
		TPolyMarker *pm = (TPolyMarker*)(splittrackshist->GetListOfFunctions()->FindObject("TPolyMarker");
		*/
		hittimes = new Double_t[numhitsthisevent];
		hittimebranch->SetAddress(hittimes);
		
		TSpectrum *sp1 = new TSpectrum();
		sp1->Search(splittrackshist,0.,"nobackground",0.01);	// Search(TH1*,width, option, threshold)
		numtracksthisevent = sp1->GetNPeaks();
    		Float_t *trackhittimes = sp1->GetPositionX();
    		cout<<numtracksthisevent<<" tracks this event"<<endl;
		hittimebranch->GetEntry(entry);
		for(Int_t j=0; j<numtracksthisevent; j++){
			cout<<"Track struck MRD at = "<<trackhittimes[j]<<"ms in event "<<entry<<endl;
			Int_t binnumber = splittrackshist->FindBin(trackhittimes[j]);
			while(splittrackshist->GetBinContent(binnumber+1)!=0){
				binnumber++;
			}	// Search will find peak in hit time, but we need to encapsulate all hits, so scan forward from
				// the peak until we have an empty bin. 
			Double_t trackmaxtime = splittrackshist->GetXaxis()->GetBinUpEdge(binnumber);
			cout<<"Setting track time threshold at "<<trackmaxtime<<endl;
			for(Int_t thishit=0;thishit<numhitsthisevent;thishit++){
				if(tracknumthisevent[thishit]<0&&hittimes[thishit]<trackmaxtime){tracknumthisevent[thishit]=j;}
			}
		}
		for(Int_t k=0;k<numhitsthisevent;k++){if(tracknumthisevent[k]==-1){cout<<"*****unbinned hit!"<<k<<" "<<hittimes[k]<<endl;}}
		tracksplittree->Fill();
		delete sp1;
		delete[] hittimes;
	}	// next hit in entry.
}	// next entry in mrdtree

//tracksplittree->Write();
mrdtree->AddFriend("tracksplittree");

delete[] tracknumthisevent;
mrdtree->ResetBranchAddresses();

}		

/*void DoReconstruction(TTree* mrdtree){
	
	TTree* recotree = new TTree
	recotree->Branch(
	
	
	recotree->Fill();
}	*/	

		


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
    TString particleName="";
    
    if(particleCode==100) particleName="opticalphoton";
    if(particleCode==-11) particleName="e+";
    if(particleCode==11) particleName="e-";
    if(particleCode==22) particleName="gamma";
    if(particleCode==-13) particleName="mu+";
    if(particleCode==13) particleName="mu-";
    if(particleCode==111) particleName="pi0";
    if(particleCode==211) particleName="pi+";
    if(particleCode==-211) particleName="pi-";
    if(particleCode==2112) particleName ="neutron";
    if(particleCode==2212) particleName = "proton";
    if(particleCode==14) particleName="nu_mu";
    if(particleCode==12) particleName="nu_e";
    if(particleCode==-12) particleName="anti_nu_e";
    if(particleCode==-14) particleName="anti_nu_mu";
    
    if(particleName==""){Char_t buffer[80];
        		sprintf(buffer,"%d",particleCode);
			particleName=buffer;}
    
    return particleName;
}
//------------------------------------------------------------------
//------------------------------------------------------------------
TString ConvertProcessCodeToName(Int_t processCode){
    TString processName="";

    if(processCode==0) processName="Transportation";
    if(processCode==1) processName="Scintillation";
    if(processCode==2) processName="Cerenkov";
    if(processCode==3) processName="phot";
    if(processCode==4) processName="OpAbsorption";
    if(processCode==5) processName="eIoni";
    if(processCode==6) processName="muIoni";
    if(processCode==7) processName="hIoni";
    if(processCode==8) processName="eBrem";
    if(processCode==9) processName="muBrem";
    if(processCode==10) processName ="HadronElastic";
    if(processCode==11) processName = "nCapture";
    if(processCode==12) processName="compt";
    if(processCode==13) processName="Decay";
    if(processCode==14) processName="muMinusCaptureAtRest";
    if(processCode==15) processName="NeutronInelastic";
    if(processCode==16) processName="ProtonInelastic";
    if(processCode==17) processName="conv";
    if(processCode==18) processName="annihil";
    if(processCode==19) processName="PiMinusAbsorptionAtRest"; 
    if(processCode==20) processName="PionMinusInelastic"; 
    if(processCode==21) processName="PionPlusInelastic"; 
    
    if(processName==""){Char_t buffer[80];
    			sprintf(buffer,"%d",processCode);
    			processName=buffer;}
   
    return processName;
}
//------------------------------------------------------------------
//------------------------------------------------------------------


Int_t ConvertParticleNameToCode(TString particleName){
    Int_t particleCode = -1;
    if(particleName=="opticalphoton") particleCode=100;
    if(particleName=="e+") particleCode=-11;
    if(particleName=="e-") particleCode=11;
    if(particleName=="gamma") particleCode=22;
    if(particleName=="mu+") particleCode=-13;
    if(particleName=="mu-") particleCode=13;
    if(particleName=="pi0") particleCode=111;
    if(particleName=="pi+") particleCode=211;
    if(particleName=="pi-") particleCode=-211;
    if(particleName=="neutron") particleCode = 2112 ;
    if(particleName=="proton") particleCode = 2212 ;
    if(particleName=="nu_mu") particleCode=14;
    if(particleName=="nu_e") particleCode=12;
    if(particleName=="anti_nu_e") particleCode=-12;
    if(particleName=="anti_nu_mu") particleCode=-14;
    if(particleName=="alpha") particleCode=3328;
    if(particleName=="deuteron") particleCode=3329;
    if(particleName=="triton") particleCode=3330;
    if(particleName=="C10[0.0]") particleCode=3331;
    if(particleName=="C12[0.0]") particleCode=3332;
    if(particleName=="N15[0.0]") particleCode=3333;
    if(particleName=="N16[0.0]") particleCode=3334;
    if(particleName=="O16[0.0]") particleCode=3335;
    if(particleName=="Gd158[0.0]") particleCode=3336;
    if(particleName=="Gd156[0.0]") particleCode=3337;
    if(particleName=="Gd157[0.0]") particleCode=3338;
    if(particleName=="Gd155[0.0]") particleCode=3339;
    
    return particleCode;
}

Int_t ConvertProcessNameToCode(TString processName){
    Int_t processCode=-1;
    if(processName=="Transportation") processCode=0;
    if(processName=="Scintillation") processCode=1;
    if(processName=="Cerenkov") processCode=2;
    if(processName=="phot") processCode=3;
    if(processName=="OpAbsorption") processCode=4;
    if(processName=="eIoni") processCode=5;
    if(processName=="muIoni") processCode=6;
    if(processName=="hIoni") processCode=7;
    if(processName=="eBrem") processCode=8;
    if(processName=="muBrem") processCode=9;
    if(processName=="HadronElastic") processCode=10;
    if(processName=="nCapture") processCode=11;
    if(processName=="compt") processCode=12;
    if(processName=="Decay") processCode=13;
    if(processName=="muMinusCaptureAtRest") processCode=14;
    if(processName=="NeutronInelastic") processCode=15;
    if(processName=="ProtonInelastic") processCode=16;
    if(processName=="conv") processCode=17;
    if(processName=="annihil") processCode=18;
    if(processName=="PiMinusAbsorptionAtRest") processCode=19;
    if(processName=="PionMinusInelastic") processCode=20;
    if(processName=="PionPlusInelastic") processCode=21;
    
    return processCode;
}

