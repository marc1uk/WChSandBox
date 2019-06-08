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

using namespace std;

void coincidence_search(const Char_t* infile="/home/marc/LinuxSystemFiles/WChSandBox/build/MRDEvents.root"){

cout<<"opening file"<<endl;
TFile* mrdfile = TFile::Open(infile, "READ");
if(mrdfile==0||mrdfile->IsZombie()){cout<<"input file is a zombie arrrrgh!"<<endl; return;}
cout<<"getting trees"<<endl;
TTree* mrdtree = (TTree*)mrdfile->Get("MRDTree");
TTree* mrdpmttree = (TTree*)mrdfile->Get("MRDPMTTree");
TTree* facctree = (TTree*)mrdfile->Get("FACCTree");
TTree* faccpmttree = (TTree*)mrdfile->Get("FACCPMTTree");
if(mrdtree==0||mrdpmttree==0||facctree==0||faccpmttree==0){cout<<"trees are zombies aargh!"<<endl; return;}

TFile* mrdanalysis = new TFile("mrdanalysis","RECREATE");
TFile* vetoanalysis = new TFile("vetoanalysis","RECREATE");

cout<<"getting entries"<<endl;
Int_t mrdhits = mrdtree->GetEntriesFast();
Int_t mrdpmthits = mrdpmttree->GetEntriesFast();
Int_t facchits = facctree->GetEntriesFast();
Int_t faccpmthits = faccpmttree->GetEntriesFast();
cout<<"entries: "<<mrdhits<<", "<<mrdpmthits<<", "<<facchits<<", "<<faccpmthits<<endl;

cout<<"setting branch addresses"<<endl;
TBranch* nummrdhitsthiseventb;
Int_t nummrdhitsthisevent;
mrdtree->SetBranchAddress("mrd_numhits",&nummrdhitsthisevent,&nummrdhitsthiseventb);
TBranch* numfacchitsthiseventb;
Int_t numfacchitsthisevent;
facctree->SetBranchAddress("facc_numhits",&numfacchitsthisevent,&numfacchitsthiseventb);
TBranch* nummrdpmthitsthiseventb;
Int_t nummrdpmthitsthisevent;
mrdpmttree->SetBranchAddress("mrdpmt_numhits",&nummrdpmthitsthisevent,&nummrdpmthitsthiseventb);
TBranch* numfaccpmthitsthiseventb;
Int_t numfaccpmthitsthisevent;
faccpmttree->SetBranchAddress("faccpmt_numhits",&numfaccpmthitsthisevent,&numfaccpmthitsthiseventb);

TBranch* mrdtrackbranch = mrdtree->GetBranch("mrdhit_trackID");
Int_t* mrdtrackids=0;
Int_t mrdtrackid;
TBranch* mrdpmttrackbranch = mrdpmttree->GetBranch("mrdpmthit_parentID");
Int_t* mrdpmttrackids=0;
Int_t mrdpmttrackid;
TBranch* facctrackbranch = facctree->GetBranch("facchit_trackID");
Int_t* facctrackids=0;
Int_t facctrackid;
TBranch* faccpmttrackbranch = faccpmttree->GetBranch("faccpmthit_parentID");
Int_t* faccpmttrackids=0;
Int_t faccpmttrackid;

if(nummrdhitsthiseventb==0||numfacchitsthiseventb==0||nummrdpmthitsthiseventb==0||numfaccpmthitsthiseventb==0||mrdtrackbranch==0||mrdpmttrackbranch==0||facctrackbranch==0||faccpmttrackbranch==0){cout<<"branches are zombies aargh!"<<endl; return;}

// multimap of events vs track ids that had hits
std::multimap<int,int> mrdhittracks, facchittracks, mrdpmthittracks, faccpmthittracks;
std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> tracksinthisevent, tracksinthisevent2;
bool tracknoted;
int numhits;

// &&&&&&&&&&&&&&&&&&&&&&
/* 
THINK ABOUT WHAT YOU WANT OUT
*/
// &&&&&&&&&&&&&&&&&&&&&&

cout<<"getting mrd hits"<<endl;
numhits=0;
for(int entry=0;entry<mrdhits;entry++){
	//cout<<"getting num hits in entry "<<entry<<endl;
	nummrdhitsthiseventb->GetEntry(entry);
	//cout<<nummrdhitsthisevent<<" hits this entry"<<endl;
	//cout<<"resetting hit track ids branch address"<<endl;
	mrdtrackbranch->ResetAddress();
	//cout<<"cleaning up last array of trackids"<<endl;
	if(mrdtrackids){delete mrdtrackids;}
	//cout<<"creating new array of trackids"<<endl;
	mrdtrackids = new Int_t[nummrdhitsthisevent];
	//cout<<"setting branch address"<<endl;
	mrdtrackbranch->SetAddress(mrdtrackids);
	//cout<<"getting track ids for this entry"<<endl;
	mrdtrackbranch->GetEntry(entry);
	//cout<<"trackIDs got"<<endl;
	for(int innerentry=0;innerentry<nummrdhitsthisevent;innerentry++){
		//cout<<"getting track id "<<innerentry<<" of this entry"<<endl;
		mrdtrackid=mrdtrackids[innerentry];
		//cout<<"checking if we have this trackid recorded for this entry"<<endl;
		tracksinthisevent = mrdhittracks.equal_range(entry);
		tracknoted=false;
		for(std::multimap<int,int>::iterator it=tracksinthisevent.first;it!=tracksinthisevent.second;it++){
			if((*it).second==mrdtrackid){
				tracknoted=true;
				break;
			}
		}
		if(!tracknoted){
			mrdhittracks.insert(std::pair<int,int>(entry,mrdtrackid));
		}
	}
}
cout<<"found "<<mrdhittracks.size()<<" tracks with hits in the MRD"<<endl;

for(int entry=0;entry<mrdpmthits;entry++){
	nummrdpmthitsthiseventb->GetEntry(entry);
	mrdpmttrackbranch->ResetAddress();
	if(mrdpmttrackids){delete mrdpmttrackids;}
	mrdpmttrackids = new Int_t[nummrdpmthitsthisevent];
	mrdpmttrackbranch->SetAddress(mrdpmttrackids);
	mrdpmttrackbranch->GetEntry(entry);
	//if(nummrdpmthitsthisevent!=0){cout<<nummrdpmthitsthisevent<<" mrd pmt hits this event"<<endl;}
	for(int innerentry=0;innerentry<nummrdpmthitsthisevent;innerentry++){
		mrdpmttrackid=mrdpmttrackids[innerentry];
		//cout<<"track "<<mrdpmttrackid;
		tracksinthisevent = mrdpmthittracks.equal_range(entry);
		tracknoted=false;
		for(std::multimap<int,int>::iterator it=tracksinthisevent.first;it!=tracksinthisevent.second;it++){
			if((*it).second==mrdpmttrackid){
				tracknoted=true;
				//cout<<" found"<<endl;
				break;
			}
		}
		if(!tracknoted){
			//cout<<" not found"<<endl;
			mrdpmthittracks.insert(std::pair<int,int>(entry,mrdpmttrackid));
		}
	}
}
cout<<"found "<<mrdpmthittracks.size()<<" tracks with detected hits"<<endl;

cout<<"getting facc hits"<<endl;
for(int entry=0;entry<facchits;entry++){
	numfacchitsthiseventb->GetEntry(entry);
	facctrackbranch->ResetAddress();
	if(facctrackids){delete facctrackids;}
	facctrackids = new Int_t[numfacchitsthisevent];
	facctrackbranch->SetAddress(facctrackids);
	facctrackbranch->GetEntry(entry);
	for(int innerentry=0;innerentry<numfacchitsthisevent;innerentry++){
		facctrackid=facctrackids[innerentry];
		tracksinthisevent = facchittracks.equal_range(entry);
		tracknoted=false;
		for(std::multimap<int,int>::iterator it=tracksinthisevent.first;it!=tracksinthisevent.second;it++){
			if((*it).second==facctrackid){
				tracknoted=true;
				break;
			}
		}
		if(!tracknoted){
			facchittracks.insert(std::pair<int,int>(entry,facctrackid));
		}
	}
}
cout<<"found "<<facchittracks.size()<<" tracks with hits on the FACC"<<endl;

for(int entry=0;entry<faccpmthits;entry++){
	numfaccpmthitsthiseventb->GetEntry(entry);
	faccpmttrackbranch->ResetAddress();
	if(faccpmttrackids){delete faccpmttrackids;}
	faccpmttrackids = new Int_t[numfaccpmthitsthisevent];
	faccpmttrackbranch->SetAddress(faccpmttrackids);
	faccpmttrackbranch->GetEntry(entry);
	//if(numfaccpmthitsthisevent!=0){cout<<numfaccpmthitsthisevent<<" facc pmt hits this event"<<endl;}
	for(int innerentry=0;innerentry<numfaccpmthitsthisevent;innerentry++){
		faccpmttrackid=faccpmttrackids[innerentry];
		//cout<<"track "<<faccpmttrackid;
		tracksinthisevent = faccpmthittracks.equal_range(entry);
		tracknoted=false;
		for(std::multimap<int,int>::iterator it=tracksinthisevent.first;it!=tracksinthisevent.second;it++){
			if((*it).second==faccpmttrackid){
				tracknoted=true;
				//cout<<" found"<<endl;
				break;
			}
		}
		if(!tracknoted){
			//cout<<"track "<<faccpmttrackid<<" added"<<endl;
			faccpmthittracks.insert(std::pair<int,int>(entry,faccpmttrackid));
		}
	}
}
cout<<"found "<<faccpmthittracks.size()<<" tracks with detected hits"<<endl;

cout<<mrdhittracks.size()<<" tracks hit the MRD, of which "<<mrdpmthittracks.size()<<" left hits on PMTs"<<endl;
cout<<facchittracks.size()<<" tracks hit the FACC, of which "<<faccpmthittracks.size()<<" left hits on PMTs"<<endl;

std::map<int,std::vector<int>> tracksinbothits;
cout<<"looking for coincidences in hits"<<endl;
// loop over events
for(std::multimap<int,int>::iterator eventit=mrdhittracks.begin();eventit!=mrdhittracks.end();++eventit){
	Int_t eventnum = (*eventit).first;
	// see if there are any tracks from this event with hits also in the facc
	if(facchittracks.count(eventnum)>0){
		// look through tracks in each of the detectors within this event and see if any of the tracks match
		tracksinthisevent = mrdhittracks.equal_range(eventnum);
		tracksinthisevent2 = facchittracks.equal_range(eventnum);
		// loop over mrd tracks
		for(std::multimap<int,int>::iterator mrdtrackit=tracksinthisevent.first;mrdtrackit!=tracksinthisevent.second;mrdtrackit++){
			int themrdtracknumber = (*mrdtrackit).second;
			// and loop over facc tracks for a match
			for(std::multimap<int,int>::iterator facctrackit=tracksinthisevent2.first;facctrackit!=tracksinthisevent2.second;facctrackit++){
				int thefacctracknumber = (*facctrackit).second;
				if(themrdtracknumber==thefacctracknumber){
					if((tracksinbothits.count(eventnum))>0){
						std::vector<int> tempvec = tracksinbothits.at(eventnum);
						tempvec.push_back(themrdtracknumber);
						tracksinbothits.at(eventnum)=tempvec;
					} else {
						std::vector<int> tempvec (eventnum);
						tracksinbothits.insert(std::pair<int,std::vector<int>>(eventnum,tempvec));
					}
				}
			}
		}
	}
}

std::map<int,std::vector<int>> tracksinbothpmts;
cout<<"looking for coincidences in pmt hits"<<endl;
// loop over events
for(std::multimap<int,int>::iterator eventit=mrdpmthittracks.begin();eventit!=mrdpmthittracks.end();++eventit){
	Int_t eventnum = (*eventit).first;
	// see if there are any tracks from this event with hits also in the facc
	if(faccpmthittracks.count(eventnum)>0){
		// look through tracks in each of the detectors within this event and see if any of the tracks match
		tracksinthisevent = mrdpmthittracks.equal_range(eventnum);
		tracksinthisevent2 = faccpmthittracks.equal_range(eventnum);
		// loop over mrd tracks
		for(std::multimap<int,int>::iterator mrdtrackit=tracksinthisevent.first;mrdtrackit!=tracksinthisevent.second;mrdtrackit++){
			int themrdtracknumber = (*mrdtrackit).second;
			// and loop over facc tracks for a match
			for(std::multimap<int,int>::iterator facctrackit=tracksinthisevent2.first;facctrackit!=tracksinthisevent2.second;facctrackit++){
				int thefacctracknumber = (*facctrackit).second;
				if(themrdtracknumber==thefacctracknumber){
					if((tracksinbothpmts.count(eventnum))>0){
						std::vector<int> tempvec = tracksinbothpmts.at(eventnum);
						tempvec.push_back(themrdtracknumber);
						tracksinbothpmts.at(eventnum)=tempvec;
					} else {
						std::vector<int> tempvec (eventnum);
						tracksinbothpmts.insert(std::pair<int,std::vector<int>>(eventnum,tempvec));
					}
				}
			}
		}
	}
}

cout<<tracksinbothits.size()<<" tracks left hits in both SDs"<<endl;
cout<<tracksinbothpmts.size()<<" tracks left hits in both SD's PMTs"<<endl;

cout<<"resetting trees"<<endl;
mrdtree->ResetBranchAddresses();
facctree->ResetBranchAddresses();
mrdpmttree->ResetBranchAddresses();
faccpmttree->ResetBranchAddresses();
cout<<"closing the file"<<endl;
mrdfile->Close();
cout<<"deleting the file"<<endl;
delete mrdfile;
return;

}


