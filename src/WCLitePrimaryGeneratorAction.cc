#include "WCLitePrimaryGeneratorAction.hh"
#include "WCLiteDetectorConstruction.hh"
#include "WCLitePrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TRandom3.h"
#include "G4OpticalPhoton.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "G4SystemOfUnits.hh"

#include <iostream>
#include "TChain.h"
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using std::vector;
using std::string;
using std::fstream;
using namespace std;

enum EEntryStatus {
       kEntryValid = 0, ///< data read okay
       kEntryNotLoaded, ///< no entry has been loaded yet
       kEntryNoTree, ///< the tree does not exist
       kEntryNotFound, ///< the tree entry number does not exist
       kEntryChainSetupError, ///< problem in accessing a chain element, e.g. file without the tree
       kEntryChainFileError, ///< problem in opening a chain's file
       kEntryDictionaryError, ///< problem reading dictionary info from tree
       kEntryLast, ///< last entry was reached
};

enum  ESetupStatus {
  kSetupNotSetup = -7, kSetupTreeDestructed = -8, kSetupMakeClassModeMismatch = -7, kSetupMissingCounterBranch = -6,
  kSetupMissingBranch = -5, kSetupInternalError = -4, kSetupMissingCompiledCollectionProxy = -3, kSetupMismatch = -2,
  kSetupClassMismatch = -1, kSetupMatch = 0, kSetupMatchBranch = 0, kSetupMatchConversion,
  kSetupMatchConversionCollection, kSetupMakeClass, kSetupVoidPtr, kSetupNoCheck,
  kSetupMatchLeaf
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// function declarations

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  inFile.getline(inBuf,lineSize);
  return tokenize(" $", inBuf);	
}

inline float atof( const string& str ) {return std::atof( str.c_str() );}
inline int   atoi( const string& str ) {return std::atoi( str.c_str() );}

//void LoadFiles(const char* dirname, TTree* tf);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// constructor

WCLitePrimaryGeneratorAction::WCLitePrimaryGeneratorAction(WCLiteDetectorConstruction* myDC) : myDetector(myDC), useNuanceTextFormat(true)
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  if(useLaserEvt){ MyGPS = new G4GeneralParticleSource();}

  if(useNuanceTextFormat){
	  const int kmaxtrk=30;
	  const int kmaxMRDhits=50;
	  mpid = new G4int[kmaxtrk];
	  mpx = new G4double[kmaxtrk];
	  mpy = new G4double[kmaxtrk];
	  mpz = new G4double[kmaxtrk];
	  mKE = new G4double[kmaxtrk];

	  mMRDhitlayer = new G4int[kmaxMRDhits];
	  mMRDhitorientation = new G4int[kmaxMRDhits];
	  mMRDhitEdep = new G4double[kmaxMRDhits];
	  mMRDhitEchdep = new G4double[kmaxMRDhits];
	  mMRDhitx = new G4double[kmaxMRDhits];
	  mMRDhity = new G4double[kmaxMRDhits];
	  mMRDhitz = new G4double[kmaxMRDhits];
	  
	  
	  tcardfile = new TTree("tcardfile","tcardfile");

	  tcardfile->Branch("evt",&mevt); 
	  tcardfile->Branch("mode",&mmode);
	  tcardfile->Branch("neutrino_E",&mE);
	  tcardfile->Branch("neutrino_id",&mbeam_id);
	  tcardfile->Branch("neutrino_px",&mbeam_px);
	  tcardfile->Branch("neutrino_py",&mbeam_py);
	  tcardfile->Branch("neutrino_pz",&mbeam_pz);
	  tcardfile->Branch("ntrks",&ntrks);
	  tcardfile->Branch("nneutrons",&nneutrons);
	  tcardfile->Branch("vtxx",&mvtxx); 
	  tcardfile->Branch("vtxy",&mvtxy); 
	  tcardfile->Branch("vtxz",&mvtxz);

	  tcardfile->Branch("nMRDlayers",&nMRDlayers); 
	  tcardfile->Branch("mMRDtotEdep",&mMRDtotEdep); 
	  tcardfile->Branch("mMRDtotEchdep",&mMRDtotEchdep);

	  tcardfile->Branch("mpid",mpid,"mpid[ntrks]/I");
	  tcardfile->Branch("px",mpx,"px[ntrks]/D"); 
	  tcardfile->Branch("py",mpy,"py[ntrks]/D"); 
	  tcardfile->Branch("pz",mpz,"pz[ntrks]/D");
	  tcardfile->Branch("KE",mKE,"KE[ntrks]/D");

	  tcardfile->Branch("mMRDhitlayer",mMRDhitlayer,"mMRDhitlayer[nMRDlayers]/I");
	  tcardfile->Branch("mMRDhitorientation",mMRDhitorientation,"mMRDhitorientation[nMRDlayers]/I");
	  tcardfile->Branch("mMRDhitEdep",mMRDhitEdep,"mMRDhitEdep[nMRDlayers]/D");
	  tcardfile->Branch("mMRDhitEchdep",mMRDhitEchdep,"mMRDhitEchdep[nMRDlayers]/D");
	  tcardfile->Branch("mMRDhitx",mMRDhitx,"mMRDhitx[nMRDlayers]/D");
	  tcardfile->Branch("mMRDhity",mMRDhity,"mMRDhity[nMRDlayers]/D");
	  tcardfile->Branch("mMRDhitz",mMRDhitz,"mMRDhitz[nMRDlayers]/D");
	//
	////  tcardfile->Branch("issum",&msum); 
	////  tcardfile->Branch("Nneut",&nneut);
	////  tcardfile->Branch("np",&np);
	
	  mR = new TRandom3();
  }


  // Initialize to zero
  mode = 0;
  vtxvol = 0;
  vtx = G4ThreeVector(0.,0.,0.);
  nuEnergy = 0.;
  _counterRock=0; 	// counter for generated in Rock
  _counterCublic=0; 	// counter generated
 
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*m,0.*m,0.*m));
    
  messenger = new WCLitePrimaryGeneratorMessenger(this);
  useMulineEvt = false; //true for old numu_center.txt;
  useOpticalPhotEvt = false;
  useGeantinoEvt = false; 
  useNormalEvt = false;
  useBeamEvt = true;	// true for new geant4 dirt simulations
  
  if(useBeamEvt){
	inputdata = new TChain("tankflux");
	inputdata->Add("fluxesandtables/annie_tank_flux.*.root");
	inputdata->LoadTree(0);
	
	inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
	inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
	vtxxBranch=inputdata->GetBranch("vx");
	vtxyBranch=inputdata->GetBranch("vy");
	vtxzBranch=inputdata->GetBranch("vz");
	vtxtBranch=inputdata->GetBranch("vt");
	pxBranch=inputdata->GetBranch("px");
	pyBranch=inputdata->GetBranch("py");
	pzBranch=inputdata->GetBranch("pz");
	EBranch=inputdata->GetBranch("E");
	KEBranch=inputdata->GetBranch("kE");
	pdgBranch=inputdata->GetBranch("pdgtank");
	
	if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0){
		G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
	}
	
	inputEntry=900;
	runBranch->GetEntry(inputEntry);
	G4cout<<"first run: "<<runbranchval<<G4endl;
	treeNumber=inputdata->GetTreeNumber();
	
	
//	inputfilereader = new TTreeReader(inputdata);
//	Int_t numinputevents = inputfilereader->GetEntries(true);
//	G4cout<<numinputevents<<" events in the input file."<<G4endl;
//	TTreeReader::EEntryStatus entrystats = inputfilereader->SetEntry(0);
//	G4cout<<"entrystat: "<<entrystats<<G4endl;
//	Bool_t nextsuccess = inputfilereader->Next();
//	if(nextsuccess){G4cout<<"success"<<G4endl;} else {G4cout<<"failure"<<G4endl;}
//	
//	TTreeReaderValue<Int_t> rdrrun((*inputfilereader), "run");
//	TTreeReaderValue<Int_t> rdrentry((*inputfilereader), "entry");
//	TTreeReaderValue<Int_t> rdriter((*inputfilereader), "iter");
//	TTreeReaderValue<Int_t> rdrniter((*inputfilereader), "niter");
//	TTreeReaderValue<Int_t> rdrnupdg((*inputfilereader), "nupdg");
//	TTreeReaderValue<Double_t> rdrnuvtxx((*inputfilereader), "nuvtxx");
//	TTreeReaderValue<Double_t> rdrnuvtxy((*inputfilereader), "nuvtxy");
//	TTreeReaderValue<Double_t> rdrnuvtxz((*inputfilereader), "nuvtxz");
//	TTreeReaderValue<Double_t> rdrnuvtxt((*inputfilereader), "nuvtxt");
//	TTreeReaderValue<Int_t> rdrintank((*inputfilereader), "intank");
//	TTreeReaderValue<Int_t> rdrinhall((*inputfilereader), "inhall");
//	TTreeReaderValue<Char_t> rdrvtxvol((*inputfilereader), "vtxvol");
//	TTreeReaderValue<Char_t> rdrvtxmat((*inputfilereader), "vtxmat");
//	TTreeReaderValue<Int_t> rdrntank((*inputfilereader), "ntank");
//	TTreeReaderArray<Int_t> rdrpdgtank((*inputfilereader), "pdgtank");
//	TTreeReaderArray<Int_t> rdrprimary((*inputfilereader), "primary");
//	TTreeReaderArray<Double_t> rdrvx((*inputfilereader), "vx");
//	TTreeReaderArray<Double_t> rdrvy((*inputfilereader), "vy");
//	TTreeReaderArray<Double_t> rdrvz((*inputfilereader), "vz");
//	TTreeReaderArray<Double_t> rdrvt((*inputfilereader), "vt");
//	TTreeReaderArray<Double_t> rdrpx((*inputfilereader), "px");
//	TTreeReaderArray<Double_t> rdrpy((*inputfilereader), "py");
//	TTreeReaderArray<Double_t> rdrpz((*inputfilereader), "pz");
//	TTreeReaderArray<Double_t> rdrE((*inputfilereader), "E");
//	TTreeReaderArray<Double_t> rdrkE((*inputfilereader), "kE");
//	
//	entrystats = inputfilereader->SetEntry(0);
//	G4cout<<"entrystat: "<<entrystats<<G4endl;
//	nextsuccess = inputfilereader->Next();
//	if(nextsuccess){G4cout<<"success"<<G4endl;} else {G4cout<<"failure"<<G4endl;}
//	
//	if(rdrrun.GetSetupStatus()<0) G4cout<<rdrrun.GetSetupStatus()<<G4endl;
//	if(rdrentry.GetSetupStatus()<0) G4cout<<rdrentry.GetSetupStatus()<<G4endl;
//	if(rdriter.GetSetupStatus()<0) G4cout<<rdriter.GetSetupStatus()<<G4endl;
//	if(rdrniter.GetSetupStatus()<0) G4cout<<rdrniter.GetSetupStatus()<<G4endl;
//	if(rdrnupdg.GetSetupStatus()<0) G4cout<<rdrnupdg.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxx.GetSetupStatus()<0) G4cout<<rdrnuvtxx.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxy.GetSetupStatus()<0) G4cout<<rdrnuvtxy.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxz.GetSetupStatus()<0) G4cout<<rdrnuvtxz.GetSetupStatus()<<G4endl;
//	if(rdrnuvtxt.GetSetupStatus()<0) G4cout<<rdrnuvtxt.GetSetupStatus()<<G4endl;
//	if(rdrintank.GetSetupStatus()<0) G4cout<<rdrintank.GetSetupStatus()<<G4endl;
//	if(rdrinhall.GetSetupStatus()<0) G4cout<<rdrinhall.GetSetupStatus()<<G4endl;
//	if(rdrvtxvol.GetSetupStatus()<0) G4cout<<rdrvtxvol.GetSetupStatus()<<G4endl;
//	if(rdrvtxmat.GetSetupStatus()<0) G4cout<<rdrvtxmat.GetSetupStatus()<<G4endl;
//	if(rdrntank.GetSetupStatus()<0) G4cout<<rdrntank.GetSetupStatus()<<G4endl;
//	if(rdrpdgtank.GetSetupStatus()<0) G4cout<<rdrpdgtank.GetSetupStatus()<<G4endl;
//	if(rdrprimary.GetSetupStatus()<0) G4cout<<rdrprimary.GetSetupStatus()<<G4endl;
//	if(rdrvx.GetSetupStatus()<0) G4cout<<rdrvx.GetSetupStatus()<<G4endl;
//	if(rdrvy.GetSetupStatus()<0) G4cout<<rdrvy.GetSetupStatus()<<G4endl;
//	if(rdrvz.GetSetupStatus()<0) G4cout<<rdrvz.GetSetupStatus()<<G4endl;
//	if(rdrvt.GetSetupStatus()<0) G4cout<<rdrvt.GetSetupStatus()<<G4endl;
//	if(rdrpx.GetSetupStatus()<0) G4cout<<rdrpx.GetSetupStatus()<<G4endl;
//	if(rdrpy.GetSetupStatus()<0) G4cout<<rdrpy.GetSetupStatus()<<G4endl;
//	if(rdrpz.GetSetupStatus()<0) G4cout<<rdrpz.GetSetupStatus()<<G4endl;
//	if(rdrE.GetSetupStatus()<0) G4cout<<rdrE.GetSetupStatus()<<G4endl;
//	if(rdrkE.GetSetupStatus()<0) G4cout<<rdrkE.GetSetupStatus()<<G4endl;
//	
//	rdrrunp = &rdrrun;
//	rdrentryp = &rdrentry;
//	rdriterp = &rdriter;
//	rdrniterp = &rdrniter;
//	rdrnupdgp = &rdrnupdg;
//	rdrnuvtxxp = &rdrnuvtxx;
//	rdrnuvtxyp = &rdrnuvtxy;
//	rdrnuvtxtp = &rdrnuvtxt;
//	rdrintankp = &rdrintank;
//	rdrinhallp = &rdrinhall;
//	rdrvtxvolp = &rdrvtxvol;
//	rdrntankp = &rdrntank;
//	rdrpdgtankp = &rdrpdgtank;
//	rdrprimaryp = &rdrprimary;
//	rdrvxp = &rdrvx;
//	rdrvyp = &rdrvy;
//	rdrvzp = &rdrvz;
//	rdrvtp = &rdrvt;
//	rdrpxp = &rdrpx;
//	rdrpyp = &rdrpy;
//	rdrpzp = &rdrpz;
//	rdrEp = &rdrE;
//	rdrkEp = &rdrkE;
	
  }
  
  if(useOpticalPhotEvt){
	  // optical photon gun for testing MRD optical response
	  OpPhotonGun = new G4ParticleGun(G4OpticalPhoton::Definition());
	  OpPhotonGun->SetParticlePosition(G4ThreeVector(3.*cm,3.*cm,(-1.545)*m));
	  // (3*cm,3*cm, 1.569*m) is in first scint panel.
	  // (3*cm,3*cm,-1.545*m) is in the second veto panel.
	  OpPhotonGun->SetParticleMomentumDirection(G4ThreeVector(0.5,-0.5,0.));		// fire in x direction
	  OpPhotonGun->SetParticleEnergy(3*eV);
	  OpPhotonGun->SetParticlePolarization(G4ThreeVector(0.,0.,1.));
  }
  
  if(useGeantinoEvt){
	  // geantino gun for testing geometry
	  GeantinoGun = new G4ParticleGun();
	  //GeantinoGun->SetParticleDefinition( G4Geantino::Geantino() ); 			// default, don't need to set
	  //GeantinoGun->SetParticlePosition(G4ThreeVector(3.0*cm, 3.0*cm, 0.0*cm));
	  //GeantinoGun->SetParticleMomentumDirection( G4ThreeVector(0.0,0.0,1.0) );		// needs to be off-centre due to inter-paddle gaps
	  GeantinoGun->SetParticlePosition(G4ThreeVector(3.0*cm, 3.0*cm,(-1.545)*m));		// in veto paddle
	  GeantinoGun->SetParticleMomentumDirection( G4ThreeVector(0.0,1.0,0.0) );		// fire toward WLS ribbon in +'ve y
	  GeantinoGun->SetParticleEnergy( 1.0*GeV );
	  GeantinoGun->SetParticleTime( 0.0*ns );
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// destructor

WCLitePrimaryGeneratorAction::~WCLitePrimaryGeneratorAction()
{
  if (IsGeneratingVertexInRock()){
    G4cout << "Fraction of Rock volume is : " << G4endl;
      G4cout << " Random number generated in Rock / in Cublic = " 
             << _counterRock << "/" << _counterCublic 
             << " = " << _counterRock/(G4double)_counterCublic << G4endl;
  }
  inputFile.close();
  
  //useMulineEvt = true; //true for old numu_center.txt;
  
  if(useOpticalPhotEvt){delete OpPhotonGun;}
  if(useGeantinoEvt) {delete GeantinoGun;}
  if(useLaserEvt){delete MyGPS;}   //T. Akiri: Delete the GPS variable
  if(useNuanceTextFormat){
	//  TFile *ttttf = new TFile("generatorcardfile.root","RECREATE");
	//  tcardfile->Write();
	//  ttttf->Close();
	//  delete ttttf;
	delete tcardfile;
	// following deleted automatically by closing file/tree?
  	delete[] mpid;
  	delete[] mpx;
  	delete[] mpy;
  	delete[] mpz;
  	delete[] mKE;
  	delete[] mMRDhitlayer;
  	delete[] mMRDhitorientation ;
  	delete[] mMRDhitEdep;
  	delete[] mMRDhitEchdep;
  	delete[] mMRDhitx;
  	delete[] mMRDhity;
  	delete[] mMRDhitz;
  	delete mR;
  }
  if(useBeamEvt){
    inputdata->ResetBranchAddresses();
    delete inputdata;
    delete inputfilereader;
  }
  
  delete particleGun;
  
  delete messenger;
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// generate primaries main method

void WCLitePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Temporary kludge to turn on/off vector text format 

  useNuanceTextFormat = true;
  G4bool fullsim = true;
  G4int vtxchooser=1;

  // in cm
  G4double voffsetx = 0.0;
  G4double voffsety = 0.0;
  G4double voffsetz = -50.0;
  G4double vscalex = 150.00;
  G4double vscaley = 150.00;
  G4double vscalez = 50.00;

  if(vtxchooser==0) G4cout<<"Vertices are forced to 0!!!"<<G4endl;
  if(vtxchooser==1) G4cout<<"Using vertices from cardfile"<<G4endl;
  if(vtxchooser==2) {G4cout<<"Randomized vtx.Offsets:" << voffsetx
			  <<" "<<voffsety<<" "<<voffsetz<<" scales: "
			  <<vscalex<<" "<<vscaley<<" "<<vscalez<<G4endl;}
    
  // Do for every event

  if (useMulineEvt)
  {

    if ( !inputFile.is_open() )
    {
      G4cout << "Set a vector file using the command /mygen/vecfile name"
	     << G4endl;
      return;
    }


    vector< vector<double> > particles;
    vector< vector<double> > layers;

    vector<G4int> idvec;
    vector<G4ThreeVector> dirvec;
    vector<G4double> Evec;
   
    if (useNuanceTextFormat)
      {
	const int lineSize=100;
	char      inBuf[lineSize];
	vector<string> token(1);
	
	token = readInLine(inputFile, lineSize, inBuf);
	  
        if (token.size() == 0) 
	  {
	    G4cout << "end of nuance vector file!" << G4endl;
	  }
	else if (token[0] != "begin")
	  {
	    G4cout << "unexpected line begins with " << token[0] << G4endl;
	  }
	else   // normal parsing begins here
	  {
	    // Read the nuance line (ignore value now)

	    mevt = atoi(token[1]); 

	    token = readInLine(inputFile, lineSize, inBuf);
	    mode = atoi(token[1]);
	    mmode=mode;

	    // Read the info line, basically a dummy
	    // token=readInLine(inputFile, lineSize, inBuf);
	    //G4cout << "Vector File Record Number " << token[2] << G4endl;
            // vecRecNumber = atoi(token[2]);
	    
	    // Now read the outgoing particles
	    // These we will simulate.

	    ntrks=0;
	    nneutrons=0;

	    while ( token=readInLine(inputFile, lineSize, inBuf),
		    ((token[0] == "track")||(token[0] == "vertex")) )
	      {
		// We are only interested in the particles
		// that leave the nucleus, tagged by "0"
		
		if(token[0]=="track"){

		  if( atoi(token[1]) == 14 ){	//muon neutrino
		    //read in neutrino line
		    beampdg = atoi(token[1]);
		    beamenergy = atof(token[2])*MeV;
		    beamdir = G4ThreeVector(atof(token[3]),
					    atof(token[4]),
					    atof(token[5]));

		    mE=beamenergy;
		    mbeam_id = beampdg;
		    mbeam_px = beamdir.x();
		    mbeam_py = beamdir.y();
		    mbeam_pz = beamdir.z();
		    
		    // G4cout<<"neutrino line: "<<token[1]<<" "<<token[2]<<" "<<token[3]<<G4endl;
		  }

		  if( atoi(token[1]) == 8016 ){  // oxygen nucleus
		    //read in target line
		    targetpdg = atoi(token[1]);
		    targetenergy = atof(token[2])*MeV;
		    targetdir = G4ThreeVector(atof(token[3]),
					      atof(token[4]),
					      atof(token[5]));
		    
		    //G4cout<<"target line: "<<token[1]<<" "<<token[2]<<" "<<token[3]<<" "<<token[6]<<" "<<mevt<<G4endl;
		  }


		  //		  if( token[6]=="0" )
		  if( (atoi(token[6]))==0 && (atoi(token[1])!=8016) )	
		  // NUANCE format: The last number on the line is either -1 (initial state particle), -2 (final state particle before interactions) or 0 (the final state particle after interactions). Normally your detector simulation should only track particles with a status code of 0. 
		    {
//		      G4int pdgid = atoi(token[1]);
		      energy = atof(token[2])*MeV; //G4double
		      G4ThreeVector dir = G4ThreeVector(atof(token[3]),
							atof(token[4]),
							atof(token[5]));

		      idvec.push_back( atoi(token[1]) );
		      dirvec.push_back(dir);
		      Evec.push_back(energy);
		      
		      ntrks++;
		    }
		}else {
		  G4cout<<"vertex: "<<token[1]<<" "<<token[2]<<" "<<token[3]<<G4endl;

		  //generate a random vertex



		  //use vertex from cardfile

		  if(vtxchooser==0){
		    //force vertex to 0
		    mvtxx = 0.;    mvtxy = 0.;    mvtxz = 0.;
		  } 

		  if(vtxchooser==1){
		    //vertex from cardfile
		    mvtxx = atof(token[1]);     mvtxy = atof(token[2]);    
		    mvtxz = atof(token[3]);
		  }

		  if(vtxchooser==2){
		    // force random vertex with fractional offset and scale
		    mvtxx = (((mR->Rndm())-0.5)*vscalex + voffsetx)/10.;    
		    mvtxy = (((mR->Rndm())-0.5)*vscaley + voffsety)/10.;   
		    mvtxz = (((mR->Rndm())-0.5)*vscalez + voffsetz)/10.;
		  }

		  //G4cout<<"test1 "<<mvtxx<<" "<<mvtxy<<" "<<mvtxz<<G4endl;

		  vtx = G4ThreeVector(mvtxx*cm,
				      mvtxy*cm,
				      mvtxz*cm);
		  
		  /*
 		  vtx = G4ThreeVector(mvtxx,
				      mvtxy,
				      mvtxz);
		  */

		  //G4cout<<"test2 "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<G4endl;

		}
	      }


	    
//	    G4double npl=-555.;
	    
	    for(int pgi=0; pgi<ntrks; pgi++){
	             
	      mpid[pgi] = idvec.at(pgi);
		    
	      particleGun->
		SetParticleDefinition( particleTable->
				       FindParticle(idvec.at(pgi)) );

	      G4double mass = 
		particleGun->GetParticleDefinition()->GetPDGMass();
	      
	      energy = Evec.at(pgi);

	      G4double ekin = energy - mass;
	      mKE[pgi] = ekin;	     	      		      

	      G4ThreeVector idir = dirvec.at(pgi);
	      
	      mpx[pgi] = idir.x();
	      mpy[pgi] = idir.y();
	      mpz[pgi] = idir.z();
	      
	      if(idvec.at(pgi)==2112){
		nneutrons++; 
	      }


	      if(fullsim){
		// G4cout<<"Final Loop: "<<mpid[pgi]<<" "<<mpx[pgi]<<" "<<mpy[pgi]<<" "<<mpz[pgi]<<" "<<vtx.x()<<" "<<vtx.y()<<" "<<vtx.z()<<" "<<ekin<<G4endl;
		particleGun->SetParticleEnergy(ekin);
		particleGun->SetParticlePosition(vtx);
		//particleGun->SetParticleTime(0.0)
		particleGun->SetParticleMomentumDirection( dirvec.at(pgi) );
		particleGun->GeneratePrimaryVertex(anEvent);
	      
		// true : Generate vertex in Rock , false : Generate vertex in WC tank
		SetGenerateVertexInRock(false);
	      }
	    }
	 


	    // Read in MRD data if bobfile=true
	    //	    token = readInLine(inputFile, lineSize, inBuf);
	    if(token[0]=="headerend"){

	      // read in summary line - total depth, total E, total Cherenk E
	      token = readInLine(inputFile, lineSize, inBuf);
	      // G4cout<<token[0]<<G4endl;
	      if(token[0]=="T"){
		nMRDlayers =  atoi(token[1]);
		mMRDtotEdep = atof(token[2]);
		mMRDtotEchdep = atof(token[2]);
		// G4cout<<"number of layers: "<<numlayers<<G4endl;
	      } else {
		G4cout<<"PROBLEM WITH FILE!!! No T line"<<G4endl;
	      }
	      
	      G4int layercount=0;

	      while ( token=readInLine(inputFile, lineSize, inBuf),
		      ((token[0] == "V")||(token[0] == "H")) )
		{

		  if(token[0]=="V") mrdstrp_orient = 0;
		  if(token[0]=="H") mrdstrp_orient = 1;

		  mMRDhitlayer[layercount] = atoi(token[1]);
		  mMRDhitorientation[layercount] = mrdstrp_orient;

		  mMRDhitEdep[layercount]=atof(token[2]);
		  mMRDhitEchdep[layercount]=atof(token[3]);
		  mMRDhitx[layercount]=1.0;
		  mMRDhity[layercount]=1.0;
		  mMRDhitz[layercount]=1.0;

		  // G4cout<<token[0]<<" "<<token[1]<<" "<<token[2]<<G4endl;


		  layercount++;
		}

	      nMRDlayers = layercount;


	      //token=readInLine(inputFile, lineSize, inBuf);

	      if(layercount!=nMRDlayers){
		G4cout<<"PROBLEM WITH FILE!!! numlayers don't match "<<layercount<<" "<<numlayers<<G4endl;
	      }

	    } else {
	      G4cout<<"PROBLEM WITH FILE!!! No end Header"<<G4endl;
	    }


	  }
	tcardfile->Fill();
      }
    else 
      {    // old muline format  
	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
		  >> xDir >> yDir >> zDir;
	
	//	G4double random_z = ((myDetector->GetWaterTubePosition())
	//			     - .5*(myDetector->GetWaterTubeLength()) 
	//			     + 1.*m + 15.0*m*G4UniformRand())/m;
	
	G4double random_z = 55.;
	zPos = random_z;
	vtx = G4ThreeVector(xPos, yPos, random_z);
	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

	particleGun->SetParticleEnergy(energy*MeV);
	particleGun->SetParticlePosition(vtx);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
  }

  else if (useNormalEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double mass       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

//     particleGun->SetParticleEnergy(E);
//     particleGun->SetParticlePosition(vtx);
//     particleGun->SetParticleMomentumDirection(dir);

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  
  else if (useLaserEvt)
  {
    //T. Akiri: Create the GPS LASER event
    MyGPS->GeneratePrimaryVertex(anEvent);
    
    G4ThreeVector P    = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    vtx                = anEvent->GetPrimaryVertex()->GetPosition();
    G4int pdg          = anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
    
    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P)));
    
    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  
  else if (useOpticalPhotEvt)
  {	// manual optical photon firing
   OpPhotonGun->GeneratePrimaryVertex(anEvent);
   
    G4ThreeVector P   = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    vtx = anEvent->GetPrimaryVertex()->GetPosition();
    G4double mass        = anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg         = anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  
  else if (useGeantinoEvt)
  {	// manual geantino firing
   GeantinoGun->GeneratePrimaryVertex(anEvent);
   
    G4ThreeVector P   = anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4double mass     = anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg         = anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
    vtx               = anEvent->GetPrimaryVertex()->GetPosition();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(mass*mass));

    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  
  else if (useBeamEvt)
  {     // read in from ROOT file output of rob hatcher's annie dirt simulations
//  	if(!inputfilereader->Next()){
//  		G4cout<<" END OF INPUT FILE!!!!!!"<<G4endl;
//  		inputfilereader->SetEntry(0);
//  	}
//  
//      if (inputfilereader->GetEntryStatus() == kEntryValid) {
//         G4cout << "Loaded entry " << inputfilereader->GetCurrentEntry() << '\n';
//      } else {
//         switch (inputfilereader->GetEntryStatus()) {
//         kEntryValid:
//            // Handled above.
//            break;
//         kEntryNotLoaded:
//            G4cerr << "Error: TTreeReader has not loaded any data yet!\n";
//            break;
//         kEntryNoTree:
//            G4cerr << "Error: TTreeReader cannot find a tree names \"MyTree\"!\n";
//            break;
//         kEntryNotFound:
//            // Can't really happen as TTreeReader::Next() knows when to stop.
//            G4cerr << "Error: The entry number doe not exist\n";
//            break;
//         kEntryChainSetupError:
//            G4cerr << "Error: TTreeReader cannot access a chain element, e.g. file without the tree\n";
//            break;
//         kEntryChainFileError:
//            G4cerr << "Error: TTreeReader cannot open a chain element, e.g. missing file\n";
//            break;
//         kEntryDictionaryError:
//            G4cerr << "Error: TTreeReader cannot find the dictionary for some data\n";
//            break;
//         }
//         G4Exception("WCLitePrimaryGeneratorAction::GeneratePrimaries","Expl01",FatalException,"Error reading input file");
//         return;
//      }
//  	G4cout<<"getting entry"<<G4endl;
//  	for(int i=0;i<(**rdrntankp);i++){
//  	// generate each of the primaries from this event
//  		G4cout<<"entry got"<<G4endl;
//  		
//  		G4ThreeVector rdrvtx = G4ThreeVector((*rdrvxp)[i], (*rdrvyp)[i], (*rdrvzp)[i]);
//  		G4ThreeVector rdrdir = G4ThreeVector((*rdrpxp)[i], (*rdrpyp)[i], (*rdrpzp)[i]);
//  		rdrdir.unit();		// normalise to unit vector
//  		G4cout<<"setting particle gun"<<G4endl;
//  		particleGun->SetParticleDefinition( particleTable->FindParticle( (*rdrpdgtankp)[i] ) );	//FindParticle finds by either name or pdg.
//  		particleGun->SetParticleEnergy((*rdrkEp)[i]);
//  		particleGun->SetParticlePosition(rdrvtx);
//  		particleGun->SetParticleTime((*rdrvtp)[i]);
//  		particleGun->SetParticleMomentumDirection(rdrdir);
//  		particleGun->GeneratePrimaryVertex(anEvent);	//anEvent is provided by geant4 when invoking the method
//  		
//  	}
		
		Long64_t localEntry = inputdata->LoadTree(inputEntry);
		if(localEntry<0){
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@REACHED END OF INPUT FILE!!!#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			G4cout<<"@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#"<<G4endl;
			inputEntry=0;
			localEntry = inputdata->LoadTree(inputEntry);
		}
		
		Int_t nextTreeNumber = inputdata->GetTreeNumber();
		if(treeNumber!=nextTreeNumber){
			G4cout<<"Reached end of Tree. Last entries' tree number was "<<treeNumber<<", this entries' tree number is "<<nextTreeNumber<<G4endl;
			G4cout<<"Getting new tree branches"<<G4endl;
//			inputdata->SetBranchAddress("run",&runbranchval,&runBranch);
//			inputdata->SetBranchAddress("ntank",&ntankbranchval,&nTankBranch);
			vtxxBranch=inputdata->GetBranch("vx");
			vtxyBranch=inputdata->GetBranch("vy");
			vtxzBranch=inputdata->GetBranch("vz");
			vtxtBranch=inputdata->GetBranch("vt");
			pxBranch=inputdata->GetBranch("px");
			pyBranch=inputdata->GetBranch("py");
			pzBranch=inputdata->GetBranch("pz");
			EBranch=inputdata->GetBranch("E");
			KEBranch=inputdata->GetBranch("kE");
			pdgBranch=inputdata->GetBranch("pdgtank");
			
			if(runBranch==0||nTankBranch==0||vtxxBranch==0||vtxyBranch==0||vtxzBranch==0||vtxtBranch==0||pxBranch==0||pyBranch==0||pzBranch==0||EBranch==0||KEBranch==0||pdgBranch==0){
				G4cout<<"BRANCHES ARE ZOMBIES ARGH!"<<G4endl;
			} else { G4cout<<"entries in this tree: "<<vtxxBranch->GetEntries()<<G4endl; }
			
			treeNumber=nextTreeNumber;
		}
		
		G4cout<<"Loading primaries from entry "<<inputEntry<<", localentry "<<localEntry<<G4endl;
		//runBranch->GetEntry(localEntry);
		//G4cout<<"Run is "<<runbranchval<<G4endl;
		nTankBranch->GetEntry(localEntry);
		G4cout<<"This entry has "<<ntankbranchval<<" primaries"<<G4endl;
		
		if(vtxxbranchval){delete[] vtxxbranchval;}
		if(vtxybranchval){delete[] vtxybranchval;}
		if(vtxzbranchval){delete[] vtxzbranchval;}
		if(vtxtbranchval){delete[] vtxtbranchval;}
		if(pxbranchval){delete[] pxbranchval;}
		if(pybranchval){delete[] pybranchval;}
		if(pzbranchval){delete[] pzbranchval;}
		if(ebranchval){delete[] ebranchval;}
		if(kebranchval){delete[] kebranchval;}
		if(pdgbranchval){delete[] pdgbranchval;}
	
		vtxxbranchval = new Double_t[ntankbranchval];
		vtxybranchval = new Double_t[ntankbranchval];
		vtxzbranchval = new Double_t[ntankbranchval];
		vtxtbranchval = new Double_t[ntankbranchval];
		pxbranchval = new Double_t[ntankbranchval];
		pybranchval = new Double_t[ntankbranchval];
		pzbranchval = new Double_t[ntankbranchval];
		ebranchval = new Double_t[ntankbranchval];
		kebranchval = new Double_t[ntankbranchval];
		pdgbranchval = new Int_t[ntankbranchval];
		
		if(vtxxbranchval==0||vtxybranchval==0||vtxzbranchval==0||vtxtbranchval==0||pxbranchval==0||pybranchval==0||pzbranchval==0||ebranchval==0||kebranchval==0||pdgbranchval==0){
			G4cout<<"Arrays are zombies!"<<G4endl;
		}
		
		//G4cout<<"Setting branch addresses"<<G4endl;
		vtxxBranch->SetAddress(vtxxbranchval);
		vtxyBranch->SetAddress(vtxybranchval);
		vtxzBranch->SetAddress(vtxzbranchval);
		vtxtBranch->SetAddress(vtxtbranchval);
		pxBranch->SetAddress(pxbranchval);
		pyBranch->SetAddress(pybranchval);
		pzBranch->SetAddress(pzbranchval);
		EBranch->SetAddress(ebranchval);
		KEBranch->SetAddress(kebranchval);
		pdgBranch->SetAddress(pdgbranchval);
	
		//G4cout<<"Getting primary arrays"<<G4endl;
		vtxxBranch->GetEntry(localEntry);
		vtxyBranch->GetEntry(localEntry);
		vtxzBranch->GetEntry(localEntry);
		vtxtBranch->GetEntry(localEntry);
		pxBranch->GetEntry(localEntry);
		pyBranch->GetEntry(localEntry);
		pzBranch->GetEntry(localEntry);
		EBranch->GetEntry(localEntry);
		KEBranch->GetEntry(localEntry);
		pdgBranch->GetEntry(localEntry);
		
		//G4cout<<"Looping over primaries"<<G4endl;
		for(int i=0;i<ntankbranchval;i++){
			//G4cout<<"Loading details of primary "<<i<<G4endl;
			vtxxval=vtxxbranchval[i]*cm;
			vtxyval=vtxybranchval[i]*cm;
			vtxzval=vtxzbranchval[i]*cm;
			vtxtval=vtxtbranchval[i]*ns;
			pxval=pxbranchval[i]*GeV;
			pyval=pybranchval[i]*GeV;
			pzval=pzbranchval[i]*GeV;
			eval=ebranchval[i]*GeV;
			keval=kebranchval[i]*GeV;
			pdgval=pdgbranchval[i];
			//G4cout<<"particle code is "<<pdgval<<G4endl;
			//G4ParticleDefinition* parttype = particleTable->FindParticle(pdgval);
			//if(parttype==0){G4cout<<"particle type not found"<<G4endl;} else {G4cout<<"Particle is a "<<parttype->GetParticleName()<<G4endl;}
			
			G4ThreeVector thevtx = G4ThreeVector(vtxxval, vtxyval, vtxzval);
			G4ThreeVector thepdir = G4ThreeVector(pxval, pyval, pzval);
			thepdir.unit();		// normalise to unit vector
			particleGun->SetParticleDefinition( particleTable->FindParticle(pdgval) );	//FindParticle finds by either name or pdg.
			particleGun->SetParticleEnergy(eval*MeV);
			particleGun->SetParticlePosition(thevtx);
			particleGun->SetParticleTime(vtxtval);
			particleGun->SetParticleMomentumDirection(thepdir);
			particleGun->GeneratePrimaryVertex(anEvent);	//anEvent is provided by geant4 when invoking the method
			//G4cout<<"Vertex set"<<G4endl;
		}
		
		inputEntry++;
		
	}
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// function definitions

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  int startToken = 0, endToken; // Pointers to the token pos unsigned 
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {
    
      while( startToken < input.size() )
	{
	  // Find the start of token
	  startToken = input.find_first_not_of( separators, startToken );
      
	  // If found...
	  if( startToken != input.npos ) 	// always false due to limited range of datatype!!
	    {
	      // Find end of token
	      endToken = input.find_first_of( separators, startToken );
	      if( endToken == input.npos )
		// If there was no end of token, assign it to the end of string
		endToken = input.size();
        
	      // Extract token
	      tokens.push_back( input.substr( startToken, endToken - startToken ) );
        
	      // Update startToken
	      startToken = endToken;
	    }
	}
    }
  
  return tokens;
}

//void LoadFiles(const char* dirname, TTree* tf){

//   const char* ext = ".root";
//   TSystemDirectory dir(dirname, dirname);
//   TList *files = dir.GetListOfFiles();
//   
//   if (files) {
//      TSystemFile *file;
//      TString fname;
//      TIter next(files);
//      while ((file=(TSystemFile*)next())) {
//         fname = file->GetName();
//         if (!file->IsDirectory() && fname.EndsWith(ext)) {
//            //cout << fname.Data() << endl;
//            try{
//                    std::smatch submatches;
//                    std::regex theexpression ("annie_tank_flux.([0-9]+).root");
//                    int filematched = std::regex_match ((std::string)fname, submatches, theexpression);
//                    if(filematched){
////                            std::string runstring = (std::string)submatches[1];
////                            std::string subrunstring = (std::string)submatches[2];
////                            std::string partstring = (std::string)submatches[3];
////                            int run = stoi(runstring);
////                            int subrun = stoi(subrunstring);
////                            int part = stoi(partstring);
////
////                            if(run>0)){
////                                //cout<<"adding: "<<dirname+fname<<endl;
//                                tf->Add(dirname+fname);
////                            } else {
////                            //	cout<<"skipping run "<<run<<" file"<<endl;
////                            }
//                    }
//            } catch(const std::out_of_range& oor){continue;}
//         }
//      }
//  }
//  
//}

//}

