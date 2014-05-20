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

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

using std::vector;
using std::string;
using std::fstream;

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
  // Read in line break it up into tokens
  inFile.getline(inBuf,lineSize);
  return tokenize(" $", inBuf);
}

inline float atof( const string& s ) {return std::atof( s.c_str() );}
inline int   atoi( const string& s ) {return std::atoi( s.c_str() );}

WCLitePrimaryGeneratorAction::WCLitePrimaryGeneratorAction(
					  WCLiteDetectorConstruction* myDC)
  :myDetector(myDC)
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

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

  //  tcardfile->Branch("issum",&msum); 
  //  tcardfile->Branch("Nneut",&nneut);
  //  tcardfile->Branch("np",&np);


  mR = new TRandom3();

  // Initialize to zero
  mode = 0;
  vtxvol = 0;
  vtx = G4ThreeVector(0.,0.,0.);
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
 
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
  useMulineEvt = true;
  useNormalEvt = false;
}

WCLitePrimaryGeneratorAction::~WCLitePrimaryGeneratorAction()
{
  if (IsGeneratingVertexInRock()){
    G4cout << "Fraction of Rock volume is : " << G4endl;
      G4cout << " Random number generated in Rock / in Cublic = " 
             << _counterRock << "/" << _counterCublic 
             << " = " << _counterRock/(G4double)_counterCublic << G4endl;
  }
  inputFile.close();

  TFile *ttttf = new TFile("generatorcardfile.root","RECREATE");
  tcardfile->Write();
  ttttf->Close();

  delete tcardfile;
  delete particleGun;
  delete MyGPS;   //T. Akiri: Delete the GPS variable
  delete messenger;
}

void WCLitePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // We will need a particle table
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // Temporary kludge to turn on/off vector text format 

  G4bool useNuanceTextFormat = true;
  G4bool fullsim = true;
  G4int vtxchooser=2;

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

		  if( atoi(token[1]) == 14 ){
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

		  if( atoi(token[1]) == 8016 ){
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
		    {
		      G4int pdgid = atoi(token[1]);
		      G4double energy = atof(token[2])*MeV;
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


	    
	    G4double npl=-555.;
	    
	    for(int pgi=0; pgi<ntrks; pgi++){
	             
	      mpid[pgi] = idvec.at(pgi);
		    
	      particleGun->
		SetParticleDefinition( particleTable->
				       FindParticle(idvec.at(pgi)) );

	      G4double mass = 
		particleGun->GetParticleDefinition()->GetPDGMass();
	      
	      G4double energy = Evec.at(pgi);

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
	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
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
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(m*m));

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
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P)));
      
      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
}

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  unsigned int startToken = 0, endToken; // Pointers to the token pos
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {
    
      while( startToken < input.size() )
	{
	  // Find the start of token
	  startToken = input.find_first_not_of( separators, startToken );
      
	  // If found...
	  if( startToken != input.npos ) 
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

