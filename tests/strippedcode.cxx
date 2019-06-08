#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
//#include "TObject.h"
#include "TNtuple.h"

void addmrdplots(const Char_t* mrdinfile="MRDEvents.root"){

TFile mrdfile = TFile(mrdinfile);
TTree* mrdtree = (TTree*)mrdfile.Get("MRDTree");
TCanvas aCanv = TCanvas("aCanv","Title");
aCanv.Divide(4,3);
Int_t mrdnumlayers=12;
Int_t mrdpaddlesperpanel=15;
Double_t* totedepinpanels = new Double_t[mrdnumlayers];
mrdtree->SetAlias("mrdhit_panelnum",Form("TMath::Floor(mrdhit_copynum/%d)",mrdpaddlesperpanel));
mrdtree->SetAlias("mrdhit_paddlenum",Form("mrdhit_copynum-(mrdhit_panelnum*%d)",mrdpaddlesperpanel));
for(Int_t panelnum=0;panelnum<mrdnumlayers;panelnum++){
	aCanv.cd(panelnum+1);
	mrdtree->Draw(Form("mrdhit_paddlenum>>edepinpanel%d(%d,1,%d)",panelnum,mrdpaddlesperpanel,mrdpaddlesperpanel+1),Form("mrdhit_edep*(mrdhit_panelnum==%d)",panelnum));
	TCanvas* currentcanvas = (TCanvas*)aCanv.GetPrimitive(Form("aCanv_%d",panelnum+1));
	TH1F* edepinpanel = (TH1F*)currentcanvas->GetPrimitive(Form("edepinpanel%d",panelnum));
	//edepinpanel->Write();
	Double_t totedepinpanel = edepinpanel->Integral();
	totedepinpanels[panelnum]=totedepinpanel;
}
aCanv.Modified();
aCanv.SaveAs("deposition in panels.png");

mrdfile.Close();
cout<<"marco"<<endl;

}
