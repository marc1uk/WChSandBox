{
  TString fname;
  fname+="test.root";
  
  TFile* tf = new TFile(fname,"READ");
  TH1D* qval = (TH1D*) tf->Get("q");
  TH1D* qval_mrd = (TH1D*) tf->Get("qmrd");

  TH1D* eff = qval_mrd->Clone();
  TH1D* effd = qval->Clone();
  effd->Sumw2();
  eff->Sumw2();

  eff->Divide(effd);

  qval->Draw();
  qval_mrd->SetLineColor(2);
  qval_mrd->Draw("SAME");

  TCanvas* tc = new TCanvas();
  eff->Draw("E");

    


  }
