{
  // TheLooks::set();
  TString f_name_nom  = "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasNominal-69200.root";  
  TString f_name_up = "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasUp-72383.root";
  TString f_name_dn  = "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasDown-66017.root";
  
  TFile* f_nom = new TFile(f_name_nom,"READ");
  TH1* h_nom = (TH1*) f_nom->Get("pileup");
  h_nom->SetDirectory(0);
  h_nom->SetName("nom");
  h_nom->UseCurrentStyle();
  h_nom->SetLineColor(kBlack);
  h_nom->SetLineWidth(3);
  h_nom->Scale(1./h_nom->Integral());
  h_nom->GetXaxis()->SetTitle("N(PU_{true})");
  h_nom->SetTitle("DataPileupHistogram_Run2018_294927-306462_13TeV_PromptReco_Collisions18");

  TFile* f_up  = new TFile(f_name_up ,"READ");
  TH1* h_up = (TH1*) f_up->Get("pileup");
  h_up->SetDirectory(0);
  h_up->SetName("up");
  h_up->SetLineColor(kBlue);
  h_up->SetLineWidth(3);
  h_up->Scale(1./h_up->Integral());

  TFile* f_dn  = new TFile(f_name_dn ,"READ");
  TH1* h_dn = (TH1*) f_dn->Get("pileup");
  h_dn->SetDirectory(0);
  h_dn->SetName("dn");
  h_dn->SetLineColor(kBlue);
  h_dn->SetLineWidth(3);
  h_dn->Scale(1./h_dn->Integral());

  TCanvas* can = new TCanvas("can","",600,600);
  can->cd();
  h_nom->Draw("HIST");
  h_up->Draw("HISTsame");
  h_dn->Draw("HISTsame");
  h_nom->Draw("HISTsame");

  TLegend* leg = new TLegend(0.5,0.78,0.85,0.9);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(h_nom,"nominal #sigma_{mb} = 69.2 mb","L");
  leg->AddEntry(h_up,"up/down #sigma_{mb} +/-4.6%","L");
  leg->Draw("same");

  gPad->RedrawAxis();
  can->SaveAs("DataPileupHistogram_Run2018_294927-306462_13TeV_PromptReco_Collisions18.pdf");
}
