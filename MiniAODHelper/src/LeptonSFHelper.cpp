#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

//PUBLIC
LeptonSFHelper::LeptonSFHelper( ){

  //std::cout << "Initializing Lepton scale factors" << std::endl;

  SetElectronHistos( );
  SetMuonHistos( );
  SetElectronElectronHistos( );
  SetMuonMuonHistos( );
  SetElectronMuonHistos( );

  electronLowPtRangeCut=20.0;
  electronMaxPt = 150.0;
  electronMinPt = 20.0;
  electronMinPtLowPt = 10;
  electronMaxPtLowPt = 19.9; 
  electronMaxPtHigh= 201.0;
  electronMaxPtHigher= 499.0;  //TH2 histos from Fall17 are binned up to 500 GeV
  electronMaxEta=2.49;
  electronMaxEtaLow=2.19;
  
  
  
  muonMaxPt = 119.0;
  muonMaxPtHigh = 1199.;       //TH2 Trigger SF histos from Fall17 are binned up to 1200 GEV 
  muonMinPt = 20.0;
  muonMinPtHigh = 29.0;
  
  muonMaxEta = 2.39;

}

LeptonSFHelper::~LeptonSFHelper( ){

}

std::map< std::string, float >  LeptonSFHelper::GetLeptonSF( const std::vector< pat::Electron >& Electrons,
							     const std::vector< pat::Muon >& Muons  ) {


  std::map< std::string , float > ScaleFactorMap;

  float ElectronTriggerSF = 1.0;
  float ElectronTriggerSF_Up = 1.0;
  float ElectronTriggerSF_Down = 1.0;
  
  float MuonTriggerSF = 1.0;
  float MuonTriggerSF_Up = 1.0;
  float MuonTriggerSF_Down = 1.0;
  
  float ElectronElectronTriggerSF = 1.0;
  float MuonMuonTriggerSF = 1.0;
  float ElectronMuonTriggerSF = 1.0;

  double pt = 0;
  double eta = 0;
  for (auto Electron: Electrons){ //Electron is of type pat::Electron
      
    if(Electron.hasUserFloat("ptBeforeRun2Calibration")) {
      pt = Electron.userFloat("ptBeforeRun2Calibration");
      eta = Electron.superCluster()->eta();
      ElectronTriggerSF = ElectronTriggerSF * GetLeptonTriggerSF(pt, eta, 0, h_ele_TRIGGER_abseta_pt_ratio);
      ElectronTriggerSF_Up = ElectronTriggerSF_Up  * GetLeptonTriggerSF(pt, eta, 1, h_ele_TRIGGER_abseta_pt_ratio);
      ElectronTriggerSF_Down = ElectronTriggerSF_Down * GetLeptonTriggerSF(pt, eta, -1, h_ele_TRIGGER_abseta_pt_ratio);
    }
    
    else {
        
        // ElectronTriggerSF = ElectronTriggerSF * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 0, "Trigger");
        // ElectronTriggerSF_Up = ElectronTriggerSF_Up  * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 1, "Trigger");
        // ElectronTriggerSF_Down = ElectronTriggerSF_Down * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), -1, "Trigger");
        
        std::cerr << "ERROR: could not get ElectronTriggerSF because muon lacks property 'ptBeforeRun2Calibration'!" << std::endl;
        throw std::exception();
    }


  }
  
  for (auto Muon: Muons){ //Muon is of type pat::Muon
    // trigger SFs are derived without corrections
    if(Muon.hasUserFloat("PtbeforeRC")) {
      pt = Muon.userFloat("PtbeforeRC");
      eta = fabs(Muon.eta());
      MuonTriggerSF = MuonTriggerSF * GetLeptonTriggerSF(pt, eta, 0, h_mu_TRIGGER_abseta_pt);
      MuonTriggerSF_Up = MuonTriggerSF_Up  * GetLeptonTriggerSF(pt, eta, 1, h_mu_TRIGGER_abseta_pt);
      MuonTriggerSF_Down = MuonTriggerSF_Down * GetLeptonTriggerSF(pt, eta, -1, h_mu_TRIGGER_abseta_pt);
    }
    
    else {
        // MuonTriggerSF = MuonTriggerSF * GetMuonSF(Muon.pt(), Muon.eta(), 0);
        // MuonTriggerSF_Up = MuonTriggerSF_Up  * GetMuonSF(Muon.pt(), Muon.eta(), 1);
        // MuonTriggerSF_Down = MuonTriggerSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1);
      std::cerr << "ERROR: could not get MuonTriggerSF because muon lacks property 'PtbeforeRC'!" << std::endl;
      throw std::exception();
        
    }

  }

  //std::cout << "Anzahl El plus Mu " << Electrons.size()+Muons.size() << std::endl;

  if(Electrons.size()+Muons.size()==2) {
    if(Electrons.size()==2) {

	//std::cout << "zwei Elektronen ! " << std::endl;
	//ElectronElectronTriggerSF = ElectronElectronTriggerSF * GetElectronElectronSF(Electrons.at(0).eta(), Electrons.at(1).eta(), 0, "Trigger");


    }
    else if(Muons.size()==2) {

	//std::cout << "zwei Muonen ! " << std::endl;
	//MuonMuonTriggerSF = MuonMuonTriggerSF * GetMuonMuonSF(Muons.at(0).eta(), Muons.at(1).eta(), 0, "Trigger");


    }
    else {

	//std::cout << "ein Ele und ein Mu ! " << std::endl;
	//ElectronMuonTriggerSF = ElectronMuonTriggerSF * GetElectronMuonSF(Electrons.at(0).eta(), Muons.at(0).eta(), 0, "Trigger");


    }
  }

  //std::cout << ElectronElectronTriggerSF << "  " << MuonMuonTriggerSF << "  " << ElectronMuonTriggerSF << std::endl;

  ScaleFactorMap["ElectronSFTrigger"] = ElectronTriggerSF;
  ScaleFactorMap["ElectronSFTrigger_Up"] = ElectronTriggerSF_Up;
  ScaleFactorMap["ElectronSFTrigger_Down"] = ElectronTriggerSF_Down;

  // ScaleFactorMap["MuonSFTrigger"] = MuonTriggerSF;
  // ScaleFactorMap["MuonSFTrigger_Up"] = MuonTriggerSF_Up;
  // ScaleFactorMap["MuonSFTrigger_Down"] = MuonTriggerSF_Down;
  
  ScaleFactorMap["MuonMuonTriggerSF"] = MuonMuonTriggerSF;

  ScaleFactorMap["ElectronMuonTriggerSF"] = ElectronMuonTriggerSF;

  ScaleFactorMap["ElectronSF"]= ElectronTriggerSF;
  ScaleFactorMap["ElectronSF_Up"]= ElectronTriggerSF_Up;
  ScaleFactorMap["ElectronSF_Down"]= ElectronTriggerSF_Down;

  ScaleFactorMap["MuonSF"]=  MuonTriggerSF;
  ScaleFactorMap["MuonSF_Up"]= MuonTriggerSF_Up;
  ScaleFactorMap["MuonSF_Down"]= MuonTriggerSF_Down;

  ScaleFactorMap["LeptonSF"]= ScaleFactorMap["ElectronSF"] * ScaleFactorMap["MuonSF"];
  ScaleFactorMap["LeptonSF_Up"]= ScaleFactorMap["ElectronSF_Up"] * ScaleFactorMap["MuonSF_Up"];
  ScaleFactorMap["LeptonSF_Down"]= ScaleFactorMap["ElectronSF_Down"] * ScaleFactorMap["MuonSF_Down"];

  return ScaleFactorMap;
}

float LeptonSFHelper::GetLeptonTriggerSF(  const double& lepton_pt , const double& lepton_eta , const int& syst, TH2F* h_SFs){
  if ( lepton_pt == 0.0 ){ return 1.0; }

  int thisBin=0;
  float nomval = 0;
  float error = 0;
  float upval = 0;
  float downval= 0;
  double eta = 0;
  double pt = 0;
  if ( h_SFs ){
    // determine the ranges of the given TH2Fs
    auto xmin = h_SFs->GetXaxis()->GetXmin();
    auto xmax = h_SFs->GetXaxis()->GetXmax();
    auto ymin = h_SFs->GetYaxis()->GetXmin();
    auto ymax = h_SFs->GetYaxis()->GetXmax();
    
    // make sure to stay within the range ot the histograms
    eta = std::max(xmin+0.1,lepton_eta);
    eta = std::min(xmax-0.1,lepton_eta);
    pt = std::max(ymin+0.1,lepton_pt);
    pt = std::min(ymax-0.1,lepton_pt);


    thisBin = h_SFs->FindBin(  pt, eta  );
    nomval=h_SFs->GetBinContent( thisBin );
    error=h_SFs->GetBinError( thisBin );
    upval=( nomval+error );
    downval=( nomval-error );
  }
  
  else {
    std::cerr << "ERROR: Could not load histogram for lepton trigger SFs" << std::endl;
    throw std::exception();

  }


  if ( syst==-1 ){ return downval; }
  else if ( syst==1 ){ return upval; }
  else { return nomval; }

}

float LeptonSFHelper::GetElectronSF(  float electronPt , float electronEta , int syst ) {
  if ( electronPt == 0.0 ){ return 1.0; }

  int thisBin=0;

  // restrict electron eta 
  float searchEta=electronEta;
  if(searchEta<0 and searchEta<=-electronMaxEta){searchEta=-electronMaxEta;}
  if(searchEta>0 and searchEta>=electronMaxEta){searchEta=electronMaxEta;}
  
  if(searchEta<0 and searchEta<=-electronMaxEtaLow){searchEta=-electronMaxEtaLow;}
  if(searchEta>0 and searchEta>=electronMaxEtaLow){searchEta=electronMaxEtaLow;}
    
  // restrict electron pT
  float searchPt=electronPt;
  if(searchPt>electronLowPtRangeCut){
    if(searchPt>=electronMaxPtHigher) {searchPt=electronMaxPtHigher;}; // if e_pt >= 500 go to last bin by setting searchPt to 499
    if(searchPt<electronMinPt){searchPt=electronMinPt;}; // if e_pt < 20 go to first bin by setting searchPt to 20
  }
  else{
  // these are now for the low pt Reco SF 
  if(searchPt>=electronMaxPtLowPt) {searchPt=electronMaxPtLowPt;}; // if e_pt >= 500 go to last bin by setting searchPtLowPt to 499
  if(searchPt<electronMinPtLowPt){searchPt=electronMinPtLowPt;}; // if e_pt < 20 go to first bin by setting searchPtLowPt to 20
  }

  float nomval = 0;
  float error = 0;
  float upval = 0;
  float downval= 0;
  float nomvalBtoF = 0;
  float errorBtoF = 0;
  float upvalBtoF = 0;
  float downvalBtoF= 0;

  if ( h_ele_TRIGGER_abseta_pt_ratio ){
    // std::cout << "getting Trigger SF\n";
    thisBin = h_ele_TRIGGER_abseta_pt_ratio->FindBin( searchPt, searchEta );
    nomval=h_ele_TRIGGER_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_TRIGGER_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error;
    downval=nomval-error;

  }
  else {

    std::cerr << "ERROR: Could not load histogram for electron trigger SFs" << std::endl;
    throw std::exception();

  }

  if ( syst==-1 ){ return downval; }
  else if ( syst==1 ){ return upval; }
  else { return nomval; }

}

float LeptonSFHelper::GetMuonSF(  const double& muonPt , const double& muonEta , const int& syst){
  if ( muonPt == 0.0 ){ return 1.0; }

  int thisBin=0;

  auto eta=fabs( muonEta );
  auto pt= muonPt;
  auto h_SFs = h_mu_TRIGGER_abseta_pt; //TODO: change this!

  float nomval = 0;
  //float error = 0;
  float upval = 0;
  float downval= 0;
  float nomvalBtoF = 0;
  float errorBtoF = 0;
  float upvalBtoF = 0;
  float downvalBtoF= 0;
  if ( h_SFs ){
    // determine the ranges of the given TH2Fs
    auto xmin = h_SFs->GetXaxis()->GetXmin();
    auto xmax = h_SFs->GetXaxis()->GetXmax();
    auto ymin = h_SFs->GetYaxis()->GetXmin();
    auto ymax = h_SFs->GetYaxis()->GetXmax();
    
    // make sure to stay within the range ot the histograms
    eta = std::max(xmin+0.1,eta);
    eta = std::min(xmax-0.1,eta);
    pt = std::max(ymin+0.1,pt);
    pt = std::min(ymax-0.1,pt);


    thisBin = h_SFs->FindBin(  pt, eta  );
    nomvalBtoF=h_SFs->GetBinContent( thisBin );
    errorBtoF=h_SFs->GetBinError( thisBin );
    upvalBtoF=( nomvalBtoF+errorBtoF );
    downvalBtoF=( nomvalBtoF-errorBtoF );

    nomval=nomvalBtoF;
    upval=upvalBtoF;
    downval=downvalBtoF;
    
  }
  
  else {

    std::cerr << "ERROR: Could not load histogram for muon trigger SFs" << std::endl;
    throw std::exception();

  }


  if ( syst==-1 ){ return downval; }
  else if ( syst==1 ){ return upval; }
  else { return nomval; }


}
float LeptonSFHelper::GetElectronElectronSF(  float electronEta1 , float electronEta2 , int syst , std::string type  ) {

  int thisBin=0;

  float searchEta1=fabs(electronEta1);
  float searchEta2=fabs(electronEta2);

  float nomval = 0;
  if ( type == "Trigger" ){

    thisBin = h_ele_ele_TRIGGER_abseta_abseta->FindBin( searchEta1 , searchEta2 );
    nomval = h_ele_ele_TRIGGER_abseta_abseta->GetBinContent( thisBin );

  }
  return nomval;
}
float LeptonSFHelper::GetMuonMuonSF(  float muonEta1 , float muonEta2 , int syst , std::string type  ) {
  int thisBin=0;

  float searchEta1=fabs(muonEta1);
  float searchEta2=fabs(muonEta2);

  float nomval = 0;
  if ( type == "Trigger" ){

    thisBin = h_mu_mu_TRIGGER_abseta_abseta->FindBin( searchEta1 , searchEta2 );
    nomval = h_mu_mu_TRIGGER_abseta_abseta->GetBinContent( thisBin );

  }
  return nomval;
}
float LeptonSFHelper::GetElectronMuonSF(  float electronEta , float muonEta , int syst , std::string type  ) {
  int thisBin=0;

  float searchEta1=fabs(electronEta);
  float searchEta2=fabs(muonEta);

  float nomval = 0;
  if ( type == "Trigger" ){

    thisBin = h_ele_mu_TRIGGER_abseta_abseta->FindBin( searchEta1 , searchEta2 );
    nomval= h_ele_mu_TRIGGER_abseta_abseta->GetBinContent( thisBin );

  }
  return nomval;
}

//PRIVATE

void LeptonSFHelper::SetElectronHistos( ){

  
  TFile *f_TRIGGERSF = new TFile(electron_TRIGGERinputFile.c_str(),"READ");
  
  h_ele_TRIGGER_abseta_pt_ratio = (TH2F*)f_TRIGGERSF->Get(electron_TRIGGERhistname.c_str());

}

void LeptonSFHelper::SetMuonHistos( ){
  
  TFile *f_TRIGGERSF = new TFile(muon_TRIGGERinputFile.c_str(),"READ");
  
  h_mu_TRIGGER_abseta_pt= (TH2F*)f_TRIGGERSF->Get(muon_TRIGGERhistname.c_str());  
}

void LeptonSFHelper::SetElectronElectronHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/oct202017/" + "triggerSummary_ee_ReReco2016_ttH.root";

  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");

  h_ele_ele_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}

void LeptonSFHelper::SetMuonMuonHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/oct202017/" + "triggerSummary_mumu_ReReco2016_ttH.root";

  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");

  h_mu_mu_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}

void LeptonSFHelper::SetElectronMuonHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/oct202017/" + "triggerSummary_emu_ReReco2016_ttH.root";

  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");

  h_ele_mu_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}

int LeptonSFHelper::findPoint(TGraphAsymmErrors& graph,float& x_) {
    double x=0.;
    double y=0.;
    for(int i=0;i<graph.GetN();i++) {
        graph.GetPoint(i,x,y);
        double l=0.;
        double r=0.;
        l = x-graph.GetErrorXlow(i);
        r = x+graph.GetErrorXhigh(i);
        if((l<=x_) && (x_<r)) {return i;}
    }
    return -1;
}

float LeptonSFHelper::getValue(TGraphAsymmErrors& graph,float& x_,int syst) {
    int i = findPoint(graph,x_);
    if(i<0) {std::cerr << "x-value " << x_ << " cannot be assigned to a valid point" << std::endl;}
    double x=0.;
    double y=0.;
    graph.GetPoint(i,x,y);
    float y_=y;
    if(syst==1) {return y_+graph.GetErrorYhigh(i);}
    else if(syst==-1) {return y_-graph.GetErrorYlow(i);}
    else {return y_;}
}
