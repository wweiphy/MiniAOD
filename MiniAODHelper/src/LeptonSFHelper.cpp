#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

//PUBLIC
LeptonSFHelper::LeptonSFHelper( const edm::ParameterSet& iConfig){

  //std::cout << "Initializing Lepton scale factors" << std::endl;
  if( iConfig.existsAs<edm::ParameterSet>("leptonTriggerSFInfos",true) ) {
    const edm::ParameterSet leptonTriggerSFInfos = iConfig.getParameter<edm::ParameterSet>("leptonTriggerSFInfos");
    electron_TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/" + leptonTriggerSFInfos.getParameter<std::string>("elecFileName");
    electron_TRIGGERhistname = leptonTriggerSFInfos.getParameter<std::string>("elecHistName");
    muon_TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/" + leptonTriggerSFInfos.getParameter<std::string>("muonFileName");
    muon_TRIGGERhistname = leptonTriggerSFInfos.getParameter<std::string>("muonHistName");

    SetElectronHistos( );
    SetMuonHistos( );
    SetElectronElectronHistos( );
    SetMuonMuonHistos( );
    SetElectronMuonHistos( );
    
  } 

}

LeptonSFHelper::~LeptonSFHelper( ){
  if(f_muon_TRIGGERSF) delete f_muon_TRIGGERSF;
  if(f_electron_TRIGGERSF) delete f_electron_TRIGGERSF;
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
      //pt = Electron.userFloat("ptBeforeRun2Calibration");
      //eta = Electron.superCluster()->eta();
      //ElectronTriggerSF = ElectronTriggerSF * GetLeptonTriggerSF(pt, eta, 0, h_ele_TRIGGER_abseta_pt_ratio);
      //ElectronTriggerSF_Up = ElectronTriggerSF_Up  * GetLeptonTriggerSF(pt, eta, 1, h_ele_TRIGGER_abseta_pt_ratio);
      //ElectronTriggerSF_Down = ElectronTriggerSF_Down * GetLeptonTriggerSF(pt, eta, -1, h_ele_TRIGGER_abseta_pt_ratio);
    }
    
    else {
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

  //ScaleFactorMap["ElectronTriggerSF"] = ElectronTriggerSF;
  //ScaleFactorMap["ElectronTriggerSF_Up"] = ElectronTriggerSF_Up;
  //ScaleFactorMap["ElectronTriggerSF_Down"] = ElectronTriggerSF_Down;

  ScaleFactorMap["MuonTriggerSF"] = MuonTriggerSF;
  ScaleFactorMap["MuonTriggerSF_Up"] = MuonTriggerSF_Up;
  ScaleFactorMap["MuonTriggerSF_Down"] = MuonTriggerSF_Down;
  
  ScaleFactorMap["MuonMuonTriggerSF"] = MuonMuonTriggerSF;

  ScaleFactorMap["ElectronMuonTriggerSF"] = ElectronMuonTriggerSF;

  // ScaleFactorMap["ElectronSF"]= ElectronTriggerSF;
  // ScaleFactorMap["ElectronSF_Up"]= ElectronTriggerSF_Up;
  // ScaleFactorMap["ElectronSF_Down"]= ElectronTriggerSF_Down;

  // ScaleFactorMap["MuonSF"]=  MuonTriggerSF;
  // ScaleFactorMap["MuonSF_Up"]= MuonTriggerSF_Up;
  // ScaleFactorMap["MuonSF_Down"]= MuonTriggerSF_Down;

  //ScaleFactorMap["LeptonTriggerSF"]= ScaleFactorMap["ElectronTriggerSF"] * ScaleFactorMap["MuonTriggerSF"];
  //ScaleFactorMap["LeptonTriggerSF_Up"]= ScaleFactorMap["ElectronTriggerSF_Up"] * ScaleFactorMap["MuonTriggerSF_Up"];
  //ScaleFactorMap["LeptonTriggerSF_Down"]= ScaleFactorMap["ElectronTriggerSF_Down"] * ScaleFactorMap["MuonTriggerSF_Down"];

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
    pt = std::max(xmin+0.1,lepton_pt);
    pt = std::min(xmax-0.1,pt);
    eta = std::max(ymin+0.1,lepton_eta);
    eta = std::min(ymax-0.1,eta);


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

  
  f_electron_TRIGGERSF = new TFile(electron_TRIGGERinputFile.c_str(),"READ");
  
  h_ele_TRIGGER_abseta_pt_ratio = (TH2F*)f_electron_TRIGGERSF->Get(electron_TRIGGERhistname.c_str());

}

void LeptonSFHelper::SetMuonHistos( ){
  
  f_muon_TRIGGERSF = new TFile(muon_TRIGGERinputFile.c_str(),"READ");
  
  h_mu_TRIGGER_abseta_pt= (TH2F*)f_muon_TRIGGERSF->Get(muon_TRIGGERhistname.c_str());  
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
