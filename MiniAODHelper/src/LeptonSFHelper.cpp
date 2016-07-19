#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

//PUBLIC
LeptonSFHelper::LeptonSFHelper( ){

  //std::cout << "Initializing Lepton scale factors" << std::endl;

  SetElectronHistos( );
  SetMuonHistos( );
  SetElectronElectronHistos( );
  SetMuonMuonHistos( );
  SetElectronMuonHistos( );

  electronMaxPt = 150;
  muonMaxPt = 115;

}

LeptonSFHelper::~LeptonSFHelper( ){

}

std::map< std::string, float >  LeptonSFHelper::GetLeptonSF( const std::vector< pat::Electron >& Electrons,
							     const std::vector< pat::Muon >& Muons  ) {
       

  std::map< std::string , float > ScaleFactorMap;

  float ElectronIDSF = 1.0;
  float ElectronIDSF_Up = 1.0;
  float ElectronIDSF_Down = 1.0;
  float ElectronIsoSF = 1.0;
  float ElectronIsoSF_Up = 1.0;
  float ElectronIsoSF_Down = 1.0;
  float ElectronTriggerSF = 1.0;
  float ElectronTriggerSF_Up = 1.0;
  float ElectronTriggerSF_Down = 1.0;
  float MuonIDSF = 1.0;
  float MuonIDSF_Up = 1.0;
  float MuonIDSF_Down = 1.0;
  float MuonIsoSF = 1.0;
  float MuonIsoSF_Up = 1.0;
  float MuonIsoSF_Down = 1.0;
  float MuonTriggerSF = 1.0;
  float MuonTriggerSF_Up = 1.0;
  float MuonTriggerSF_Down = 1.0;
  float ElectronElectronTriggerSF = 1.0;
  float MuonMuonTriggerSF = 1.0;
  float ElectronMuonTriggerSF = 1.0;

 
  for (auto Electron: Electrons){ //Electron is of type pat::Electron

    ElectronIDSF = ElectronIDSF * GetElectronSF(Electron.pt(), Electron.eta(), 0, "ID");
    ElectronIDSF_Up = ElectronIDSF_Up *GetElectronSF(Electron.pt(), Electron.eta(), 1, "ID");
    ElectronIDSF_Down = ElectronIDSF_Down * GetElectronSF(Electron.pt(), Electron.eta(), -1, "ID");
    
    ElectronIsoSF = ElectronIsoSF * GetElectronSF(Electron.pt(), Electron.eta(), 0, "Iso");
    ElectronIsoSF_Up = ElectronIsoSF_Up  * GetElectronSF(Electron.pt(), Electron.eta(), 1, "Iso");
    ElectronIsoSF_Down = ElectronIsoSF_Down * GetElectronSF(Electron.pt(), Electron.eta(), -1, "Iso");

    ElectronTriggerSF = ElectronTriggerSF * GetElectronSF(Electron.pt(), Electron.eta(), 0, "Trigger");
    ElectronTriggerSF_Up = ElectronTriggerSF_Up  * GetElectronSF(Electron.pt(), Electron.eta(), 1, "Trigger");
    ElectronTriggerSF_Down = ElectronTriggerSF_Down * GetElectronSF(Electron.pt(), Electron.eta(), -1, "Trigger");

  } for (auto Muon: Muons){ //Muon is of type pat::Muon
    
    MuonIDSF = MuonIDSF * GetMuonSF(Muon.pt(), Muon.eta(), 0, "ID");
    MuonIDSF_Up = MuonIDSF_Up * GetMuonSF(Muon.pt(), Muon.eta(), 1, "ID");
    MuonIDSF_Down = MuonIDSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1, "ID");

    MuonIsoSF = MuonIsoSF * GetMuonSF(Muon.pt(), Muon.eta(), 0, "Iso");
    MuonIsoSF_Up = MuonIsoSF_Up  * GetMuonSF(Muon.pt(), Muon.eta(), 1, "Iso");
    MuonIsoSF_Down = MuonIsoSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1, "Iso");

    MuonTriggerSF = MuonTriggerSF * GetMuonSF(Muon.pt(), Muon.eta(), 0, "Trigger");
    MuonTriggerSF_Up = MuonTriggerSF_Up  * GetMuonSF(Muon.pt(), Muon.eta(), 1, "Trigger");
    MuonTriggerSF_Down = MuonTriggerSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1, "Trigger");

  }
  
  //std::cout << "Anzahl El plus Mu " << Electrons.size()+Muons.size() << std::endl;
  
  if(Electrons.size()+Muons.size()==2) {
    if(Electrons.size()==2) {
      
	//std::cout << "zwei Elektronen ! " << std::endl;
	ElectronElectronTriggerSF = ElectronElectronTriggerSF * GetElectronElectronSF(Electrons.at(0).eta(), Electrons.at(1).eta(), 0, "Trigger");

      
    }
    else if(Muons.size()==2) {
      
	//std::cout << "zwei Muonen ! " << std::endl;
	MuonMuonTriggerSF = MuonMuonTriggerSF * GetMuonMuonSF(Muons.at(0).eta(), Muons.at(1).eta(), 0, "Trigger");

      
    }
    else {
      
	//std::cout << "ein Ele und ein Mu ! " << std::endl;
	ElectronMuonTriggerSF = ElectronMuonTriggerSF * GetElectronMuonSF(Electrons.at(0).eta(), Muons.at(0).eta(), 0, "Trigger");

      
    }
  }
  
  //std::cout << ElectronElectronTriggerSF << "  " << MuonMuonTriggerSF << "  " << ElectronMuonTriggerSF << std::endl;

  ScaleFactorMap["ElectronSFID"] = ElectronIDSF;
  ScaleFactorMap["ElectronSFID_Up"] = ElectronIDSF_Up;
  ScaleFactorMap["ElectronSFID_Down"] = ElectronIDSF_Down;
  ScaleFactorMap["ElectronSFIso"] = ElectronIsoSF;
  ScaleFactorMap["ElectronSFIso_Up"] = ElectronIsoSF_Up;
  ScaleFactorMap["ElectronSFIso_Down"] = ElectronIsoSF_Down;
  ScaleFactorMap["ElectronSFTrigger"] = ElectronTriggerSF;
  ScaleFactorMap["ElectronSFTrigger_Up"] = ElectronTriggerSF_Up;
  ScaleFactorMap["ElectronSFTrigger_Down"] = ElectronTriggerSF_Down;
  ScaleFactorMap["ElectronElectronTriggerSF"] =ElectronElectronTriggerSF;

  ScaleFactorMap["MuonSFID"] = MuonIDSF;
  ScaleFactorMap["MuonSFID_Up"] = MuonIDSF_Up;
  ScaleFactorMap["MuonSFID_Down"] = MuonIDSF_Down;
  ScaleFactorMap["MuonSFIso"] = MuonIsoSF;
  ScaleFactorMap["MuonSFIso_Up"] = MuonIsoSF_Up;
  ScaleFactorMap["MuonSFIso_Down"] = MuonIsoSF_Down;
  ScaleFactorMap["MuonSFTrigger"] = MuonTriggerSF;
  ScaleFactorMap["MuonSFTrigger_Up"] = MuonTriggerSF_Up;
  ScaleFactorMap["MuonSFTrigger_Down"] = MuonTriggerSF_Down;
  ScaleFactorMap["MuonMuonTriggerSF"] = MuonMuonTriggerSF;

  ScaleFactorMap["ElectronMuonTriggerSF"] = ElectronMuonTriggerSF;

  ScaleFactorMap["ElectronSF"]= ElectronIDSF * ElectronIsoSF * ElectronTriggerSF;
  ScaleFactorMap["ElectronSF_Up"]= ElectronIDSF_Up * ElectronIsoSF_Up * ElectronTriggerSF_Up;
  ScaleFactorMap["ElectronSF_Down"]= ElectronIDSF_Down * ElectronIsoSF_Down * ElectronTriggerSF_Down;

  ScaleFactorMap["MuonSF"]= MuonIDSF * MuonIsoSF * MuonTriggerSF;
  ScaleFactorMap["MuonSF_Up"]= MuonIDSF_Up * MuonIsoSF_Up * MuonTriggerSF_Up;
  ScaleFactorMap["MuonSF_Down"]= MuonIDSF_Down * MuonIsoSF_Down * MuonTriggerSF_Down;

  ScaleFactorMap["LeptonSF"]= ScaleFactorMap["ElectronSF"] * ScaleFactorMap["MuonSF"];
  ScaleFactorMap["LeptonSF_Up"]= ScaleFactorMap["ElectronSF_Up"] * ScaleFactorMap["MuonSF_Up"];
  ScaleFactorMap["LeptonSF_Down"]= ScaleFactorMap["ElectronSF_Down"] * ScaleFactorMap["MuonSF_Down"];

  return ScaleFactorMap;
}
float LeptonSFHelper::GetElectronSF(  float electronPt , float electronEta , int syst , std::string type  ) {
  if ( electronPt == 0.0 ){ return 1.0; }

  int thisBin=0;
  
  float searchEta=electronEta;
  float searchPt=TMath::Min( electronPt , electronMaxPt );
  
  float nomval = 0;
  float error = 0;
  float upval = 0;
  float downval= 0;
  

  if ( type == "ID" ){

    thisBin = h_ele_ID_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_ele_ID_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_ID_abseta_pt_ratio->GetBinError( thisBin );
    
    upval=nomval+error;
    downval=nomval-error;
   
  }
  else if ( type == "Trigger" ){

    thisBin = h_ele_TRIGGER_abseta_pt_ratio->FindBin( searchPt , searchEta );
    nomval=h_ele_TRIGGER_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_TRIGGER_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error;
    downval=nomval-error;

  }
  else if ( type == "Iso" ){

    thisBin = h_ele_ISO_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_ele_ISO_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_ISO_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error;
    downval=nomval-error;
    
  }
  else {
   
    std::cout << "Unknown Type. Supported Types are: ID, Trigger, Iso" << std::endl;
    nomval = -1;
    upval = -1;
    downval= -1;

  }

  if ( syst==-1 ){ return downval; }
  else if ( syst==1 ){ return upval; }
  else { return nomval; }

}

float LeptonSFHelper::GetMuonSF(  float muonPt , float muonEta , int syst , std::string type  ){
  if ( muonPt == 0.0 ){ return 1.0; }

  int thisBin=0;

  float searchEta=fabs( muonEta );
  float searchPt=TMath::Min( muonPt , muonMaxPt );

  float nomval = 0;
  float error = 0;
  float upval = 0;
  float downval= 0;

  if ( type == "ID" ){

    thisBin = h_mu_ID_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_mu_ID_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_mu_ID_abseta_pt_ratio->GetBinError( thisBin );
    upval=( nomval+error );
    downval=( nomval-error );
    upval=upval*( 1.0+0.01 );
    downval=downval*( 1.0-0.01 );



  }
  else if ( type == "Trigger" ){

    float mult4p2 = 0.2843;
    float mult4p3 = 0.716;

    thisBin = h_mu_TRIGGER_abseta_pt_ratio4p3->FindBin(searchEta,searchPt);
    float nomval4p3=h_mu_TRIGGER_abseta_pt_ratio4p3->GetBinContent(thisBin);
    float error4p3=h_mu_TRIGGER_abseta_pt_ratio4p3->GetBinError(thisBin);
    float upval4p3=nomval4p3+error4p3;
    float downval4p3=nomval4p3-error4p3;
    thisBin = h_mu_TRIGGER_abseta_pt_ratio4p2->FindBin(searchEta,searchPt);
    float nomval4p2=h_mu_TRIGGER_abseta_pt_ratio4p2->GetBinContent(thisBin);
    float error4p2=h_mu_TRIGGER_abseta_pt_ratio4p2->GetBinError(thisBin);
    float upval4p2=nomval4p2+error4p2;
    float downval4p2=nomval4p2-error4p2;
    
    nomval = mult4p2*nomval4p2 + mult4p3*nomval4p3;
    upval = mult4p2*upval4p2 + mult4p3*upval4p3;
    downval = mult4p2*downval4p2 + mult4p3*downval4p3;
    upval=upval*(1.0+0.005);
    downval=downval*(1.0-0.005);
    
  }
  else if ( type == "Iso" ){

    thisBin = h_mu_ISO_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_mu_ISO_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_mu_ISO_abseta_pt_ratio->GetBinError( thisBin );
    upval=( nomval+error );
    downval=( nomval-error );
    upval=upval*( 1.0+0.005 );
    downval=downval*( 1.0-0.005 );
    
      
  }
  else {
    
    std::cout << "Unknown Type. Supported Types are: ID, Trigger, Iso" << std::endl;
    nomval = -1;
    upval = -1;
    downval= -1;

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

  std::string IDinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root";
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "eleTrig_SF.root";
  std::string ISOinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "Isolation_SF.root";

  TFile *f_IDSF = new TFile(std::string(IDinputFile).c_str(),"READ");
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  TFile *f_ISOSF = new TFile(std::string(ISOinputFile).c_str(),"READ");

  h_ele_ID_abseta_pt_ratio = (TH2F*)f_IDSF->Get("EGamma_SF2D"); 
  h_ele_TRIGGER_abseta_pt_ratio = (TH2F*)f_TRIGGERSF->Get("h_eleTrig_SF");
  h_ele_ISO_abseta_pt_ratio = (TH2F*)f_ISOSF->Get("IsolationSF");

}

void LeptonSFHelper::SetMuonHistos( ){
  
  std::string IDinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "MuonID_Z_2016runB_2p6fb.root";
  std::string TRIGGERinputFile =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root";
  std::string ISOinputFile =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "MuonISO_Z_2016runB_2p6fb.root";

  TFile *f_IDSF = new TFile(std::string(IDinputFile).c_str(),"READ");
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  TFile *f_ISOSF = new TFile(std::string(ISOinputFile).c_str(),"READ");

  h_mu_ID_abseta_pt_ratio = (TH2F*)f_IDSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio"); 
  h_mu_TRIGGER_abseta_pt_ratio4p3 = (TH2F*)f_TRIGGERSF->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio");
  h_mu_TRIGGER_abseta_pt_ratio4p2 = (TH2F*)f_TRIGGERSF->Get("runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio");
  h_mu_ISO_abseta_pt_ratio = (TH2F*)f_ISOSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  
}

void LeptonSFHelper::SetElectronElectronHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "triggerSummary_ee.root";
  
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  
  h_ele_ele_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}

void LeptonSFHelper::SetMuonMuonHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "triggerSummary_mumu.root";
  
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  
  h_mu_mu_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}

void LeptonSFHelper::SetElectronMuonHistos( ){
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/" + "triggerSummary_emu.root";
  
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  
  h_ele_mu_TRIGGER_abseta_abseta = (TH2F*)f_TRIGGERSF->Get("scalefactor_eta2d_with_syst");
}
