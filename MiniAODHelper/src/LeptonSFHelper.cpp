#include "MiniAOD/MiniAODHelper/interface/LeptonSFHelper.h"

//PUBLIC
LeptonSFHelper::LeptonSFHelper( ){

  //std::cout << "Initializing Lepton scale factors" << std::endl;

  SetElectronHistos( );
  SetMuonHistos( );
  SetElectronElectronHistos( );
  SetMuonMuonHistos( );
  SetElectronMuonHistos( );

  electronMaxPt = 199.0;
  electronMaxPtHigh=250.0;
  muonMaxPt = 119.0;
  muonMaxPtHigh = 499.0;
  
  ljets_mu_BtoF_lumi=19691.782;
  ljets_mu_GtoH_lumi=16226.452;


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
  float ElectronGFSSF = 1.0;
  float ElectronGFSSF_Up = 1.0;
  float ElectronGFSSF_Down = 1.0;
  float MuonIDSF = 1.0;
  float MuonIDSF_Up = 1.0;
  float MuonIDSF_Down = 1.0;

  float MuonHIPSF = 1.0;
  float MuonHIPSF_Up = 1.0;
  float MuonHIPSF_Down = 1.0;

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

    ElectronIDSF = ElectronIDSF * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 0, "ID");
    ElectronIDSF_Up = ElectronIDSF_Up *GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 1, "ID");
    ElectronIDSF_Down = ElectronIDSF_Down * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), -1, "ID");

    ElectronIsoSF = ElectronIsoSF * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 0, "Iso");
    ElectronIsoSF_Up = ElectronIsoSF_Up  * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 1, "Iso");
    ElectronIsoSF_Down = ElectronIsoSF_Down * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), -1, "Iso");

    ElectronTriggerSF = ElectronTriggerSF * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 0, "Trigger");
    ElectronTriggerSF_Up = ElectronTriggerSF_Up  * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 1, "Trigger");
    ElectronTriggerSF_Down = ElectronTriggerSF_Down * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), -1, "Trigger");

    ElectronGFSSF = ElectronGFSSF * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 0, "GFS");
    ElectronGFSSF_Up = ElectronGFSSF_Up *GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), 1, "GFS");
    ElectronGFSSF_Down = ElectronGFSSF_Down * GetElectronSF(Electron.pt(), Electron.superCluster()->eta(), -1, "GFS");


  } for (auto Muon: Muons){ //Muon is of type pat::Muon

    MuonIDSF = MuonIDSF * GetMuonSF(Muon.pt(), Muon.eta(), 0, "ID");
    MuonIDSF_Up = MuonIDSF_Up * GetMuonSF(Muon.pt(), Muon.eta(), 1, "ID");
    MuonIDSF_Down = MuonIDSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1, "ID");

    MuonHIPSF = MuonHIPSF * GetMuonSF(Muon.pt(), Muon.eta(), 0, "HIP");
    MuonHIPSF_Up = MuonHIPSF_Up * GetMuonSF(Muon.pt(), Muon.eta(), 1, "HIP");
    MuonHIPSF_Down = MuonHIPSF_Down * GetMuonSF(Muon.pt(), Muon.eta(), -1, "HIP");

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
  ScaleFactorMap["ElectronSFGFS"] = ElectronGFSSF;
  ScaleFactorMap["ElectronSFGFS_Up"] = ElectronGFSSF_Up;
  ScaleFactorMap["ElectronSFGFS_Down"] = ElectronGFSSF_Down;
  ScaleFactorMap["ElectronElectronTriggerSF"] =ElectronElectronTriggerSF;

  ScaleFactorMap["MuonSFID"] = MuonIDSF;
  ScaleFactorMap["MuonSFID_Up"] = MuonIDSF_Up;
  ScaleFactorMap["MuonSFID_Down"] = MuonIDSF_Down;

  ScaleFactorMap["MuonSFHIP"] = MuonHIPSF;
  ScaleFactorMap["MuonSFHIP_Up"] = MuonHIPSF_Up;
  ScaleFactorMap["MuonSFHIP_Down"] = MuonHIPSF_Down;

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
  if (type=="Trigger"){
    searchPt=TMath::Min( electronPt , electronMaxPtHigh );
  }

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

    thisBin = h_ele_TRIGGER_abseta_pt_ratio->FindBin( searchPt, searchEta );
    nomval=h_ele_TRIGGER_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_TRIGGER_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error;
    downval=nomval-error;

  }
  else if ( type == "Iso" ){

    thisBin = h_ele_ISO_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_ele_ISO_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_ISO_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error;  //DANGERZONE need to add pT depnednet 1% uncertainty
    downval=nomval-error;

  }
  else if ( type == "GFS" ){

    thisBin = h_ele_GFS_abseta_pt_ratio->FindBin( searchEta , searchPt );
    nomval=h_ele_GFS_abseta_pt_ratio->GetBinContent( thisBin );
    error=h_ele_GFS_abseta_pt_ratio->GetBinError( thisBin );
    upval=nomval+error; //DANGERZONE need to add pT depnednet 1% uncertainty
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
  if (type=="Trigger"){
    searchPt=TMath::Min( muonPt , muonMaxPtHigh );
  }
  float nomval = 0;
  float error = 0;
  float upval = 0;
  float downval= 0;
  float nomvalBtoF = 0;
  float errorBtoF = 0;
  float upvalBtoF = 0;
  float downvalBtoF= 0;
  float nomvalGtoH = 0;
  float errorGtoH = 0;
  float upvalGtoH = 0;
  float downvalGtoH= 0;
  

  if ( type == "ID" ){

    thisBin = h_mu_ID_abseta_pt_ratioBtoF->FindBin(  searchPt, searchEta  );
    nomvalBtoF=h_mu_ID_abseta_pt_ratioBtoF->GetBinContent( thisBin );
    errorBtoF=h_mu_ID_abseta_pt_ratioBtoF->GetBinError( thisBin );
    upvalBtoF=( nomvalBtoF+errorBtoF );
    downvalBtoF=( nomvalBtoF-errorBtoF );
    upvalBtoF=upvalBtoF*( 1.0+sqrt(0.01*0.01+0.005*0.005) );
    downvalBtoF=downvalBtoF*( 1.0-sqrt(0.01*0.01+0.005*0.005) );
    
    thisBin = h_mu_ID_abseta_pt_ratioGtoH->FindBin(  searchPt, searchEta  );
    nomvalGtoH=h_mu_ID_abseta_pt_ratioGtoH->GetBinContent( thisBin );
    errorGtoH=h_mu_ID_abseta_pt_ratioGtoH->GetBinError( thisBin );
    upvalGtoH=( nomvalGtoH+errorGtoH );
    downvalGtoH=( nomvalGtoH-errorGtoH );
    upvalGtoH=upvalGtoH*( 1.0+sqrt(0.01*0.01+0.005*0.005) );
    downvalGtoH=downvalGtoH*( 1.0-sqrt(0.01*0.01+0.005*0.005) );

    nomval=(ljets_mu_BtoF_lumi*nomvalBtoF + ljets_mu_GtoH_lumi * nomvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    upval=(ljets_mu_BtoF_lumi*upvalBtoF + ljets_mu_GtoH_lumi * upvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    downval=(ljets_mu_BtoF_lumi*downvalBtoF + ljets_mu_GtoH_lumi * downvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);

  }
  else if ( type == "Trigger" ){

    thisBin = h_mu_TRIGGER_abseta_ptBtoF->FindBin(  searchPt, searchEta  );
    nomvalBtoF=h_mu_TRIGGER_abseta_ptBtoF->GetBinContent( thisBin );
    errorBtoF=h_mu_TRIGGER_abseta_ptBtoF->GetBinError( thisBin );
    upvalBtoF=( nomvalBtoF+errorBtoF );
    downvalBtoF=( nomvalBtoF-errorBtoF );
    
    thisBin = h_mu_TRIGGER_abseta_ptGtoH->FindBin(  searchPt, searchEta  );
    nomvalGtoH=h_mu_TRIGGER_abseta_ptGtoH->GetBinContent( thisBin );
    errorGtoH=h_mu_TRIGGER_abseta_ptGtoH->GetBinError( thisBin );
    upvalGtoH=( nomvalGtoH+errorGtoH );
    downvalGtoH=( nomvalGtoH-errorGtoH );

    nomval=(ljets_mu_BtoF_lumi*nomvalBtoF + ljets_mu_GtoH_lumi * nomvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    upval=(ljets_mu_BtoF_lumi*upvalBtoF + ljets_mu_GtoH_lumi * upvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    downval=(ljets_mu_BtoF_lumi*downvalBtoF + ljets_mu_GtoH_lumi * downvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    
  }
  else if ( type == "Iso" ){
    
    
    thisBin = h_mu_ISO_abseta_pt_ratioBtoF->FindBin(  searchPt, searchEta  );
    nomvalBtoF=h_mu_ISO_abseta_pt_ratioBtoF->GetBinContent( thisBin );
    errorBtoF=h_mu_ISO_abseta_pt_ratioBtoF->GetBinError( thisBin );
    upvalBtoF=( nomvalBtoF+errorBtoF );
    downvalBtoF=( nomvalBtoF-errorBtoF );
    upvalBtoF=upvalBtoF*(1.0+0.005  );
    downvalBtoF=downvalBtoF*(1.0-0.005 );
    
    thisBin = h_mu_ISO_abseta_pt_ratioGtoH->FindBin(  searchPt, searchEta  );
    nomvalGtoH=h_mu_ISO_abseta_pt_ratioGtoH->GetBinContent( thisBin );
    errorGtoH=h_mu_ISO_abseta_pt_ratioGtoH->GetBinError( thisBin );
    upvalGtoH=( nomvalGtoH+errorGtoH );
    downvalGtoH=( nomvalGtoH-errorGtoH );
    upvalGtoH=upvalGtoH*(1.0+0.005  );
    downvalGtoH=downvalGtoH*( 1.0-0.005 );

    nomval=(ljets_mu_BtoF_lumi*nomvalBtoF + ljets_mu_GtoH_lumi * nomvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    upval=(ljets_mu_BtoF_lumi*upvalBtoF + ljets_mu_GtoH_lumi * upvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    downval=(ljets_mu_BtoF_lumi*downvalBtoF + ljets_mu_GtoH_lumi * downvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);

  }

  else if ( type == "HIP" ){

    thisBin = h_mu_HIP_eta_ratioBtoF->FindBin( searchEta );
    nomvalBtoF=h_mu_HIP_eta_ratioBtoF->GetBinContent( thisBin );
    errorBtoF=h_mu_HIP_eta_ratioBtoF->GetBinError( thisBin );
    upvalBtoF=( nomvalBtoF+errorBtoF );
    downvalBtoF=( nomvalBtoF-errorBtoF );
    
    thisBin = h_mu_HIP_eta_ratioGtoH->FindBin( searchEta );
    nomvalGtoH=h_mu_HIP_eta_ratioGtoH->GetBinContent( thisBin );
    errorGtoH=h_mu_HIP_eta_ratioGtoH->GetBinError( thisBin );
    upvalGtoH=( nomvalGtoH+errorGtoH );
    downvalGtoH=( nomvalGtoH-errorGtoH );
    
    nomval=(ljets_mu_BtoF_lumi*nomvalBtoF + ljets_mu_GtoH_lumi * nomvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    upval=(ljets_mu_BtoF_lumi*upvalBtoF + ljets_mu_GtoH_lumi * upvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
    downval=(ljets_mu_BtoF_lumi*downvalBtoF + ljets_mu_GtoH_lumi * downvalGtoH)/(ljets_mu_BtoF_lumi+ljets_mu_GtoH_lumi);
   
//     upval=upval*( 1.0+0.005 );
//     downval=downval*( 1.0-0.005 );


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

  std::string IDinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "ele_ID_SF.root";
  std::string TRIGGERinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "ele_TriggerSF_Run2016All_v1.root";
  std::string ISOinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "ele_Reco_EGM2D.root"; // DANGERZONE: no iso SF yet??
  std::string GFSinputFile = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "ele_Reco_EGM2D.root";

  TFile *f_IDSF = new TFile(std::string(IDinputFile).c_str(),"READ");
  TFile *f_TRIGGERSF = new TFile(std::string(TRIGGERinputFile).c_str(),"READ");
  TFile *f_ISOSF = new TFile(std::string(ISOinputFile).c_str(),"READ");
  TFile *f_GFSSF = new TFile(std::string(GFSinputFile).c_str(),"READ");

  h_ele_ID_abseta_pt_ratio = (TH2F*)f_IDSF->Get("EGamma_SF2D");
  h_ele_TRIGGER_abseta_pt_ratio = (TH2F*)f_TRIGGERSF->Get("Ele27_WPTight_Gsf");
  h_ele_ISO_abseta_pt_ratio = (TH2F*)f_ISOSF->Get("EGamma_SF2D");
  h_ele_GFS_abseta_pt_ratio = (TH2F*)f_GFSSF->Get("EGamma_SF2D");

}

void LeptonSFHelper::SetMuonHistos( ){

  std::string IDinputFileBtoF = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_ID_EfficienciesAndSF_BCDEF.root";
  std::string IDinputFileGtoH = std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_ID_EfficienciesAndSF_GH.root";

  std::string TRIGGERinputFileBtoF =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_TRIGGER_BtoF.root";
  std::string TRIGGERinputFileGtoH =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_TRIGGER_GtoH.root";

  std::string ISOinputFileBtoF =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_ISO_EfficienciesAndSF_BCDEF.root";
  std::string ISOinputFileGtoH =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "mu_ISO_EfficienciesAndSF_GH.root";
  
  std::string HIPinputFileBtoF =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "HIP_BCDEF_histos.root";
  std::string HIPinputFileGtoH =  std::string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/leptonSF/feb160317/" + "HIP_GH_histos.root";


  TFile *f_IDSFBtoF = new TFile(std::string(IDinputFileBtoF).c_str(),"READ");
  TFile *f_IDSFGtoH = new TFile(std::string(IDinputFileGtoH).c_str(),"READ");
  
  TFile *f_HIPSFBtoF = new TFile(std::string(HIPinputFileBtoF).c_str(),"READ");
  TFile *f_HIPSFGtoH = new TFile(std::string(HIPinputFileGtoH).c_str(),"READ");

  
  TFile *f_TRIGGERSFBtoF = new TFile(std::string(TRIGGERinputFileBtoF).c_str(),"READ");
  TFile *f_TRIGGERSFGtoH = new TFile(std::string(TRIGGERinputFileGtoH).c_str(),"READ");

  TFile *f_ISOSFBtoF = new TFile(std::string(ISOinputFileBtoF).c_str(),"READ");
  TFile *f_ISOSFGtoH = new TFile(std::string(ISOinputFileGtoH).c_str(),"READ");

  h_mu_ID_abseta_pt_ratioBtoF = (TH2F*)f_IDSFBtoF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
  h_mu_ID_abseta_pt_ratioGtoH = (TH2F*)f_IDSFGtoH->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");

  h_mu_HIP_eta_ratioBtoF = (TH1D*)f_HIPSFBtoF->Get("ratio_eff_aeta_dr030e030_corr");
  h_mu_HIP_eta_ratioGtoH = (TH1D*)f_HIPSFGtoH->Get("ratio_eff_aeta_dr030e030_corr");

  h_mu_TRIGGER_abseta_ptBtoF= (TH2F*)f_TRIGGERSFBtoF->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");
  h_mu_TRIGGER_abseta_ptGtoH= (TH2F*)f_TRIGGERSFGtoH->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");

  h_mu_ISO_abseta_pt_ratioBtoF = (TH2F*)f_ISOSFBtoF->Get("TightISO_TightID_pt_eta/pt_abseta_ratio");
  h_mu_ISO_abseta_pt_ratioGtoH = (TH2F*)f_ISOSFGtoH->Get("TightISO_TightID_pt_eta/pt_abseta_ratio");

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
