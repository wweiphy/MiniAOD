#ifndef _LeptonSFHelper_h
#define _LeptonSFHelper_h

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <map>

#include "TMath.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class LeptonSFHelper {

 public:
  LeptonSFHelper(const edm::ParameterSet& iConfig);
  ~LeptonSFHelper( );

  std::map< std::string, float>  GetLeptonSF( const std::vector< pat::Electron >& Electrons,
					      const std::vector< pat::Muon >& Muons);

  float GetLeptonTriggerSF(  const double& lepton_pt , const double& lepton_eta , const int& syst, TH2F* h_SFs);
  float GetElectronElectronSF( float electronEta1, float electronEta2, int syst , std::string type);
  float GetMuonMuonSF( float muonEta1, float muonEta2, int syst , std::string type);
  float GetElectronMuonSF( float electronEta, float muonEta, int syst , std::string type);

 private:

  void SetElectronHistos( );
  void SetMuonHistos( );
  void SetElectronElectronHistos( );
  void SetMuonMuonHistos( );
  void SetElectronMuonHistos( );
  int findPoint(TGraphAsymmErrors& graph,float& x_);
  float getValue(TGraphAsymmErrors& graph,float& x_,int syst);

  TH2F *h_ele_TRIGGER_abseta_pt_ratio;

  TH2F *h_mu_TRIGGER_abseta_pt;
  
  TH2F *h_ele_ele_TRIGGER_abseta_abseta;
  TH2F *h_mu_mu_TRIGGER_abseta_abseta;
  TH2F *h_ele_mu_TRIGGER_abseta_abseta;

  TFile *f_muon_TRIGGERSF = nullptr;
  std::string muon_TRIGGERinputFile;
  std::string muon_TRIGGERhistname;

  TFile *f_electron_TRIGGERSF = nullptr;

  std::string electron_TRIGGERinputFile;
  std::string electron_TRIGGERhistname;

};

#endif
