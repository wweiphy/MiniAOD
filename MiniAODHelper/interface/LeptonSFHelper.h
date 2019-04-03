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


class LeptonSFHelper {

 public:
  LeptonSFHelper( );
  ~LeptonSFHelper( );

  std::map< std::string, float>  GetLeptonSF( const std::vector< pat::Electron >& Electrons,
					      const std::vector< pat::Muon >& Muons);

  float GetLeptonTriggerSF(  const double& lepton_pt , const double& lepton_eta , const int& syst, TH2F* h_SFs);
  float GetElectronSF(  float electronPt , float electronEta , int syst);
  float GetMuonSF( const double& muonPt , const double& muonEta , const int& syst);
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

  float electronLowPtRangeCut;
  float electronMaxPt;
  float electronMinPt;
  float electronMinPtLowPt;
  float electronMaxPtLowPt;
  
  float electronMaxPtHigh;
  float electronMaxPtHigher;
  float electronMaxEta;
  float electronMaxEtaLow;
  
  
  float muonMaxPt;
  float muonMaxPtHigh;
  float muonMinPt;
  float muonMinPtHigh;
  float muonMaxEta;

  TFile *f_muon_TRIGGERSF = nullptr;
  std::string muon_TRIGGERinputFile;
  std::string muon_TRIGGERhistname;

  TFile *f_electron_TRIGGERSF = nullptr;

  std::string electron_TRIGGERinputFile;
  std::string electron_TRIGGERhistname;
  
  // float ljets_mu_BtoF_lumi;
  // //float ljets_mu_GtoH_lumi;
  // float ljets_ele_BtoF_lumi;
  //float ljets_ele_GtoH_lumi;  

};

#endif
