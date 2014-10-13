#include "../interface/MiniAODHelper.h"

using namespace std;

// Constructor
MiniAODHelper::MiniAODHelper(){

  isSetUp = false;

  vertexIsSet = false;
  rhoIsSet = false;
  jetcorrectorIsSet = false;
  factorizedjetcorrectorIsSet = false;

  CSVLwp = 0.244;
  CSVMwp = 0.679;
  CSVTwp = 0.898;

  samplename = "blank";

}

// Destructor
MiniAODHelper::~MiniAODHelper(){

}

// Set up parameters one by one
void MiniAODHelper::SetUp(string iEra, int iSampleNumber, const analysisType::analysisType iAnalysis, bool iIsData){
  // Make sure we don't set up more than once
  if(isSetUp){ ThrowFatalError("Trying to set up 'BEANhelper' for the second time. Check your code."); }

  // Bring in the external values
  era          = iEra;
  sampleNumber = iSampleNumber;
  analysis     = iAnalysis;
  isData       = iIsData;

  // Error checking here
  if((era != "2011") && (era != "2012_52x") && (era != "2012_53x")){ ThrowFatalError("era set to '" + era + "' but it has to be either 2011, 2012_52x, or 2012_53x"); }
  if(sampleNumber==0){ ThrowFatalError("'sampleNumber' cannot be '0'."); }

  // Setup PU reweighing
  //SetUpPUreweighing(iCollisionDS);

  // Setup CSV reweighting
  //SetUpCSVreweighting();

  // Setup jet efficiency scale factors
  //SetUpJetSF();

  // Setup lepton efficiency scale factors
  //SetUpLeptonSF();

  // Awknowledge setup
  isSetUp = true;

}

// Set up parameters one by one
void MiniAODHelper::SetVertex(const reco::Vertex& inputVertex){

  vertex = inputVertex;

  vertexIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetRho(double inputRho){

  useRho = inputRho;

  rhoIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetJetCorrector(const JetCorrector* iCorrector){

  corrector = iCorrector;

  jetcorrectorIsSet = true;
}

// Set up parameters one by one
void MiniAODHelper::SetJetCorrectorUncertainty(){

  std::string inputJECfile = ( isData ) ? string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt" : string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt";

  jecUnc_ = new JetCorrectionUncertainty(inputJECfile);

}

// Set up parameters one by one
void MiniAODHelper::SetFactorizedJetCorrector(){

  // Create the JetCorrectorParameter objects, the order does not matter.
  //JetCorrectorParameters *ResJetPar = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/"); 
  JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L3Absolute_AK4PFchs.txt");
  JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L2Relative_AK4PFchs.txt");
  JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/PLS170_V7AN1_L1FastJet_AK4PFchs.txt");
  //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  std::vector<JetCorrectorParameters> vPar;
  vPar.push_back(*L1JetPar);
  vPar.push_back(*L2JetPar);
  vPar.push_back(*L3JetPar);
  //vPar->push_back(ResJetPar);

  useJetCorrector = new FactorizedJetCorrector(vPar);

  std::string inputJECfile = ( isData ) ? string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_DATA_Uncertainty_AK5PFchs.txt" : string(getenv("CMSSW_BASE")) + "/src/MiniAOD/MiniAODHelper/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt";

  jecUnc_ = new JetCorrectionUncertainty(inputJECfile);

  factorizedjetcorrectorIsSet = true;
}



std::vector<pat::Muon> 
MiniAODHelper::GetSelectedMuons(const std::vector<pat::Muon>& inputMuons, const float iMinPt, const muonID::muonID iMuonID){

  CheckSetUp();

  std::vector<pat::Muon> selectedMuons;

  for( std::vector<pat::Muon>::const_iterator it = inputMuons.begin(), ed = inputMuons.end(); it != ed; ++it ){
    if( isGoodMuon(*it,iMinPt,iMuonID) ) selectedMuons.push_back(*it);
  }

  return selectedMuons;
}


std::vector<pat::Electron> 
MiniAODHelper::GetSelectedElectrons(const std::vector<pat::Electron>& inputElectrons, const float iMinPt, const electronID::electronID iElectronID){

  CheckSetUp();

  std::vector<pat::Electron> selectedElectrons;

  for( std::vector<pat::Electron>::const_iterator it = inputElectrons.begin(), ed = inputElectrons.end(); it != ed; ++it ){
    if( isGoodElectron(*it,iMinPt,iElectronID) ) selectedElectrons.push_back(*it);
  }

  return selectedElectrons;
}


std::vector<pat::Tau> 
MiniAODHelper::GetSelectedTaus(const std::vector<pat::Tau>& inputTaus, const float iMinPt, const tauID::tauID iTauID){

  CheckSetUp();

  std::vector<pat::Tau> selectedTaus;

  for( std::vector<pat::Tau>::const_iterator it = inputTaus.begin(), ed = inputTaus.end(); it != ed; ++it ){
    if( isGoodTau(*it,iMinPt,iTauID) ) selectedTaus.push_back(*it);
  }

  return selectedTaus;
}



std::vector<pat::Jet> 
MiniAODHelper::GetSelectedJets(const std::vector<pat::Jet>& inputJets, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVwp){

  CheckSetUp();

  std::vector<pat::Jet> selectedJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    if( isGoodJet(*it, iMinPt, iMaxAbsEta, iJetID, iCSVwp) ) selectedJets.push_back(*it);
  }

  return selectedJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetUncorrectedJets(const std::vector<pat::Jet>& inputJets){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    jet.setP4( it->correctedJet(0).p4() );
    outputJets.push_back(jet);
  }

  return outputJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const edm::Event& event, const edm::EventSetup& setup, const sysType::sysType iSysType){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    double scale = 1.;

    if( jetcorrectorIsSet ) scale = corrector->correction(*it, event, setup);
    else std::cout << " !! ERROR !! Trying to use Full Framework GetCorrectedJets without setting jet corrector !" << std::endl;

    jet.scaleEnergy( scale );

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      jecUnc_->setJetEta(jet.eta());
      jecUnc_->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = jecUnc_->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = jecUnc_->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }

    outputJets.push_back(jet);
  }

  return outputJets;
}


std::vector<pat::Jet> 
MiniAODHelper::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, const sysType::sysType iSysType ){

  CheckSetUp();

  std::vector<pat::Jet> outputJets;

  if( !(factorizedjetcorrectorIsSet && rhoIsSet) ){
   std::cout << " !! ERROR !! Trying to use FWLite Framework GetCorrectedJets without setting factorized jet corrector !" << std::endl;

     return inputJets;
  }

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    pat::Jet jet = (*it);
    double scale = 1.;

    useJetCorrector->setJetEta(it->eta());
    useJetCorrector->setJetPt(it->pt());
    useJetCorrector->setJetA(it->jetArea());
    useJetCorrector->setRho(useRho); 

    scale = useJetCorrector->getCorrection();

    jet.scaleEnergy( scale );

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      jecUnc_->setJetEta(jet.eta());
      jecUnc_->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = jecUnc_->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = jecUnc_->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }

    outputJets.push_back(jet);
  }

  return outputJets;
}

bool 
MiniAODHelper::isGoodMuon(const pat::Muon& iMuon, const float iMinPt, const muonID::muonID iMuonID){

  CheckVertexSetUp();

  double minMuonPt = iMinPt;

  float maxLooseMuonAbsEta = 2.5;

  float maxTightMuonAbsEta = 2.1;


  // Be skeptical about this muon making it through
  bool passesKinematics	= false;
  bool passesIso        = false;
  bool passesID         = false;
  bool isPFMuon         = false;
  bool passesTrackerID  = false;

  bool passesGlobalTrackID   = false;
  bool passesMuonBestTrackID = false;
  bool passesInnerTrackID    = false;
  bool passesTrackID         = false;


  switch(iMuonID){
  case muonID::muonSide:
  case muonID::muonSideLooseMVA:
  case muonID::muonSideTightMVA:
  case muonID::muonPtOnly:
  case muonID::muonPtEtaOnly:
  case muonID::muonPtEtaIsoOnly:
  case muonID::muonPtEtaIsoTrackerOnly:
  case muonID::muonNoCuts:
  case muonID::muonLoose:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxLooseMuonAbsEta));
    passesIso        = (GetMuonRelIso(iMuon) < 0.200);
    isPFMuon         = true;
    passesID         = (( iMuon.isGlobalMuon() || iMuon.isTrackerMuon() ) && isPFMuon);
    break;
  case muonID::muonTight:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxTightMuonAbsEta));
    passesIso        = (GetMuonRelIso(iMuon) < 0.120);
    isPFMuon         = true;

    if( iMuon.globalTrack().isAvailable() ){
      passesGlobalTrackID = ( (iMuon.globalTrack()->normalizedChi2() < 10.) 
			      && (iMuon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0)
			      );
    }
    if( iMuon.muonBestTrack().isAvailable() ){
      passesMuonBestTrackID = ( (fabs(iMuon.muonBestTrack()->dxy(vertex.position())) < 0.2)
				&& (fabs(iMuon.muonBestTrack()->dz(vertex.position())) < 0.5)
				);
    }
    if( iMuon.innerTrack().isAvailable() )
      passesInnerTrackID = (iMuon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0);
    if( iMuon.track().isAvailable() )
      passesTrackID = (iMuon.track()->hitPattern().trackerLayersWithMeasurement() > 5);

    passesTrackerID = ( passesGlobalTrackID && passesMuonBestTrackID && passesInnerTrackID && passesTrackID && (iMuon.numberOfMatchedStations() > 1) );

    passesID        = ((iMuon.isGlobalMuon() || iMuon.isTrackerMuon()) && isPFMuon && passesTrackerID);
    break;
  }

  return (passesKinematics && passesIso && passesID);
}



bool 
MiniAODHelper::isGoodElectron(const pat::Electron& iElectron, const float iMinPt, const electronID::electronID iElectronID){

  CheckVertexSetUp();

  double minElectronPt = iMinPt;

  float maxLooseElectronAbsEta = 2.5;
  float maxTightElectronAbsEta = 2.5;


  // Be skeptical about this electron making it through
  bool passesKinematics	= false;
  bool passesIso        = false;
  bool passesID         = false;

  bool inCrack = false;
  if( iElectron.superCluster().isAvailable() )
    inCrack = ( fabs(iElectron.superCluster()->position().eta())>1.4442 && fabs(iElectron.superCluster()->position().eta())<1.5660 );


  bool myTrigPresel = true;

  double eleID      = iElectron.electronID("eidTight");
  bool passMVAId53x = ( eleID>0.5 );  // For 2012_53x, tighter selection

  bool d02 = false; 
  bool d04 = false;
  bool dZ  = false;
  bool no_exp_inner_trkr_hits = false;
  if( iElectron.gsfTrack().isAvailable() ){
    d02 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.02 );
    d04 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.04 );
    no_exp_inner_trkr_hits = ( iElectron.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 );
    dZ = ( fabs(iElectron.gsfTrack()->dz(vertex.position())) < 1. );
  }


  bool notConv = ( iElectron.passConversionVeto() );
  bool id      = ( passMVAId53x && d02 && dZ && notConv );


  switch(iElectronID){
  case electronID::electronSide:
  case electronID::electronSideLooseMVA:
  case electronID::electronSideTightMVA:
  case electronID::electronLooseMinusTrigPresel:
  case electronID::electronNoCuts:
  case electronID::electronLoose:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxLooseElectronAbsEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.200);
    passesID         = ( passMVAId53x && no_exp_inner_trkr_hits && d04 && notConv && myTrigPresel );
    break;
  case electronID::electronTightMinusTrigPresel:
  case electronID::electronTight:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxTightElectronAbsEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.100);
    passesID         = ( id && no_exp_inner_trkr_hits && myTrigPresel );
    break;
  }

  return (passesKinematics && passesIso && passesID);
}

bool 
MiniAODHelper::isGoodTau(const pat::Tau& iTau, const float iMinPt, const tauID::tauID iTauID){

  CheckVertexSetUp();
 
  double minTauPt = iMinPt;
  
  bool passesKinematics = false;
  passesKinematics = (iTau.pt() >= 20) && (fabs(iTau.eta()) <= 2.1) && (iTau.pt() > minTauPt);
  return passesKinematics;
}

bool 
MiniAODHelper::isGoodJet(const pat::Jet& iJet, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVworkingPoint){

  CheckVertexSetUp();

  // Transverse momentum requirement
  if( iJet.pt() < iMinPt ) return false;

  // Absolute eta requirement
  if( fabs(iJet.eta()) > iMaxAbsEta ) return false;

  bool loose = (
		iJet.neutralHadronEnergyFraction() < 0.99 &&
		iJet.chargedEmEnergyFraction() < 0.99 &&
		iJet.neutralEmEnergyFraction() < 0.99 &&
		iJet.numberOfDaughters() > 1
		);

  if( fabs(iJet.eta())<2.4 ){
    loose = ( loose &&
	      iJet.chargedHadronEnergyFraction() > 0.0 &&
	      iJet.chargedMultiplicity() > 0
	      );
  }

  // Jet ID
  switch(iJetID){
  case jetID::none:
  case jetID::jetMinimal:
  case jetID::jetLooseAOD:
  case jetID::jetLoose:
  case jetID::jetTight:
    if( !loose ) return false;
    break;
  default:
    break;
  }

  if( !PassesCSV(iJet, iCSVworkingPoint) ) return false;

  return true;
}


float MiniAODHelper::GetMuonRelIso(const pat::Muon& iMuon) const
{
  float result = 9999; 

  double pfIsoCharged = iMuon.pfIsolationR04().sumChargedHadronPt;
  double pfIsoNeutral = iMuon.pfIsolationR04().sumNeutralHadronEt + iMuon.pfIsolationR04().sumPhotonEt;

  double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - 0.5*iMuon.pfIsolationR04().sumPUPt );

  result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
  
  return result;
}


float MiniAODHelper::GetElectronRelIso(const pat::Electron& iElectron) const
{
  float result = 9999; 

  double pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
  double pfIsoNeutral = iElectron.pfIsolationVariables().sumNeutralHadronEt + iElectron.pfIsolationVariables().sumPhotonEt;

  double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - 0.5*iElectron.pfIsolationVariables().sumPUPt );

  result = (pfIsoCharged + pfIsoPUSubtracted)/iElectron.pt();
  
  return result;
}

bool MiniAODHelper::PassesCSV(const pat::Jet& iJet, const char iCSVworkingPoint){
  CheckSetUp();

  float csvValue = iJet.bDiscriminator("combinedSecondaryVertexBJetTags");

  // CSV b-tagging requirement
  switch(iCSVworkingPoint){
  case 'L':	if(csvValue > CSVLwp){ return true; }	break;
  case 'M':	if(csvValue > CSVMwp){ return true; }	break;
  case 'T':	if(csvValue > CSVTwp){ return true; }	break;
  case '-':	return true;                            break;
  }
  return false;
}
