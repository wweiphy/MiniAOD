#include "../interface/MiniAODHelper.h"

using namespace std;

// Constructor
MiniAODHelper::MiniAODHelper(){

  isSetUp = false;

  vertexIsSet = false;
  rhoIsSet = false;
  jetcorrectorIsSet = false;
  factorizedjetcorrectorIsSet = false;
  
  // twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging#Preliminary_working_or_operating
  // Preliminary working (or operating) points for CSVv2+IVF
  CSVLwp = 0.605;//CSVv2 0.423; // 10.1716% DUSG mistag efficiency
  CSVMwp = 0.890;//CSVv2 0.814; // 1.0623% DUSG mistag efficiency
  CSVTwp = 0.970;//CSVv2 0.941; // 0.1144% DUSG mistag efficiency

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
  if((era != "2011") && (era != "2012_52x") && (era != "2012_53x") && (era != "2015_72x") && (era != "2015_73x") && (era != "2015_74x")){ ThrowFatalError("era set to '" + era + "' but it has to be either 2011, 2012_52x, 2012_53x, 2015_72x, 2015_73x, or 2015_74x"); }
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

namespace {
  struct ByEta {
    bool operator()(const pat::PackedCandidate *c1, const pat::PackedCandidate *c2) const {
      return c1->eta() < c2->eta();
    }
    bool operator()(float c1eta, const pat::PackedCandidate *c2) const {
      return c1eta < c2->eta();
    }
    bool operator()(const pat::PackedCandidate *c1, float c2eta) const {
      return c1->eta() < c2eta;
    }
  };
}

// Set up parameters one by one
void MiniAODHelper::SetPackedCandidates(const std::vector<pat::PackedCandidate> & all, int fromPV_thresh, float dz_thresh, bool also_leptons){

  allcands_ = &all;
  charged_.clear(); neutral_.clear(); pileup_.clear();
  for (const pat::PackedCandidate &p : all) {
    if (p.charge() == 0) {
      neutral_.push_back(&p);
    } 
    else { 

      if ( (abs(p.pdgId()) == 211 ) || ( also_leptons && ((abs(p.pdgId()) == 11 ) || (abs(p.pdgId()) == 13 )) ) )  {
	
	if (p.fromPV() > fromPV_thresh && fabs(p.dz()) < dz_thresh ) {
	  charged_.push_back(&p);
	} 
	else {
	  pileup_.push_back(&p);
	}
      }
    }
  }
  //  if (weightCone_ > 0) weights_.resize(neutral_.size());
  //  std::fill(weights_.begin(), weights_.end(), -1.f);
  std::sort(charged_.begin(), charged_.end(), ByEta());
  std::sort(neutral_.begin(), neutral_.end(), ByEta());
  std::sort(pileup_.begin(),  pileup_.end(),  ByEta());
  clearVetos();
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
MiniAODHelper::GetSelectedMuons(const std::vector<pat::Muon>& inputMuons, const float iMinPt, const muonID::muonID iMuonID, const coneSize::coneSize iconeSize, const corrType::corrType icorrType, const float iMaxEta){

  CheckSetUp();

  std::vector<pat::Muon> selectedMuons;

  for( std::vector<pat::Muon>::const_iterator it = inputMuons.begin(), ed = inputMuons.end(); it != ed; ++it ){
    if( isGoodMuon(*it,iMinPt,iMaxEta,iMuonID,iconeSize,icorrType) ) selectedMuons.push_back(*it);
  }

  return selectedMuons;
}


std::vector<pat::Electron> 
MiniAODHelper::GetSelectedElectrons(const std::vector<pat::Electron>& inputElectrons, const float iMinPt, const electronID::electronID iElectronID, const float iMaxEta){

  CheckSetUp();

  std::vector<pat::Electron> selectedElectrons;

  for( std::vector<pat::Electron>::const_iterator it = inputElectrons.begin(), ed = inputElectrons.end(); it != ed; ++it ){
    if( isGoodElectron(*it,iMinPt,iMaxEta,iElectronID) ) selectedElectrons.push_back(*it);
  }

  return selectedElectrons;
}

std::vector<pat::Tau> 
MiniAODHelper::GetSelectedTaus(const std::vector<pat::Tau>& inputTaus, const float iMinPt, const tau::ID id){

  CheckSetUp();

  std::vector<pat::Tau> selectedTaus;

  for( std::vector<pat::Tau>::const_iterator it = inputTaus.begin(), ed = inputTaus.end(); it != ed; ++it ){
    if( isGoodTau(*it,iMinPt,id) ) selectedTaus.push_back(*it);
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


std::vector<pat::Jet> MiniAODHelper::GetUncorrectedJets(
	const std::vector<pat::Jet> &inputJets)
{
	CheckSetUp();
	
	std::vector<pat::Jet> outputJets;
	outputJets.reserve(inputJets.size());
	
	for (std::vector<pat::Jet>::const_iterator it = inputJets.begin(),
		ed = inputJets.end(); it != ed; ++it) {
		
		pat::Jet jet = (*it);
		jet.setP4(it->correctedJet(0).p4());
		outputJets.push_back(jet);
	}
	
	return outputJets;
}


/// Overload of GetUncorrectedJets(const std::vector<pat::Jet>&)
std::vector<pat::Jet> MiniAODHelper::GetUncorrectedJets(
	edm::Handle<pat::JetCollection> inputJets)
{
	CheckSetUp();
	
	std::vector<pat::Jet> outputJets;
	outputJets.reserve(inputJets->size());
	
	for (pat::JetCollection::const_iterator it = inputJets->begin(),
		ed = inputJets->end(); it != ed; ++it) {
		
		pat::Jet jet = (*it);
		jet.setP4(it->correctedJet(0).p4());
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

    /// JER
    double jerSF = 1.;
    if( jet.genJet() ){
      if( iSysType == sysType::JERup ){
	jerSF = getJERfactor(1, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
      }
      else if( iSysType == sysType::JERdown ){
	jerSF = getJERfactor(-1, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
      }
      else {
	jerSF = getJERfactor(0, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
      }
      // std::cout << "----->checking gen Jet pt " << jet.genJet()->pt() << ",  jerSF is" << jerSF << std::endl;
    }
    // else     std::cout << "    ==> can't find genJet" << std::endl;

    jet.scaleEnergy( jerSF );


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


std::vector<boosted::HTTTopJet> 
MiniAODHelper::GetSelectedTopJets(const std::vector<boosted::HTTTopJet>& inputJets, const float iMinFatPt, const float iMaxAbsFatEta, const float iMinSubPt, const float iMaxAbsSubEta, const jetID::jetID iJetID){

  CheckSetUp();

  std::vector<boosted::HTTTopJet> selectedJets;

  for( std::vector<boosted::HTTTopJet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    if( isGoodTopJet(*it, iMinFatPt, iMaxAbsFatEta, iMinSubPt, iMaxAbsSubEta, iJetID) ) selectedJets.push_back(*it);
  }

  return selectedJets;
}


std::vector<boosted::SubFilterJet> 
MiniAODHelper::GetSelectedHiggsJets(const std::vector<boosted::SubFilterJet>& inputJets, const float iMinFatPt, const float iMaxAbsFatEta, const float iMinSubPt, const float iMaxAbsSubEta, const jetID::jetID iJetID){

  CheckSetUp();

  std::vector<boosted::SubFilterJet> selectedJets;

  for( std::vector<boosted::SubFilterJet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    
    boosted::SubFilterJet higgsJet = *it;
    std::vector<pat::Jet> filterjets;
    
    for( std::vector<pat::Jet>::const_iterator itFilt = higgsJet.filterjets.begin(), edFilt = higgsJet.filterjets.end(); itFilt != edFilt; ++itFilt ){
      if( isGoodJet(*itFilt, iMinSubPt, iMaxAbsSubEta, iJetID, '-') ) filterjets.push_back(*itFilt);
    }
    
    higgsJet.filterjets = filterjets;
      
    if( isGoodHiggsJet(higgsJet, iMinFatPt, iMaxAbsFatEta) ) selectedJets.push_back(higgsJet);
  }

  return selectedJets;
}


bool 
MiniAODHelper::isGoodMuon(const pat::Muon& iMuon, const float iMinPt, const float iMaxEta, const muonID::muonID iMuonID, const coneSize::coneSize iconeSize, const corrType::corrType icorrType){

  CheckVertexSetUp();

  double minMuonPt = iMinPt;

  double maxMuonEta = iMaxEta;


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
  case muonID::muonPreselection:
  case muonID::muonSide:
  case muonID::muonSideLooseMVA:
  case muonID::muonSideTightMVA:
  case muonID::muonPtOnly:
  case muonID::muonPtEtaOnly:
  case muonID::muonPtEtaIsoOnly:
  case muonID::muonPtEtaIsoTrackerOnly:
  case muonID::muonRaw:
  case muonID::muonLooseMvaBased:
  case muonID::muonTightMvaBased:
  case muonID::muonLooseCutBased:
  case muonID::muonTightCutBased:
  case muonID::muonCutBased:
  case muonID::muon2lss:
  case muonID::muonLoose:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxMuonEta));
    passesIso        = (GetMuonRelIso(iMuon,iconeSize,icorrType) < 0.200);
    isPFMuon         = iMuon.isPFMuon();
    
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

    passesID        = (iMuon.isGlobalMuon() && isPFMuon && passesTrackerID);
    
    break;
  case muonID::muonTight:
    passesKinematics = ((iMuon.pt() >= minMuonPt) && (fabs(iMuon.eta()) <= maxMuonEta));
    passesIso        = (GetMuonRelIso(iMuon,iconeSize,icorrType) < 0.12);
    isPFMuon         = iMuon.isPFMuon();

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

    passesID        = (iMuon.isGlobalMuon() && isPFMuon && passesTrackerID);
    break;
  }

  return (passesKinematics && passesIso && passesID);
}



bool 
MiniAODHelper::isGoodElectron(const pat::Electron& iElectron, const float iMinPt, const float iMaxEta,const electronID::electronID iElectronID){

  CheckVertexSetUp();

  double minElectronPt = iMinPt;

  float maxElectronEta = iMaxEta;


  // Be skeptical about this electron making it through
  bool passesKinematics	= false;
  bool passesIso        = false;
  bool passesID         = false;

  double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
  double absSCeta = fabs(SCeta);

  bool inCrack = false;
  if( iElectron.superCluster().isAvailable() ) inCrack = ( absSCeta>1.4442 && absSCeta<1.5660 );


  bool myTrigPresel = true;

  //double eleID      = iElectron.electronID("eidTight");
  bool passMVAId53x = true;//( eleID>0.5 );  // For 2012_53x, tighter selection

  bool d02 = false; 
  bool d04 = false;
  bool dZ  = false;
  bool no_exp_inner_trkr_hits = true; //false; // see below
  if( iElectron.gsfTrack().isAvailable() ){
    d02 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.02 );
    d04 = ( fabs(iElectron.gsfTrack()->dxy(vertex.position())) < 0.04 );
    //no_exp_inner_trkr_hits = ( iElectron.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 ); // deprecated in 7_2_0 .. replace with ..?
    dZ = ( fabs(iElectron.gsfTrack()->dz(vertex.position())) < 1. );
  }


  bool notConv = ( iElectron.passConversionVeto() );
  bool id      = ( passMVAId53x && d02 && dZ && notConv );


  switch(iElectronID){
  case electronID::electronPreselection:
  case electronID::electronSide:
  case electronID::electronSideLooseMVA:
  case electronID::electronSideTightMVA:
  case electronID::electronLooseMinusTrigPresel:
  case electronID::electronRaw:
  case electronID::electronLooseCutBased:
  case electronID::electronTightCutBased:
  case electronID::electronCutBased:
  case electronID::electronLooseMvaBased:
  case electronID::electronTightMvaBased:
  case electronID::electron2lss:
  case electronID::electronLoose:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxElectronEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.200);
    passesID         = ( passMVAId53x && no_exp_inner_trkr_hits && d04 && notConv && myTrigPresel );
    break;
  case electronID::electronTightMinusTrigPresel:
  case electronID::electronTight:
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxElectronEta) && !inCrack);
    passesIso        = (GetElectronRelIso(iElectron) < 0.100);
    passesID         = ( id && no_exp_inner_trkr_hits && myTrigPresel );
    break;
  case electronID::electronPhys14L:
  case electronID::electronPhys14M:
  case electronID::electronPhys14T:
    id = PassElectronPhys14Id( iElectron, iElectronID );
    passesIso = id;
    passesID = id;
    passesKinematics = ((iElectron.pt() >= minElectronPt) && (fabs(iElectron.eta()) <= maxElectronEta) && !inCrack);
    break;
  }

  

  return (passesKinematics && passesIso && passesID);
}

bool
MiniAODHelper::isGoodTau(const pat::Tau& tau, const float min_pt, const tau::ID id)
{
  CheckVertexSetUp();

  double minTauPt = min_pt;

  bool passesKinematics = false;
  bool passesIsolation = false;
  bool passesID = tau.tauID("decayModeFinding") >= .5;

  switch (id) {
     case tau::nonIso:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronVLooseMVA5") >= .5;
        passesIsolation = true;
        break;
     case tau::loose:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronVLooseMVA5") >= .5;
        passesIsolation = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
     case tau::medium:
        passesID = passesID and \
                   tau.tauID("againstMuonLoose3") >= .5 and \
                   tau.tauID("againstElectronLooseMVA5") >= .5;
        passesIsolation = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
     case tau::tight:
        passesID = passesID and \
                   tau.tauID("againstMuonTight3") >= .5 and \
                   tau.tauID("againstElectronMediumMVA5") >= .5;
        passesIsolation = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") >= .5;
        break;
  }

  // systematics are only defined for p_T > 20
  passesKinematics = (tau.pt() >= 20) && (fabs(tau.eta()) <= 2.1) && (tau.pt() > minTauPt);

  return passesKinematics && passesIsolation && passesID;
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

  bool goodForMETCorrection = (
                iJet.correctedJet("Uncorrected").pt()>10.0 &&   
		(( !iJet.isPFJet() && iJet.emEnergyFraction()<0.9 ) || 
		( iJet.isPFJet() && (iJet.neutralEmEnergyFraction() + iJet.chargedEmEnergyFraction())<0.9 ))
		);

  if( fabs(iJet.eta())<2.4 ){
    loose = ( loose &&
	      iJet.chargedHadronEnergyFraction() > 0.0 &&
	      iJet.chargedMultiplicity() > 0
	      );
  }

  // Jet ID
  switch(iJetID){
  case jetID::jetMETcorrection:
    if( !goodForMETCorrection ) return false;
    break;
  case jetID::jetPU:
  case jetID::jetMinimal:
  case jetID::jetLooseAOD:
  case jetID::jetLoose:
  case jetID::jetTight:
    if( !loose ) return false;
    break;
  case jetID::none:
  default:
    break;
  }

  if( !PassesCSV(iJet, iCSVworkingPoint) ) return false;

  return true;
}


bool 
MiniAODHelper::isGoodTopJet(const boosted::HTTTopJet& iJet, const float iMinFatPt, const float iMaxAbsFatEta, const float iMinSubPt, const float iMaxAbsSubEta, const jetID::jetID iJetID){

  CheckVertexSetUp();
  
  // Fatjet requirements
  // Transverse momentum requirement
  if( iJet.fatjet.pt() < iMinFatPt ) return false;

  // Absolute eta requirement
  if( fabs(iJet.fatjet.eta()) > iMaxAbsFatEta ) return false;
  
  // Subjets requirements
  // Transverse momentum requirement
  if( iJet.nonW.pt() < iMinSubPt ) return false;
  if( iJet.W1.pt() < iMinSubPt ) return false;
  if( iJet.W2.pt() < iMinSubPt ) return false;

  // Absolute eta requirement
  if( fabs(iJet.nonW.eta()) > iMaxAbsSubEta ) return false;
  if( fabs(iJet.W1.eta()) > iMaxAbsSubEta ) return false;
  if( fabs(iJet.W2.eta()) > iMaxAbsSubEta ) return false;
  
  return true;
}


bool 
MiniAODHelper::isGoodHiggsJet(const boosted::SubFilterJet& iJet, const float iMinFatPt, const float iMaxAbsFatEta){

  CheckVertexSetUp();
  
  // Fatjet requirements
  // Transverse momentum requirement
  if( iJet.fatjet.pt() < iMinFatPt ) return false;

  // Absolute eta requirement
  if( fabs(iJet.fatjet.eta()) > iMaxAbsFatEta ) return false;
  
  // Filterjets requirements
  if( iJet.filterjets.size() < 2 ) return false;
  
  return true;
}


float MiniAODHelper::GetMuonRelIso(const pat::Muon& iMuon) const
{
  float result = 9999; 

  double pfIsoCharged = iMuon.pfIsolationR03().sumChargedHadronPt;
  double pfIsoNeutral = iMuon.pfIsolationR03().sumNeutralHadronEt + iMuon.pfIsolationR03().sumPhotonEt;

  double pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - 0.5*iMuon.pfIsolationR03().sumPUPt );

  result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
  
  return result;
}

//overloaded
float MiniAODHelper::GetMuonRelIso(const pat::Muon& iMuon,const coneSize::coneSize iconeSize, const corrType::corrType icorrType) const
{
  //rho corrections based on phys14
  //details here: https://www.dropbox.com/s/66lzhbro09diksa/effectiveareas-pog-121214.pdf?dl=0
  // !!! NOTE !!! rho used should be: fixedGridRhoFastjetAll
  float result = 9999; 
  
  double correction = 9999.;
  double EffArea = 9999.;
  double Eta = abs(iMuon.eta());
  
  double pfIsoCharged;
  double pfIsoNeutral;
  double pfIsoPUSubtracted;
  
  switch(iconeSize)
    {
    case coneSize::R04:
      pfIsoCharged = iMuon.pfIsolationR04().sumChargedHadronPt;
      pfIsoNeutral = iMuon.pfIsolationR04().sumNeutralHadronEt + iMuon.pfIsolationR04().sumPhotonEt;
      
      switch(icorrType)
	{
	case corrType::rhoEA:
	  //based on R04 Phys14_25ns_v1
	  if (Eta >= 0. && Eta < 0.8) EffArea = 0.1546;
	  else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.1325;
	  else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0913;
	  else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.1212;
	  else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.2085;
	  correction = useRho*EffArea;
	  break;
	case corrType::deltaBeta:
	  correction =  0.5*iMuon.pfIsolationR04().sumPUPt;
	  break;
	}
      
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
      break;

    case coneSize::R03:
      pfIsoCharged = iMuon.pfIsolationR03().sumChargedHadronPt;
      pfIsoNeutral = iMuon.pfIsolationR03().sumNeutralHadronEt + iMuon.pfIsolationR03().sumPhotonEt;
      
      switch(icorrType)
	{
	case corrType::rhoEA:
	  //effective area based on R03 Phys14_25ns_v1
	  if (Eta >= 0. && Eta < 0.8) EffArea = 0.0913;
	  else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.0765;
	  else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0546;
	  else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0728;
	  else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.1177;
	  correction = useRho*EffArea;
	  break;
	case corrType::deltaBeta:
	  correction = 0.5*iMuon.pfIsolationR03().sumPUPt;
	  break;
	}
      
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
      break;

    case coneSize::miniIso:
      double miniIsoR = 10.0/min(max(float(iMuon.pt()), float(50.)),float(200.));
      pfIsoCharged = isoSumRaw(charged_, iMuon, miniIsoR, 0.0001, 0.0, SelfVetoPolicy::selfVetoAll);
      pfIsoNeutral = isoSumRaw(neutral_, iMuon, miniIsoR, 0.01, 0.5, SelfVetoPolicy::selfVetoAll);
      switch(icorrType)
	{
	case corrType::rhoEA:
	  //effective area based on R03 Phys14_25ns_v1
	  if (abs(Eta) < 0.8) EffArea = 0.0913;
	  else if (abs(Eta) < 1.3) EffArea = 0.0765;
	  else if (abs(Eta) < 2.0) EffArea = 0.0546;
	  else if (abs(Eta) < 2.2) EffArea = 0.0728;
	  else EffArea = 0.1177;
	  correction = useRho*EffArea*(miniIsoR/0.3)*(miniIsoR/0.3);
	  break;
	case corrType::deltaBeta: 
	  double miniAbsIsoPU = isoSumRaw(pileup_, iMuon, miniIsoR, 0.01, 0.5, SelfVetoPolicy::selfVetoAll);
	  correction = 0.5*miniAbsIsoPU;
	  break;
	}
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction);
      result = (pfIsoCharged + pfIsoPUSubtracted)/iMuon.pt();
      break;
      
    }
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

//overloaded
float MiniAODHelper::GetElectronRelIso(const pat::Electron& iElectron,const coneSize::coneSize iconeSize, const corrType::corrType icorrType) const
{
  //rho*EA corrections based on phys14
  //details here: https://www.dropbox.com/s/66lzhbro09diksa/effectiveareas-pog-121214.pdf?dl=0
  // !!! NOTE !!! rho used should be: fixedGridRhoFastjetAll
  float result = 9999; 
  
  double correction = 9999.;
  double EffArea = 9999.;
  double Eta = abs(iElectron.eta());
  
  double pfIsoCharged;
  double pfIsoNeutral;
  double pfIsoPUSubtracted;
  
  switch(iconeSize)
    {
    case coneSize::R04:
    case coneSize::R03:
      pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
      pfIsoNeutral = iElectron.pfIsolationVariables().sumNeutralHadronEt + iElectron.pfIsolationVariables().sumPhotonEt;
      
      switch(icorrType)
	{
	case corrType::rhoEA:
	  if (Eta >= 0. && Eta < 0.8) EffArea = 0.1013;
      	  else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.0988;
      	  else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0572;
          else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0842;
          else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.1530;
	  if(!rhoIsSet) std::cout << " !! ERROR !! Trying to get rhoEffArea correction without setting rho" << std::endl;
	  correction = useRho*EffArea;
	  break;
	case corrType::deltaBeta:
	  correction = 0.5*iElectron.pfIsolationVariables().sumPUPt;
	  break;
	}
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction );
      result = (pfIsoCharged + pfIsoPUSubtracted)/iElectron.pt();
      break;
    case coneSize::miniIso:
      double innerR_ch;
      double innerR_nu;
      double miniIsoR = 10.0/min(max(float(iElectron.pt()), float(50.)),float(200.));
      if (iElectron.isEB())
	{ 
	  innerR_ch = 0.0;
	  innerR_nu = 0.0;
	} 
      else
	{
	  innerR_ch = 0.015;
	  innerR_nu = 0.08;
	}
      
      pfIsoCharged = isoSumRaw(charged_, iElectron, miniIsoR, innerR_ch, 0.0, SelfVetoPolicy::selfVetoNone);
      pfIsoNeutral = isoSumRaw(neutral_, iElectron, miniIsoR, innerR_nu, 0.0, SelfVetoPolicy::selfVetoNone, 22)+isoSumRaw(neutral_, iElectron, miniIsoR, 0.0, 0.0, SelfVetoPolicy::selfVetoNone, 130);
      switch(icorrType)
	{
	case corrType::rhoEA:
	  //effective area based on R03
	  if (Eta >= 0. && Eta < 0.8) EffArea = 0.1013;
      	  else if (Eta >= 0.8 && Eta < 1.3) EffArea = 0.0988;
          else if (Eta >= 1.3 && Eta < 2.0) EffArea = 0.0572;
          else if (Eta >= 2.0 && Eta < 2.2) EffArea = 0.0842;
          else if (Eta >= 2.2 && Eta <= 2.5) EffArea = 0.1530;
	  if(!rhoIsSet) std::cout << " !! ERROR !! Trying to get rhoEffArea correction without setting rho" << std::endl;
	  correction = useRho*EffArea*(miniIsoR/0.3)*(miniIsoR/0.3);
	  break;
	case corrType::deltaBeta: 
	  double miniAbsIsoPU = isoSumRaw(pileup_, iElectron, miniIsoR, innerR_ch, 0.0, SelfVetoPolicy::selfVetoNone);
	  correction = 0.5*miniAbsIsoPU;
	  break;
	}
      pfIsoPUSubtracted = std::max( 0.0, pfIsoNeutral - correction);
      result = (pfIsoCharged + pfIsoPUSubtracted)/iElectron.pt();
      break;
    }
  return result;
}


float MiniAODHelper::GetJetCSV(const pat::Jet& jet, const std::string taggername){
  
  float defaultFailure = -.1;
  
  float bTagVal = jet.bDiscriminator(taggername);

  if(isnan(bTagVal)) return defaultFailure;
  
  if(bTagVal > 1.) return 1.;
  if(bTagVal < 0.) return defaultFailure;
  
  return bTagVal;
}



bool MiniAODHelper::PassesCSV(const pat::Jet& iJet, const char iCSVworkingPoint){
  CheckSetUp();

  float csvValue = GetJetCSV(iJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags");

  // CSV b-tagging requirement
  switch(iCSVworkingPoint){
  case 'L':	if(csvValue > CSVLwp){ return true; }	break;
  case 'M':	if(csvValue > CSVMwp){ return true; }	break;
  case 'T':	if(csvValue > CSVTwp){ return true; }	break;
  case '-':	return true;                            break;
  }
  return false;
}


bool MiniAODHelper::PassElectronPhys14Id(const pat::Electron& iElectron, const electronID::electronID iElectronID) const{

  double SCeta = (iElectron.superCluster().isAvailable()) ? iElectron.superCluster()->position().eta() : -99;
  double absSCeta = fabs(SCeta);

  bool isEB = ( absSCeta < 1.479 );

  // double pfIsoCharged = iElectron.pfIsolationVariables().sumChargedHadronPt;
  // double pfIsoNeutralHadron = iElectron.pfIsolationVariables().sumNeutralHadronEt;
  // double pfIsoNeutralPhoton = iElectron.pfIsolationVariables().sumPhotonEt;
  // double pfIsoSumPUPt = iElectron.pfIsolationVariables().sumPUPt;

  // double relIso = (pfIsoCharged + std::max( pfIsoNeutralHadron + pfIsoNeutralPhoton - 0.5*pfIsoSumPUPt, 0.0 ))/iElectron.pt();
  double relIso = GetElectronRelIso(iElectron, coneSize::R03, corrType::rhoEA);

  double full5x5_sigmaIetaIeta = iElectron.full5x5_sigmaIetaIeta();
  double dEtaIn = fabs( iElectron.deltaEtaSuperClusterTrackAtVtx() );
  double dPhiIn = fabs( iElectron.deltaPhiSuperClusterTrackAtVtx() );
  double hOverE = iElectron.hcalOverEcal();

  double ooEmooP = -999;
  if( iElectron.ecalEnergy() == 0 ) ooEmooP = 1e30;
  else if( !std::isfinite(iElectron.ecalEnergy()) ) ooEmooP = 1e30;
  else ooEmooP = fabs(1.0/iElectron.ecalEnergy() - iElectron.eSuperClusterOverP()/iElectron.ecalEnergy() );

  double d0 = -999;
  double dZ = -999;
  double expectedMissingInnerHits = -999;
  if( iElectron.gsfTrack().isAvailable() ){
    d0 = fabs(iElectron.gsfTrack()->dxy(vertex.position()));
    dZ = fabs(iElectron.gsfTrack()->dz(vertex.position()));
    expectedMissingInnerHits = iElectron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  }

  bool passConversionVeto = ( iElectron.passConversionVeto() );

  bool pass = false;
  switch(iElectronID){
  case electronID::electronPhys14L:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.010331 &&
	       dEtaIn < 0.009277 &&
	       dPhiIn < 0.094739 &&
	       hOverE < 0.093068 &&
	       ooEmooP < 0.189968 &&
	       d0 < 0.035904 &&
	       dZ < 0.075496 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.130136
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.031838 &&
	       dEtaIn < 0.009833 &&
	       dPhiIn < 0.149934 &&
	       hOverE < 0.115754 &&
	       ooEmooP < 0.140662 &&
	       d0 < 0.099266 &&
	       dZ < 0.197897 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.163368
	       );
    }
    break;
  case electronID::electronPhys14M:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.009996 &&
	       dEtaIn < 0.008925 &&
	       dPhiIn < 0.035973 &&
	       hOverE < 0.050537  &&
	       ooEmooP < 0.091942 &&
	       d0 < 0.012235 &&
	       dZ < 0.042020 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.107587
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.030135 &&
	       dEtaIn < 0.007429 &&
	       dPhiIn < 0.067879 &&
	       hOverE < 0.086782 &&
	       ooEmooP < 0.100683 &&
	       d0 < 0.036719 &&
	       dZ < 0.138142 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.113254
	       );
    }
    break;
  case electronID::electronPhys14T:
    if( isEB ){
      pass = ( full5x5_sigmaIetaIeta < 0.009947 &&
	       dEtaIn < 0.006046 &&
	       dPhiIn < 0.028092 &&
	       hOverE < 0.045772 &&
	       ooEmooP < 0.020118 &&
	       d0 < 0.008790 &&
	       dZ < 0.021226 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.069537
	       );
    }
    else{
      pass = ( full5x5_sigmaIetaIeta < 0.028237 &&
	       dEtaIn < 0.007057 &&
	       dPhiIn < 0.030159 &&
	       hOverE < 0.067778 &&
	       ooEmooP < 0.098919 &&
	       d0 < 0.027984 &&
	       dZ < 0.133431 &&
	       expectedMissingInnerHits <= 1 &&
	       passConversionVeto &&
	       relIso < 0.078265
	       );
    }
    break;
  default:
    break;
  }

  return pass;
}

void MiniAODHelper::addVetos(const reco::Candidate &cand) {
  for (unsigned int i = 0, n = cand.numberOfSourceCandidatePtrs(); i < n; ++i) {
    const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
    if (cp.isNonnull() && cp.isAvailable()) vetos_.push_back(&*cp);
  }
}

void MiniAODHelper::clearVetos() {
  vetos_.clear();
}

float MiniAODHelper::isoSumRaw(const std::vector<const pat::PackedCandidate *> & cands, const reco::Candidate &cand, float dR, float innerR, float threshold, SelfVetoPolicy::SelfVetoPolicy selfVeto, int pdgId) const
{
  float dR2 = dR*dR, innerR2 = innerR*innerR;
  
  std::vector<const reco::Candidate *> vetos(vetos_);
  for (unsigned int i = 0, n = cand.numberOfSourceCandidatePtrs(); i < n; ++i) {
    if (selfVeto == SelfVetoPolicy::selfVetoNone) break;
    const reco::CandidatePtr &cp = cand.sourceCandidatePtr(i);
    if (cp.isNonnull() && cp.isAvailable()) {
      vetos.push_back(&*cp);
      if (selfVeto == SelfVetoPolicy::selfVetoFirst) break;
    }
  }
  
  typedef std::vector<const pat::PackedCandidate *>::const_iterator IT;
  IT candsbegin = std::lower_bound(cands.begin(), cands.end(), cand.eta() - dR, ByEta());
  IT candsend = std::upper_bound(candsbegin, cands.end(), cand.eta() + dR, ByEta());
  
  double isosum = 0;
  for (IT icharged = candsbegin; icharged < candsend; ++icharged) {
    // pdgId
    if (pdgId > 0 && abs((*icharged)->pdgId()) != pdgId) continue;
    // threshold
    if (threshold > 0 && (*icharged)->pt() < threshold) continue;
    // cone
    float mydr2 = reco::deltaR2(**icharged, cand);
    if (mydr2 >= dR2 || mydr2 < innerR2) continue;
    // veto
    if (std::find(vetos.begin(), vetos.end(), *icharged) != vetos.end()) {
      continue;
    }
    // add to sum
    isosum += (*icharged)->pt();
  }
  return isosum;
}

//// tt+X categorization -----------------------------
//// tt+b:  additionalJetEventId = 51
//// tt+2b:  additionalJetEventId = 52
//// tt+bb: additionalJetEventId = 53, 54 or 55
//// tt+c:  additionalJetEventId = 41 or 42
//// tt+cc: additionalJetEventId = 43, 44 or 45
//// tt+lf: additionalJetEventId = 0
int MiniAODHelper::ttHFCategorization(const std::vector<reco::GenJet>& genJets, const std::vector<int>& genBHadIndex, const std::vector<int>& genBHadJetIndex, const std::vector<int>& genBHadFlavour, const std::vector<int>& genBHadFromTopWeakDecay, const std::vector<reco::GenParticle>& genBHadPlusMothers, const std::vector<std::vector<int> >& genBHadPlusMothersIndices, const std::vector<int>& genBHadLeptonHadronIndex, const std::vector<int>& genBHadLeptonViaTau, const std::vector<int>& genCHadFlavour, const std::vector<int>& genCHadJetIndex, const std::vector<int>& genCHadFromTopWeakDecay, const std::vector<int>& genCHadBHadronId, const double genJetPtMin_, const double genJetAbsEtaMax_) {

    // Map <jet index, number of specific hadrons in the jet>
    // B jets with b hadrons directly from top quark decay
    std::map<int, int> bJetFromTopIds_all;
    // B jets with b hadrons directly from top quark decay
    std::map<int, int> bJetFromTopIds;
    // B jets with b hadrons after top quark decay
    std::map<int, int> bJetAfterTopIds;
    // B jets with b hadrons before top quark decay chain
    std::map<int, int> bJetBeforeTopIds;
    // C jets with c hadrons before top quark decay chain
    std::map<int, int> cJetBeforeTopIds;
    // C jets with c hadrons after top quark decay
    std::map<int, int> cJetAfterTopIds;
    
    // Counting number of specific hadrons in each b jet
    for(size_t hadronId = 0; hadronId < genBHadIndex.size(); ++hadronId) {
        // Flavour of the hadron's origin
        const int flavour = genBHadFlavour.at(hadronId);
	if(std::abs(flavour)==24)continue;
        // Whether hadron radiated before top quark decay
        // const bool fromTopDecay = genBHadFromTopWeakDecay.at(hadronId);
        // Index of a jet associated to the hadron
        const int jetIndex = genBHadJetIndex.at(hadronId);
        // Skipping hadrons which have no associated jet
        if(jetIndex < 0) continue;
        // Jet from direct top quark decay [pdgId(top)=6]
        if(std::abs(flavour) == 6) {
            if(bJetFromTopIds_all.count(jetIndex) < 1) bJetFromTopIds_all[jetIndex] = 1;
            else bJetFromTopIds_all[jetIndex]++;
        }
        // Skipping if jet is not in acceptance
        if(genJets.at(jetIndex).pt() < genJetPtMin_) continue;
        if(std::fabs(genJets.at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
        // Identifying jets with b hadrons not from top quark decay
        // Jet from direct top quark decay [pdgId(top)=6]
        if(std::abs(flavour) == 6) {
            if(bJetFromTopIds.count(jetIndex) < 1) bJetFromTopIds[jetIndex] = 1;
            else bJetFromTopIds[jetIndex]++;
        }
        // Skipping if jet is from top quark decay
        if(std::abs(flavour) == 6) continue;
        // Jet before top quark decay
        // if(!fromTopDecay) {
            if(bJetBeforeTopIds.count(jetIndex) < 1) bJetBeforeTopIds[jetIndex] = 1;
            else bJetBeforeTopIds[jetIndex]++;
        // }
        // // Jet after top quark decay but not directly from top
        // else if(fromTopDecay) {
        //     if(bJetAfterTopIds.count(jetIndex) < 1) bJetAfterTopIds[jetIndex] = 1;
        //     else bJetAfterTopIds[jetIndex]++;
        // }
    }
    
    // Counting number of specific hadrons in each c jet
    for(size_t hadronId = 0; hadronId < genCHadJetIndex.size(); ++hadronId) {
        // Skipping c hadrons that are coming from b hadrons
        if(genCHadBHadronId.at(hadronId) >= 0) continue;
        // Skipping c hadrons coming for W-dcays
        if(abs(genCHadFlavour.at(hadronId))==24) continue;
        // Index of a jet associated to the hadron
        const int jetIndex = genCHadJetIndex.at(hadronId);
        // Whether hadron radiated before top quark decay
        // const bool fromTopDecay = genCHadFromTopWeakDecay.at(hadronId);
        // Skipping hadrons which have no associated jet
        if(jetIndex < 0) continue;
        // Skipping if jet is not in acceptance
        if(genJets.at(jetIndex).pt() < genJetPtMin_) continue;
        if(std::fabs(genJets.at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
        // Jet before top quark decay
        // if(!fromTopDecay) {
            if(cJetBeforeTopIds.count(jetIndex) < 1) cJetBeforeTopIds[jetIndex] = 1;
            else cJetBeforeTopIds[jetIndex]++;
        // }
        // // Jet after top quark decay but not directly from top
        // else if(fromTopDecay) {
        //     if(cJetAfterTopIds.count(jetIndex) < 1) cJetAfterTopIds[jetIndex] = 1;
        //     else cJetAfterTopIds[jetIndex]++;
        // }
    }
    
    // Finding additional b jets (before top decay)
    std::vector<int> additionalBJetIds;
    for(std::map<int, int>::iterator it = bJetBeforeTopIds.begin(); it != bJetBeforeTopIds.end(); ++it) {
        const int jetId = it->first;
        // Skipping the jet if it contains a b hadron directly from top quark decay
        if(bJetFromTopIds.count(jetId) > 0) continue;
        additionalBJetIds.push_back(jetId);
    }
    // Finding pseudo-additional b jets (after top decay)
    std::vector<int> pseudoadditionalBJetIds;
    for(std::map<int, int>::iterator it = bJetAfterTopIds.begin(); it != bJetAfterTopIds.end(); ++it) {
        const int jetId = it->first;
        // Skipping the jet if it contains a b hadron directly from top quark decay
        if(bJetFromTopIds.count(jetId) > 0) continue;
        pseudoadditionalBJetIds.push_back(jetId);
    }
    // Finding additional c jets
    std::vector<int> additionalCJetIds;
    for(std::map<int, int>::iterator it = cJetBeforeTopIds.begin(); it != cJetBeforeTopIds.end(); ++it) {
        const int jetId = it->first;
        if(bJetFromTopIds.count(jetId) > 0) continue;
        additionalCJetIds.push_back(jetId);
    }
    // Finding pseudo-additional c jets (after top decay)
    std::vector<int> pseudoadditionalCJetIds;
    for(std::map<int, int>::iterator it = cJetAfterTopIds.begin(); it != cJetAfterTopIds.end(); ++it) {
        const int jetId = it->first;
        // Skipping the jet if it contains a b hadron directly from top quark decay
        if(bJetFromTopIds.count(jetId) > 0) continue;
        pseudoadditionalCJetIds.push_back(jetId);
    }
    
    // Categorizing event based on number of additional b/c jets 
    // and number of corresponding hadrons in each of them
    // int additionalJetEventId;
    int additionalJetEventId = bJetFromTopIds.size()*100;
    // tt + 1 additional b jet
    if (additionalBJetIds.size() == 1) {
        int nHadronsInJet = bJetBeforeTopIds[additionalBJetIds.at(0)];
        // tt + 1 additional b jet from 1 additional b hadron
        if(nHadronsInJet == 1) additionalJetEventId = 51;
        // tt + 1 additional b jet from >=2 additional b hadrons
        else additionalJetEventId = 52;
    }
    // tt + 2 additional b jets
    else if (additionalBJetIds.size() > 1) {
        int nHadronsInJet1 = bJetBeforeTopIds[additionalBJetIds.at(0)];
        int nHadronsInJet2 = bJetBeforeTopIds[additionalBJetIds.at(1)];
        // tt + 2 additional b jets each from 1 additional b hadron
        if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId = 53;
        // tt + 2 additional b jets one of which from >=2 overlapping additional b hadrons
        else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 54;
        // tt + 2 additional b jets each from >=2 additional b hadrons
        else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 55;
    }
    // tt + no additional b jets
    else if(additionalBJetIds.size() == 0) {
        // tt + >=1 pseudo-additional b jet with b hadrons after top quark decay
        if(pseudoadditionalBJetIds.size() > 0) additionalJetEventId = 56;
        // tt + 1 additional c jet
        else if(additionalCJetIds.size() == 1) {
            int nHadronsInJet = cJetBeforeTopIds[additionalCJetIds.at(0)];
            // tt + 1 additional c jet from 1 additional c hadron
            if(nHadronsInJet == 1) additionalJetEventId = 41;
            // tt + 1 additional c jet from >=2 overlapping additional c hadrons
            else additionalJetEventId = 42;
        }
        // tt + >=2 additional c jets
        else if(additionalCJetIds.size() > 1) {
            int nHadronsInJet1 = cJetBeforeTopIds[additionalCJetIds.at(0)];
            int nHadronsInJet2 = cJetBeforeTopIds[additionalCJetIds.at(1)];
            // tt + 2 additional c jets each from 1 additional c hadron
            if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId = 43;
            // tt + 2 additional c jets one of which from >=2 overlapping additional c hadrons
            else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 44;
            // tt + 2 additional c jets each from >=2 additional c hadrons
            else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId = 45;
        }
        // tt + no additional c jets
        else if(additionalCJetIds.size() == 0) {
            // tt + >=1 pseudo-additional c jet with c hadrons after top quark decay
            if(pseudoadditionalCJetIds.size() > 0) additionalJetEventId = 46;
            // tt + light jets
            else additionalJetEventId = 0;
        }
    }

    return additionalJetEventId;
}



    ///////////////////
    /// Higgs Decay ///
    ///////////////////

int MiniAODHelper::GetHiggsDecay(edm::Handle<std::vector<reco::GenParticle> >& mcparticles){

  int Hdecay = -1;
    
  if( mcparticles.isValid() ){

    Hdecay=0;
  
    for( size_t k = 0; k < mcparticles->size(); k++ ){
      const reco::Candidate & mcParticle = (*mcparticles)[k];

      int status = mcParticle.status();
      int pdgId  = mcParticle.pdgId();
      int absId  = abs( pdgId );
      int numdgt = mcParticle.numberOfDaughters();

      //// must be a Higgs and status 62(pythia 8) 
      if( absId!=25 || status!=62  ) continue;
      if (!(numdgt>1))continue;
	  
      int d0=-99, d1=-99;
      
      int ind0 = 0;
      int ind1 = 1;

      if( numdgt>2 ){
	if( mcParticle.daughter(0)->pdgId()==pdgId ){
	  ind0 = 1;
	  ind1 = 2;
	}
	if( mcParticle.daughter(1)->pdgId()==pdgId ){
	  ind0 = 0;
	  ind1 = 2;
	}
      }

      d0 = mcParticle.daughter(ind0)->pdgId();
      d1 = mcParticle.daughter(ind1)->pdgId();

      d0 = abs(d0);
      d1 = abs(d1);

      if( d0==5 && d1==5 ) Hdecay = 1;  //bb
      if( d0==24 && d1==24) Hdecay = 2;   //WW
      if( d0==15 && d1==15) Hdecay = 3;  //TauTau
      if( d0==21 && d1==21) Hdecay = 4;  //glueglue
      if( d0==4 && d1==4 ) Hdecay = 5;  //cc
      if( d0==23 && d1==23) Hdecay = 6;  //ZZ

      if( d0==22 && d1==23) Hdecay = 7;  //Zy
      if( d0==23 && d1==22) Hdecay = 7;  //Zy

      if( d0==22 && d1==22) Hdecay = 8;  //yy

      if( d0==21 && d1==22) Hdecay = 9; //gy
      if( d0==22 && d1==21) Hdecay = 9; //gy
      if( d0==3 && d1==3) Hdecay = 10; //ss
	
      if( d0==13 && d1==13) Hdecay = 11; //mumu
      
      if( (Hdecay==0) && (d0==22 || d1==22)) Hdecay =12; //?y
      if( (Hdecay==0) && (d0==21 || d1==21)) Hdecay = 13; //?g
      if( (Hdecay==0) && (d0>100 || d1>100)) Hdecay = 14; //?Hadron
      
      if( d0==1 && d1==1) Hdecay = 15; //uu
      if( d0==2 && d1==2) Hdecay = 16; //dd
      if( d0==6 && d1==6) Hdecay = 17; //tt
      if( d0==11 && d1==11) Hdecay = 18; //ee

    }
    
  }
  

  return Hdecay;
}


std::vector<pat::Jet> MiniAODHelper::GetDeltaRCleanedJets(
    const std::vector<pat::Jet> &inputJets, const std::vector<pat::Muon>& inputMuons, const std::vector<pat::Electron>& inputElectrons, const double deltaRCut)
{
	CheckSetUp();

	
	std::vector<pat::Jet> outputJets;

	for( std::vector<pat::Jet>::const_iterator iJet = inputJets.begin(); iJet!=inputJets.end(); ++iJet ){

	  bool isOverlap = false;
	  
	  TLorentzVector jet_p4;
	  jet_p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
	  
	  for( std::vector<pat::Electron>::const_iterator iEle = inputElectrons.begin(); iEle != inputElectrons.end(); iEle++ ){ 
	    TLorentzVector ele_p4;
	    ele_p4.SetPxPyPzE(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());
	    double delta_tmp = jet_p4.DeltaR(ele_p4);
	    if(delta_tmp < deltaRCut){
	      isOverlap = true;
	      break;
	    } 
	  }
	  
	  if( isOverlap ) continue;
	  
	  for( std::vector<pat::Muon>::const_iterator iMuon = inputMuons.begin(); iMuon != inputMuons.end(); iMuon++ ){ 
	    TLorentzVector muon_p4;
	    muon_p4.SetPxPyPzE(iMuon->px(),iMuon->py(),iMuon->pz(),iMuon->energy());
	    double delta_tmp = jet_p4.DeltaR(muon_p4);
	    if(delta_tmp < deltaRCut){
	      isOverlap = true;
	      break;
	    } 
	  }
	  
	  if( isOverlap ) continue;
	  
	  outputJets.push_back(*iJet);
	  
	}
	
	
	return outputJets;
}


/// JER function
double MiniAODHelper::getJERfactor( const int returnType, const double jetAbsETA, const double genjetPT, const double recojetPT){

  // CheckSetUp();
  // string samplename = GetSampleName();
  double factor = 1.;
    
  double scale_JER = 1., scale_JERup = 1., scale_JERdown = 1.;

  //// nominal SFs have changed since run1, and the new up/down SFs are still unknown???
  if( jetAbsETA<0.5 ){ 
    scale_JER = 1.079; scale_JERup = 1.105; scale_JERdown = 1.053;
  }
  else if( jetAbsETA<1.1 ){ 
    scale_JER = 1.099; scale_JERup = 1.127; scale_JERdown = 1.071;
  }
  else if( jetAbsETA<1.7 ){ 
    scale_JER = 1.121; scale_JERup = 1.150; scale_JERdown = 1.092;
  }
  else if( jetAbsETA<2.3 ){ 
    scale_JER = 1.208; scale_JERup = 1.254; scale_JERdown = 1.162;
  }
  else if( jetAbsETA<2.8 ){ 
    scale_JER = 1.254; scale_JERup = 1.316; scale_JERdown = 1.192;
  }
  else if( jetAbsETA<3.2 ){ 
    scale_JER = 1.395; scale_JERup = 1.458; scale_JERdown = 1.332;
  }
  else if( jetAbsETA<5.0 ){ 
    scale_JER = 1.056; scale_JERup = 1.247; scale_JERdown = 0.865;
  }

  double jetPt_JER = recojetPT;
  double jetPt_JERup = recojetPT;
  double jetPt_JERdown = recojetPT;

  double diff_recojet_genjet = recojetPT - genjetPT;

  if( genjetPT>10. ){
    jetPt_JER = std::max( 0., genjetPT + scale_JER * ( diff_recojet_genjet ) );
    jetPt_JERup = std::max( 0., genjetPT + scale_JERup * ( diff_recojet_genjet ) );
    jetPt_JERdown = std::max( 0., genjetPT + scale_JERdown * ( diff_recojet_genjet ) );
  }

  if( returnType==1 )       factor = jetPt_JERup/recojetPT;
  else if( returnType==-1 ) factor = jetPt_JERdown/recojetPT;
  else                      factor = jetPt_JER/recojetPT;

  if( !(genjetPT>10.) ) factor = 1.;

  return factor;
}

std::vector<pat::MET> MiniAODHelper::CorrectMET(const std::vector<pat::Jet>& oldJetsForMET, const std::vector<pat::Jet>& newJetsForMET, const std::vector<pat::MET>& pfMETs){
  // this function takes two jet collections and replaces their contribution to the Type1 correction of the MET

  std::vector<pat::MET> outputMets;

  for(std::vector<pat::MET>::const_iterator oldMET=pfMETs.begin();oldMET!=pfMETs.end();++oldMET){
    pat::MET outMET=*oldMET; 

    if(oldMET-pfMETs.begin() == 0){
    //get old MET p4
    TLorentzVector oldMETVec;
    oldMETVec.SetPxPyPzE(oldMET->p4().Px(),oldMET->p4().Py(),oldMET->p4().Pz(),oldMET->p4().E());
    // add the pT vector of the old jets with the initial correction to the MET vector
    for(std::vector<pat::Jet>::const_iterator itJet=oldJetsForMET.begin();itJet!=oldJetsForMET.end();++itJet){
      TLorentzVector oldJETVec;
      oldJETVec.SetPtEtaPhiE(itJet->pt(),itJet->eta(),itJet->phi(),itJet->energy());
      TLorentzVector PToldJETVec;
      PToldJETVec.SetPxPyPzE(oldJETVec.Px(),oldJETVec.Py(),0.0,oldJETVec.Et());
      oldMETVec+=PToldJETVec;
    }
    // now subtract the pT vectors of the clean recorrected jets
    for(std::vector<pat::Jet>::const_iterator itJet=newJetsForMET.begin();itJet!=newJetsForMET.end();++itJet){
      TLorentzVector newJETVec;
      newJETVec.SetPtEtaPhiE(itJet->pt(),itJet->eta(),itJet->phi(),itJet->energy());
      TLorentzVector PTnewJETVec;
      PTnewJETVec.SetPxPyPzE(newJETVec.Px(),newJETVec.Py(),0.0,newJETVec.Et());
      oldMETVec-=PTnewJETVec;
    }
    outMET.setP4(reco::Candidate::LorentzVector(oldMETVec.Px(),oldMETVec.Py(),oldMETVec.Pz(),oldMETVec.E()));
    }
    
    outputMets.push_back(outMET);
  }

  return outputMets;

}

