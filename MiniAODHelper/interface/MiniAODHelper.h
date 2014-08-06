#ifndef _MiniAODHelper_h
#define _MiniAODHelper_h

#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <algorithm>
#include "TVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "CommonTools/Utils/interface/normalizedPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#endif

typedef std::map<std::string, std::string> mparams;
typedef std::vector< TLorentzVector > vecTLorentzVector;
typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;
typedef std::vector<int> vint;

namespace analysisType{ enum analysisType{ LJ, DIL, TauLJ, TauDIL }; }
namespace sysType{enum sysType{NA, JERup, JERdown, JESup, JESdown, hfSFup, hfSFdown, lfSFdown, lfSFup, TESup, TESdown, CSVLFup, CSVLFdown, CSVHFup, CSVHFdown, CSVHFStats1up, CSVHFStats1down, CSVLFStats1up, CSVLFStats1down, CSVHFStats2up, CSVHFStats2down, CSVLFStats2up, CSVLFStats2down, CSVCErr1up, CSVCErr1down, CSVCErr2up, CSVCErr2down }; }
namespace jetID{		enum jetID{			none, jetMinimal, jetLooseAOD, jetLoose, jetTight }; }
namespace tauID { enum tauID{ tauNonIso, tauLoose, tauMedium, tauTight }; }
namespace muonID{		enum muonID{		muonSide, muonSideLooseMVA, muonSideTightMVA, muonLoose, muonTight, muonPtOnly, muonPtEtaOnly, muonPtEtaIsoOnly, muonPtEtaIsoTrackerOnly, muonNoCuts }; }
namespace electronID{	enum electronID{	electronSide, electronSideLooseMVA, electronSideTightMVA, electronLoose, electronTight, electronTightMinusTrigPresel, electronLooseMinusTrigPresel, electronNoCuts }; }
namespace hdecayType{	enum hdecayType{ hbb, hcc, hww, hzz, htt, hgg, hjj, hzg }; }


using namespace std;

class MiniAODHelper{

  // === Functions === //
 public: 
  // Constructor(s) and destructor
  MiniAODHelper();
  virtual ~MiniAODHelper();

  // Set up MiniAODHelper
  void SetUp(string, int, const analysisType::analysisType, bool);

  void SetVertex(const reco::Vertex&);

  void SetRho(double);

  void SetJetCorrector(const JetCorrector*);

  void SetFactorizedJetCorrector();

  std::vector<pat::Muon> GetSelectedMuons(const std::vector<pat::Muon>&, const float, const muonID::muonID);
  std::vector<pat::Electron> GetSelectedElectrons(const std::vector<pat::Electron>&, const float, const electronID::electronID);
  std::vector<pat::Jet> GetSelectedJets(const std::vector<pat::Jet>&, const float, const float, const jetID::jetID, const char);
  std::vector<pat::Jet> GetUncorrectedJets(const std::vector<pat::Jet>&);
  std::vector<pat::Jet> GetCorrectedJets(const std::vector<pat::Jet>&, const edm::Event&, const edm::EventSetup&);
  std::vector<pat::Jet> GetCorrectedJets(const std::vector<pat::Jet>&);
  bool isGoodMuon(const pat::Muon&, const float, const muonID::muonID);
  bool isGoodElectron(const pat::Electron&, const float, const electronID::electronID);
  bool isGoodJet(const pat::Jet&, const float, const float, const jetID::jetID, const char);
  float GetMuonRelIso(const pat::Muon&) const;
  float GetElectronRelIso(const pat::Electron&) const;
  bool PassesCSV(const pat::Jet&, const char);

  template <typename T, typename S> std::vector<T> RemoveOverlaps( const std::vector<S>&, const std::vector<T>& );
  template <typename T, typename S> T RemoveOverlap( const std::vector<S>&, const T& );

  template <typename T, typename S> double DeltaR( const S&, const T& );


 private:
  bool isSetUp;
  bool vertexIsSet;
  bool rhoIsSet;
  bool jetcorrectorIsSet;
  bool factorizedjetcorrectorIsSet;
  string era;
  int sampleNumber;
  bool isData;
  analysisType::analysisType analysis;
  string samplename;

  float CSVLwp, CSVMwp, CSVTwp;

  double useRho;

  reco::Vertex vertex;

  const JetCorrector* corrector;
  FactorizedJetCorrector* useJetCorrector;

  inline void ThrowFatalError(const std::string& m) const { cerr << "[ERROR]\t" << m << " Cannot continue. Terminating..." << endl; exit(1); };

  inline void CheckSetUp() const { if(!isSetUp){ ThrowFatalError("MiniAODHelper not yet set up."); } };
  inline void CheckVertexSetUp() const { if(!vertexIsSet){ ThrowFatalError("Vertex is not set."); } };

}; // End of class prototype


template <typename PATObj1, typename PATObj2> 
PATObj1 MiniAODHelper::RemoveOverlap( const std::vector<PATObj2>& other, const PATObj1& unclean ){

  unsigned int nSources1 = unclean.numberOfSourceCandidatePtrs();
  bool hasOverlaps = false;

  std::vector<reco::CandidatePtr> overlaps;

  for( typename std::vector<PATObj2>::const_iterator iobj2 = other.begin(); iobj2!=other.end(); ++iobj2 ){

    unsigned int nSources2 = iobj2->numberOfSourceCandidatePtrs();

    for( unsigned int i1=0; i1<nSources1; i1++ ){
      bool uncleanSourceHasOverlap = false;

      reco::CandidatePtr source1 = unclean.sourceCandidatePtr(i1);

      if( !(source1.isNonnull() && source1.isAvailable()) ) continue;

      for( unsigned int i2=0; i2<nSources2; i2++ ){

	reco::CandidatePtr source2 = iobj2->sourceCandidatePtr(i2);

	if( !(source2.isNonnull() && source2.isAvailable()) ) continue;

	if( source1==source2 ){
	  hasOverlaps = true;
	  overlaps.push_back(source2);
	  uncleanSourceHasOverlap = true;
	}
      }

      if( !uncleanSourceHasOverlap ){
	if( (abs((*source1).pdgId())==11 && iobj2->isElectron()) ||
	    (abs((*source1).pdgId())==13 && iobj2->isMuon()) ){

	  double deltaR = reco::deltaR((*source1).eta(), (*source1).phi(), iobj2->eta(), iobj2->phi());
	  /* std::cout << "\t deltaR = " << deltaR << std::endl; */
	  /* std::cout << "\t source: eta = " << (*source1).eta() << ",\t phi = " << (*source1).phi() << ",\t pt = " << (*source1).pt() << std::endl; */
	  /* std::cout << "\t iobj2:  eta = " << iobj2->eta() << ",\t phi = " << iobj2->phi() << ",\t pt = " << iobj2->pt() << std::endl; */
	  if( deltaR<0.01 ){
	    overlaps.push_back(source1);
	    hasOverlaps = true;
	  }
	}
      }

    }
  }// end loop over iobj22


  PATObj1 cleaned = unclean;
  if( hasOverlaps ){
    math::XYZTLorentzVector original = cleaned.p4();

    for( int iOverlap=0; iOverlap<int(overlaps.size()); iOverlap++ ){

      const reco::Candidate & cOverlap = *(overlaps[iOverlap]);
      math::XYZTLorentzVector overlaper = cOverlap.p4();

      original -= overlaper;
    }

    cleaned.setP4( original );
  }

  return cleaned;
}



template <typename PATObj1, typename PATObj2> 
std::vector<PATObj1> MiniAODHelper::RemoveOverlaps( const std::vector<PATObj2>& other, const std::vector<PATObj1>& unclean ){

  std::vector<PATObj1> cleaned;
  
  for( typename std::vector<PATObj1>::const_iterator iobj1 = unclean.begin(); iobj1!=unclean.end(); ++iobj1 ){

    PATObj1 myobj = (*iobj1);
    PATObj1 clean = RemoveOverlap(other, myobj);

    cleaned.push_back(clean);
  }

  return cleaned;
}


template <typename PATObj1, typename PATObj2> 
double MiniAODHelper::DeltaR( const PATObj2& two, const PATObj1& one ){

  double deltaR = reco::deltaR( one->eta(), one->phi(), two->eta(), two->phi() );
  return deltaR;
}

#endif // _MiniAODHelper_h
