#ifndef MINIAOD_BOOSTEDOBJECTS_BOOSTEDJET_H
#define MINIAOD_BOOSTEDOBJECTS_BOOSTEDJET_H

#include <vector>

#include "DataFormats/PatCandidates/interface/Jet.h"

namespace boosted {

  enum JetType{ Top, Higgs, NA };
  
  class BoostedJet {
    
    public:
    
      BoostedJet() :
		    fatjet(pat::Jet()),
        topjet(pat::Jet()),
		    nonW(pat::Jet()),
		    W1(pat::Jet()),
		    W2(pat::Jet()),
		    fatjetMass(-99),
        fatjetPt(-99),
        fatjetEta(-99),
        fatjetPhi(-99),
        topMass(-99),
        unfilteredMass(-99),
        prunedMass(-99),
        fRec(-99),
        massRatioPassed(-99),
        tau1Unfiltered(-99),
	      tau2Unfiltered(-99),
	      tau3Unfiltered(-99),
	      tau1Filtered(-99),
	      tau2Filtered(-99),
	      tau3Filtered(-99),
	      tau1Softdrop(-99),
	      tau2Softdrop(-99),
	      tau3Softdrop(-99),
        qWeight(-99),
        qEpsilon(-99),
        qSigmaM(-99) {};
      
      math::XYZTLorentzVector GetWJetVec() const{
        
        if(W1.pt()<=0) return math::XYZTLorentzVector();
        
        math::XYZTLorentzVector wjet = W1.p4();
        wjet += W2.p4();
        
        return wjet;
      }
      
      math::XYZTLorentzVector GetTopJetVec() const{
        return topjet.p4();
      }
      
      // Fatjet
      pat::Jet fatjet;
      
      // HTT V2 Information
      pat::Jet topjet;
	    pat::Jet nonW;
	    pat::Jet W1;
	    pat::Jet W2;

	    double fatjetMass;
      double fatjetPt;
      double fatjetEta;
      double fatjetPhi;
      
      double topMass;
      double unfilteredMass;
      double prunedMass;
      double fRec;
      double massRatioPassed;
      
      double Ropt;
      double RoptCalc;
      double ptForRoptCalc;
      
      double tau1Unfiltered;
	    double tau2Unfiltered;
	    double tau3Unfiltered;
	    double tau1Filtered;
	    double tau2Filtered;
	    double tau3Filtered;
	    float tau1Softdrop;
	    float tau2Softdrop;
	    float tau3Softdrop;
      
      double qWeight;
      double qEpsilon;
      double qSigmaM;
      
      // Subjet Filterjet Information
      std::vector<pat::Jet> subjets;
      std::vector<pat::Jet> filterjets;
      
      // Pruned Jet Information
      std::vector<pat::Jet> prunedsubjets;
      
      // Soft Drop Jet Information
      std::vector<pat::Jet> sdsubjets;
      
      // Soft Drop Z2 B1 Jet Information
      std::vector<pat::Jet> sdz2b1subjets;
  };

  typedef std::vector<BoostedJet> BoostedJetCollection;
  
}

#endif
