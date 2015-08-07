#ifndef MINIAOD_MINIAODHELPER_TOPTAGGER_H
#define MINIAOD_MINIAODHELPER_TOPTAGGER_H

#include <map>

#include "TH1F.h"
#include "TFile.h"
#include "TMVA/Reader.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"
#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

namespace TopTag{
  enum Mode{ HEP, Likelihood, TMVA };
  enum SubjetAssign{ Pt, Wmass, CSV };
}

class TopTagger{

  public:

    // Constructor & Destructor
    TopTagger(TopTag::Mode mode_=TopTag::Likelihood, TopTag::SubjetAssign subjetAssign_=TopTag::CSV, std::string filePath_ = "toplikelihoodtaggerhistos.root");
    ~TopTagger();
    
    // Return the Output of the Top Tagger
    float GetTopTaggerOutput(const boosted::HTTTopJet& topjet, bool verbose = false);

    // Comparison function for Top Tagger output
    bool FirstHasHigherTopTaggerOutput(boosted::HTTTopJet jet1, boosted::HTTTopJet jet2);

    // Sorting function for Top Tagger output
    boosted::HTTTopJetCollection GetSortedByTopTaggerOutput(const boosted::HTTTopJetCollection& topJets, bool verbose = false);

    // Find a hadronic top candidate in Top jet collection
    float GetTopHad(boosted::HTTTopJetCollection& topJets, boosted::HTTTopJet& topHadCand, bool verbose = false);

  private:

    // MiniAODHelper 
    MiniAODHelper* helper;
    
    // Top Tagger mode
    TopTag::Mode mode;

    // Mode for subjet assignment
    TopTag::SubjetAssign subjetAssign;
    
    // b-tagger name
    const char* btagger;
    
    // Likelihood Top Tagger
    // Input file containing Likelihood histograms
    TFile* file;

    // Histograms for Likelihood Top Tagger
    TH1F* mtop_top_histo;
    TH1F* mtop_nottop_histo;
    TH1F* mratio_top_histo;
    TH1F* mratio_nottop_histo;
    TH1F* atan_top_histo;
    TH1F* atan_nottop_histo;

    // TMVA Top Tagger
    // Variables and TMVA Reader for TMVA Top Tagger
    std::vector<std::string> TMVAVarNames;
    float* TMVAVars;
    TMVA::Reader* TMVAReader;

    // CSV definition for subjets: subjets[0] = b, subjets[1] = W1, subjets[2] = W2
    void TopSubjetCSVDef(std::vector<pat::Jet> &subjets);

    // Get the TMVA input variable names from the weight file
    void GetTMVAVarNames(std::string filePath_, bool verbose = false);

    // Initialize the TMVA input variables given in the weight file
    void GetTMVAVars(std::string filePath_, bool verbose = false);

    // Clear the values of the TMVA input variables
    void ResetTMVAVars();  
};


struct TopTaggerOutputComparison{

  TopTagger* topTagger;

  TopTaggerOutputComparison(TopTagger* topTagger_): topTagger(topTagger_){};

  bool operator()(boosted::HTTTopJet jet1, boosted::HTTTopJet jet2){
    return topTagger->GetTopTaggerOutput(jet1)>topTagger->GetTopTaggerOutput(jet2);
  };
};
  
#endif
