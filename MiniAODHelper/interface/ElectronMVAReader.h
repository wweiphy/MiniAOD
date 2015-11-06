#ifndef MINIAOD_MINIAODHELPER_ELECTRONMVAREADER_H
#define MINIAOD_MINIAODHELPER_ELECTRONMVAREADER_H

#include <iostream>
#include "TH1F.h"
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

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "CommonTools/Utils/interface/normalizedPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

class ElectronMVAReader{
  
  public:
    
    // Constructor & Destructor
    ElectronMVAReader(std::string WeightPath_ = "");
    ~ElectronMVAReader();
    
    // Return the Output of the reader
    float GetElectronMVAReaderOutput(const pat::Electron& iElectron, const edm::Handle< reco::ConversionCollection >& conversions, const edm::Handle< reco::BeamSpot >& theBeamSpot, bool verbose = false);
    
  private:

    std::string weightsPath;

    std::vector<std::string> TMVAVarNames;
    float* TMVAVars;
    TMVA::Reader* TMVAReader;
    
    // Get the TMVA input variable names from the weight file
    void GetTMVAVarNames(std::string filePath_, bool verbose = false);

    // Initialize the TMVA input variables given in the weight file
    void GetTMVAVars(std::string filePath_, bool verbose = false);
    
    // Clear the values of the TMVA input variables
    void ResetTMVAVars();
};

#endif
