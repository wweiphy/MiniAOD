#ifndef MINIAOD_MINIAODHELPER_ElectronMVAReader_H
#define MINIAOD_MINIAODHELPER_ElectronMVAReader_H

#include <map>

#include "TH1F.h"
#include "TFile.h"
#include "TMVA/Reader.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "MiniAOD/BoostedObjects/interface/BoostedJet.h"
#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"

class ElectronMVAReader{
  
  public:
    
    // Constructor & Destructor
    ElectronMVAReader(std::string WeightPath_ = "");
    ~ElectronMVAReader();
    
    // Return the Output of the reader
    float GetElectronMVAReaderOutput(const pat::Electron& iElectron, bool verbose = false);
    
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
