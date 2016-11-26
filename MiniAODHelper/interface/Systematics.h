#ifndef SYSTEMATICS_H
#define SYSTEMATICS_H

#include <string>

namespace sysType {
  enum sysType{
    NA,

    // total JEC uncertainties
    JESup,			
    JESdown,			

    // individual JEC uncertainties up
    JESPileUpDataMCup,
    JESPileUpPtRefup,
    JESPileUpPtBBup,			
    JESPileUpPtEC1up,
    JESPileUpPtEC2up,
    JESPileUpPtHFup,
    JESRelativeJEREC1up,
    JESRelativeJEREC2up,		
    JESRelativeJERHFup,
    JESRelativeFSRup,
    JESRelativeStatFSRup,
    JESRelativeStatEC2up,
    JESRelativeStatECup,		
    JESRelativeStatHFup,
    JESRelativePtBBup,
    JESRelativePtEC1up,
    JESRelativePtEC2up,
    JESRelativePtHFup,		
    JESTimeEtaup,
    JESAbsoluteScaleup,
    JESAbsoluteMPFBiasup,
    JESAbsoluteStatup,
    JESSinglePionECALup,		
    JESSinglePionHCALup,
    JESFragmentationup,
    JESTimePtup,
    JESFlavorQCDup,			

    // individual JEC uncertainties down
    JESPileUpDataMCdown,		
    JESPileUpPtRefdown,
    JESPileUpPtBBdown,			
    JESPileUpPtEC1down,
    JESPileUpPtEC2down,
    JESPileUpPtHFdown,
    JESRelativeJEREC1down,
    JESRelativeJEREC2down,		
    JESRelativeJERHFdown,
    JESRelativeFSRdown,
    JESRelativeStatFSRdown,
    JESRelativeStatEC2down,
    JESRelativeStatECdown,		
    JESRelativeStatHFdown,
    JESRelativePtBBdown,
    JESRelativePtEC1down,
    JESRelativePtEC2down,
    JESRelativePtHFdown,		
    JESTimeEtadown,
    JESAbsoluteScaledown,
    JESAbsoluteMPFBiasdown,
    JESAbsoluteStatdown,
    JESSinglePionECALdown,		
    JESSinglePionHCALdown,
    JESFragmentationdown,
    JESTimePtdown,
    JESFlavorQCDdown,			

    // JER uncertainty
    JERup,			
    JERdown,

    hfSFup,
    hfSFdown,
    lfSFdown,
    lfSFup,

    TESup,
    TESdown,

    CSVLFup,
    CSVLFdown,
    CSVHFup,
    CSVHFdown,
    CSVHFStats1up,
    CSVHFStats1down,
    CSVLFStats1up,
    CSVLFStats1down,
    CSVHFStats2up,
    CSVHFStats2down,
    CSVLFStats2up,
    CSVLFStats2down,
    CSVCErr1up,
    CSVCErr1down,
    CSVCErr2up,
    CSVCErr2down 
  };

  // convert between string and int representation
  sysType get(const std::string& name);
  std::string toString(const sysType type);

  // true if type is one of the JEC-related uncertainties and up
  bool isJECUncertaintyUp(const sysType type);

  // true if type is one of the JEC-related uncertainties and down
  bool isJECUncertaintyDown(const sysType type);

  // true if type is one of the JEC-related uncertainties
  bool isJECUncertainty(const sysType type);

  // return the label that is used by JetCorrectorParametersCollection
  // to label the uncertainty type. See also:
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_23/doc/html/dc/d33/classJetCorrectorParametersCollection.html#afb3d4c6fd711ca23d89e0625a22dc483 for a list of in principle valid labels. Whether the uncertainty
  std::string GetJECUncertaintyLabel(const sysType type);
};


#endif
