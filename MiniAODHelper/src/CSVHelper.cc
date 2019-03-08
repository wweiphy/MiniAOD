#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include <algorithm>
#include "FWCore/Utilities/interface/Exception.h"

#include "MiniAOD/MiniAODHelper/interface/CSVHelper.h"

CSVHelper::CSVHelper()
  : isInit_(false), nHFptBins_(0),nLFptBins_(0),nLFetaBins_(0), allowJetsOutOfBinning_(false) {}


CSVHelper::CSVHelper(const std::string& hf, const std::string& lf, const int& nHFptBins,const int& nLFptBins,const int& nLFetaBins, const std::vector<Systematics::Type>& jecsysts)
  : isInit_(false), nHFptBins_(0),nLFptBins_(0),nLFetaBins_(0), allowJetsOutOfBinning_(false) {
  init(hf,lf,nHFptBins,nLFptBins,nLFetaBins,jecsysts);
}


CSVHelper::~CSVHelper() {
  for(auto& i: h_csv_wgt_hf ) {
    for(auto& j: i) {
      if( j ) delete j;
    }
  }
  for(auto& i: c_csv_wgt_hf ) {
    for(auto& j: i) {
      if( j ) delete j;
    }
  }
  for(auto& i: h_csv_wgt_lf ) {
    for(auto& j: i) {
      for(auto& k: j) {
	if( k ) delete k;
      }
    }
  }
}


void CSVHelper::init(const std::string& hf, const std::string& lf, const int& nHFptBins,const int& nLFptBins,const int& nLFetaBins,const std::vector<Systematics::Type>& jecsysts) {
  std::cout << "Initializing b-tag scale factors"
	    << "\n  HF : " << hf << " (" << nHFptBins << " pt bins)"
	    << "\n  LF : " << lf << " (" << nLFptBins << " pt bins)" 
            << "\n  LF : " << lf << " (" << nLFetaBins << " eta bins)" <<  std::endl;

  nHFptBins_ = nHFptBins;
  nLFptBins_ = nLFptBins;
  nLFetaBins_ = nLFetaBins;
  
  //combine the vector with the csv systematics and the jec systematics into one vector
  systs.reserve(csvsysts.size()+jecsysts.size());
  systs.insert(systs.end(),jecsysts.begin(),jecsysts.end());
  systs.insert(systs.end(),csvsysts.begin(),csvsysts.end());
  
  const std::string inputFileHF = hf.size() > 0 ? hf : "data/csv_rwt_hf_IT_FlatSF.root";
  const std::string inputFileLF = lf.size() > 0 ? lf : "data/csv_rwt_lf_IT_FlatSF.root";

  TFile *f_CSVwgt_HF = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileHF).c_str());
  TFile *f_CSVwgt_LF = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/" + inputFileLF).c_str());
  fillCSVHistos(f_CSVwgt_HF, f_CSVwgt_LF,systs);
  f_CSVwgt_HF->Close();
  f_CSVwgt_LF->Close();
  delete f_CSVwgt_HF;
  delete f_CSVwgt_LF;

  isInit_ = true;
}

// fill the histograms (done once)
void
CSVHelper::fillCSVHistos(TFile *fileHF, TFile *fileLF, const std::vector<Systematics::Type>& systs)
{
  const size_t nSys = systs.size();//purity,stats1,stats2 with up/down + jec systs including nominal variation
  h_csv_wgt_hf = std::vector< std::vector<TH1*> >(nSys,std::vector<TH1*>(nHFptBins_,NULL));
  c_csv_wgt_hf = std::vector< std::vector<TH1*> >(nSys,std::vector<TH1*>(nHFptBins_,NULL));
  h_csv_wgt_lf = std::vector< std::vector< std::vector<TH1*> > >(nSys,std::vector< std::vector<TH1*> >(nLFptBins_,std::vector<TH1*>(nLFetaBins_,NULL)));
  TString syst_csv_suffix = "final";
  // loop over all the available systematics
  for (size_t iSys = 0; iSys < nSys; iSys++) {
    // some string cosmetics to search for the correct histogram in the root files  
    TString systematic = Systematics::toString(systs[iSys]);
    TString systematic_original = systematic;
    
    systematic.ReplaceAll("up","Up");
    systematic.ReplaceAll("down","Down");
    systematic.ReplaceAll("CSV","");
    std::cout << "adding histograms for systematic " << systematic << std::endl;
    if(systematic.Contains("Stats")) {
        systematic.ReplaceAll("HF","");
        systematic.ReplaceAll("LF","");
    }
    if(systematic!="") {systematic="_"+systematic;}
    
    if(systematic_original.Contains("HFStats")){
        systematic_original.ReplaceAll("HFStats","LFStats");
    }
    else if(systematic_original.Contains("LFStats")){
        systematic_original.ReplaceAll("LFStats","HFStats");
    }
    // loop over all pt bins of the different jet flavours
    for (int iPt = 0; iPt < nHFptBins_; iPt++) {
        TString name = Form("csv_ratio_Pt%i_Eta0_%s", iPt, (syst_csv_suffix+systematic).Data());
        // only read the histogram if it exits in the root file
        if(fileHF->GetListOfKeys()->Contains(name)) {
            h_csv_wgt_hf.at(iSys).at(iPt) = readHistogram(fileHF,name);
            std::cout <<"for "<<systematic_original<< " added " <<  h_csv_wgt_hf.at(iSys).at(iPt)->GetName() << " from HF file" << std::endl;
        }
        else {
            h_csv_wgt_hf.at(iSys).at(iPt) = readHistogram(fileHF,name.ReplaceAll(systematic,""));
            std::cout <<"for "<<systematic_original<< " added " <<  h_csv_wgt_hf.at(iSys).at(iPt)->GetName() << " from HF file" << std::endl;
        }
    }
    for (int iPt = 0; iPt < nHFptBins_; iPt++) {
        TString name = Form("c_csv_ratio_Pt%i_Eta0_%s", iPt, (syst_csv_suffix+systematic).Data());
        if(fileHF->GetListOfKeys()->Contains(name)) {
            c_csv_wgt_hf.at(iSys).at(iPt) = readHistogram(fileHF,name);
            std::cout <<"for "<<systematic_original<< " added " << c_csv_wgt_hf.at(iSys).at(iPt)->GetName() << " from CF(HF) file" << std::endl;
        }
        else {
            c_csv_wgt_hf.at(iSys).at(iPt) = readHistogram(fileHF,name.ReplaceAll(systematic,""));
            std::cout <<"for "<<systematic_original<< " added " << c_csv_wgt_hf.at(iSys).at(iPt)->GetName() << " from CF(HF) file" << std::endl;
        }
    }
    for (int iPt = 0; iPt < nLFptBins_; iPt++) {
        for (int iEta = 0; iEta < nLFetaBins_; iEta++) {
            TString name = Form("csv_ratio_Pt%i_Eta%i_%s", iPt, iEta, (syst_csv_suffix+systematic).Data());
            if(fileLF->GetListOfKeys()->Contains(name)) {
                h_csv_wgt_lf.at(iSys).at(iPt).at(iEta) = readHistogram(fileLF,name);
                std::cout <<"for "<<systematic_original<< " added " << h_csv_wgt_lf.at(iSys).at(iPt).at(iEta)->GetName() << " from LF file" << std::endl;
            }
            else {
                h_csv_wgt_lf.at(iSys).at(iPt).at(iEta) = readHistogram(fileLF,name.ReplaceAll(systematic,""));
                std::cout <<"for "<<systematic_original<< " added " << h_csv_wgt_lf.at(iSys).at(iPt).at(iEta)->GetName() << " from LF file" << std::endl;
            }
        }
    }
  }
}


TH1* CSVHelper::readHistogram(TFile* file, const TString& name) const {
  TH1* h = NULL;
  file->GetObject(name,h);
  if( h==NULL ) {
    throw cms::Exception("BadCSVWeightInit")
      << "Could not find CSV SF histogram '" << name
      << "' in file '" << file->GetName() << "'";
  }
  h->SetDirectory(0);
  
  return h;
}


double
CSVHelper::getCSVWeight(const std::vector<double>& jetPts,
			const std::vector<double>& jetEtas,
			const std::vector<double>& jetCSVs,
			const std::vector<int>& jetFlavors,
			const Systematics::Type syst,
			double &csvWgtHF,
			double &csvWgtLF,
			double &csvWgtCF) const
{
  if( !isInit_ ) {
    throw cms::Exception("BadCSVWeightAccess") << "CSVHelper not initialized";
  }
  // search for the position of the desired systematic in the systs vector
  const int iSys = std::find(systs.begin(),systs.end(),syst)-systs.begin();
  
  //std::cout << "Systematic index " << iSys << std::endl;
  // initialize the weight for the different jet flavours with 1
  double csvWgthf = 1.;
  double csvWgtC = 1.;
  double csvWgtlf = 1.;
  
  // loop over all jets in the event and calculate the final weight by multiplying the single jet scale factors
  for (size_t iJet = 0; iJet < jetPts.size(); iJet++) {
    const double csv = jetCSVs.at(iJet);
    const double jetPt = jetPts.at(iJet);
    const double jetAbsEta = fabs(jetEtas.at(iJet));
    const int flavor = jetFlavors.at(iJet);

    int iPt = -1;
    int iEta = -1;
    // pt binning for heavy flavour jets
    if(abs(flavor)>3) {
        if (jetPt >= 19.99 && jetPt < 30)
            iPt = 0;
        else if (jetPt >= 30 && jetPt < 50)
            iPt = 1;
        else if (jetPt >= 50 && jetPt < 70)
            iPt = 2;
        else if (jetPt >= 70 && jetPt < 100)
            iPt = 3;
        else if (jetPt >= 100 && jetPt < 160)
            iPt = 4;
        else if (jetPt >= 160)
            iPt = 5;
    }
    // pt binning for light flavour jets
    else {
        if (jetPt >= 19.99 && jetPt < 30)
            iPt = 0;
        else if (jetPt >= 30 && jetPt < 40)
            iPt = 1;
        else if (jetPt >= 40 && jetPt < 60)
            iPt = 2;
        else if (jetPt >= 60 && jetPt < 100)
            iPt = 3;
        else if (jetPt >= 100 && jetPt < 160)
            iPt = 4;
        else if (jetPt >= 160)
            iPt = 5;
    }
    // light flavour jets also have eta bins
    if (jetAbsEta >= 0 && jetAbsEta < 0.8)
      iEta = 0;
    else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6)
      iEta = 1;
    else if (jetAbsEta >= 1.6 && jetAbsEta < 2.5)
      iEta = 2;
    
    if (iPt < 0 || iEta < 0) {
      if( allowJetsOutOfBinning_ ) continue;
      throw cms::Exception("BadCSVWeightAccess") << "couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta;
    }
    
    //std::cout << "program is in front of calculating the csv weights " << std::endl;
    // b flavour jet
    if (abs(flavor) == 5) {
        //std::cout << "b flavor jet " << std::endl;
      // RESET iPt to maximum pt bin (only 5 bins for new SFs)
      if(iPt>=nHFptBins_){
	iPt=nHFptBins_-1;// [20-30], [30-50], [40-70], [70,100] and [1000-10000] only 5 Pt bins for hf
      }
      if(h_csv_wgt_hf.at(iSys).at(iPt)) {
        const int useCSVBin = (csv >= 0.) ? h_csv_wgt_hf.at(iSys).at(iPt)->FindBin(csv) : 1;
        const double iCSVWgtHF = h_csv_wgt_hf.at(iSys).at(iPt)->GetBinContent(useCSVBin);
        if (iCSVWgtHF != 0) csvWgthf *= iCSVWgtHF;
      }
    } // c flavour jet
    else if (abs(flavor) == 4) {
        //std::cout << "c flavor jet " << std::endl;
      // RESET iPt to maximum pt bin (only 5 bins for new SFs)
      if(iPt>=nHFptBins_){
	iPt=nHFptBins_-1;// [20-30], [30-50], [40-70], [70,100] and [1000-10000] only 5 Pt bins for cf
      }
      if(c_csv_wgt_hf.at(iSys).at(iPt)) {
        const int useCSVBin = (csv >= 0.) ? c_csv_wgt_hf.at(iSys).at(iPt)->FindBin(csv) : 1;
        const double iCSVWgtC = c_csv_wgt_hf.at(iSys).at(iPt)->GetBinContent(useCSVBin);
        if (iCSVWgtC != 0) csvWgtC *= iCSVWgtC;
      }
    } // light flavour jet
    else {
        //std::cout << "light flavor jet " << std::endl;
      if (iPt >= nLFptBins_) iPt = nLFptBins_-1; // [20-30], [30-40], [40-60] and [60-10000] only 4 Pt bins for lf
      if(h_csv_wgt_lf.at(iSys).at(iPt).at(iEta)) {
        const int useCSVBin = (csv >= 0.) ? h_csv_wgt_lf.at(iSys).at(iPt).at(iEta)->FindBin(csv) : 1;
        const double iCSVWgtLF = h_csv_wgt_lf.at(iSys).at(iPt).at(iEta)->GetBinContent(useCSVBin);
        if (iCSVWgtLF != 0) csvWgtlf *= iCSVWgtLF;
      }
    }
  }

  const double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}


float CSVHelper::GetJetCSV(const pat::Jet& jet, const std::string taggername){

  float defaultFailure = -.1;
  float bTagVal=0;
  if(taggername=="DeepCSV"){
    bTagVal=jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
  }
  else if(taggername=="DeepJet"){
    bTagVal=jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");    
  }
  else if (taggername == "CSVv2" or TString(taggername).BeginsWith("pfDeep")){
    bTagVal = jet.bDiscriminator(taggername);
  }
  else{
      throw cms::Exception("CSVHelper: Invalid taggername ") << "Taggername" << taggername << "not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
      bTagVal = defaultFailure;
  }

  if(isnan(bTagVal)) return defaultFailure;

  if(bTagVal > 1.) return 1.;
  if(bTagVal < 0.) return defaultFailure;

  return bTagVal;
}

float CSVHelper::GetJetCSV_DNN(const pat::Jet& jet, const std::string taggername){


  float bTagVal=0;
  if(taggername=="DeepCSV"){
    bTagVal=jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
  }
  else if(taggername=="DeepJet"){
    bTagVal=jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");    
  }
  else if (taggername == "CSVv2" or TString(taggername).BeginsWith("pfDeep")){
    bTagVal = jet.bDiscriminator(taggername);
  }
  else{
      throw cms::Exception("CSVHelper: Invalid taggername ") << "Taggername" << taggername << "not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
  }
  return bTagVal;
}

bool CSVHelper::PassesCSV(std::string dataEra, const pat::Jet& iJet, std::string taggername, const CSVwp iCSVworkingPoint){
  float csvValue = GetJetCSV(iJet, taggername);
  // CSV b-tagging requirement
  if(csvValue > GetWP(dataEra, iCSVworkingPoint, taggername)){
     return true;
  }
  else {
    return false;
  }
}

float CSVHelper::GetWP(std::string dataEra, const CSVwp iCSVworkingPoint, std::string taggername){
  if (TString(dataEra).Contains("2016")){
    if (taggername == "DeepCSV"){
      switch(iCSVworkingPoint){
      case CSVwp::Loose:	  { return 0.2217; }	break;
      case CSVwp::Medium: 	{ return 0.6321; }	break;
      case CSVwp::Tight:	  { return 0.8953; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else if (taggername == "DeepJet"){
      switch(iCSVworkingPoint){
      case CSVwp::Loose:	  { return 0.0614; }	break;
      case CSVwp::Medium: 	{ return 0.3093; }	break;
      case CSVwp::Tight:	  { return 7221; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else if (taggername == "CSVv2"){
      switch(iCSVworkingPoint){ // CSVv2 not supported for 2016 Legacy->WP are assumed to be the same ones as for 2017 by me
      case CSVwp::Loose:	  { return 0.5803; }	break;
      case CSVwp::Medium: 	{ return 0.8838; }	break;
      case CSVwp::Tight:	  { return 0.9693; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else {
      throw cms::Exception("CSVHelper: Invalid taggername ") << "Taggername" << taggername << "not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
      return 0;
    }
    return 0;
  }
  else if (TString(dataEra).Contains("2017")){
    if (taggername == "DeepCSV"){
      switch(iCSVworkingPoint){
      case CSVwp::Loose:	  { return 0.1522; }	break;
      case CSVwp::Medium: 	{ return 0.4941; }	break;
      case CSVwp::Tight:	  { return 0.8001; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else if (taggername == "DeepJet"){
      switch(iCSVworkingPoint){
      case CSVwp::Loose:	  { return 0.0521; }	break;
      case CSVwp::Medium: 	{ return 0.3033; }	break;
      case CSVwp::Tight:	  { return 0.7489; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else if (taggername == "CSVv2"){
      switch(iCSVworkingPoint){
      case CSVwp::Loose:	  { return 0.5803; }	break;
      case CSVwp::Medium: 	{ return 0.8838; }	break;
      case CSVwp::Tight:	  { return 0.9693; }	break;
      case CSVwp::None:	    return 0;  
      }
    }
    else {
      throw cms::Exception("CSVHelper: Invalid taggername ") << "Taggername" << taggername << "not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
      return 0;
    }
  }
  else {
    throw cms::Exception("CSVHelper: Invalid dataEra ") << "dataEra" << dataEra << "not recognized, only 2016/2017 data possible" << std::endl;
    return 0;
  }
  return 0;
}
