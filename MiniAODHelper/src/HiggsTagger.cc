#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"


HiggsTagger::HiggsTagger(){
}


HiggsTagger::HiggsTagger(MiniAODHelper* helper_, HiggsTag::Mode mode_, std::string filePath_):helper(helper_),mode(mode_),btagger("pfCombinedInclusiveSecondaryVertexV2BJetTags"){
  
  char* CMSSWPath = getenv("CMSSW_BASE");
  std::string filePath = CMSSWPath;
  filePath = filePath+"/src/BoostedTTH/BoostedAnalyzer/data/toptagger/"+filePath_;
  
  switch(mode){
    case HiggsTag::SecondCSV:
      break;
      
    case HiggsTag::TMVA:
      {  
        // Get absolute path for TMVA weights
        std::string weightsPath = filePath;

        // Initialize TMVA input variables
        GetTMVAVars(weightsPath);

        // Setup TMVA Reader
        TMVAReader = new TMVA::Reader("Silent");
        
        for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          TMVAReader->AddVariable(*itVarName,&TMVAVars[itVarName-TMVAVarNames.begin()]);
        }

        TMVAReader->BookMVA("TMVAHiggsTagger",(weightsPath).c_str());
      }      
      break;

    default:
      std::cout << "Error! No matching Higgs Tagger Mode found!" << std:: endl;
  }     
}


HiggsTagger::~HiggsTagger(){
}


void HiggsTagger::GetTMVAVarNames(std::string filePath_, bool verbose){
  
  TMVAVarNames.clear();
  
  std::string line;
  ifstream inFile (filePath_);
  
  if(inFile.is_open()){
  
    while(getline(inFile,line)){
      if(line.find("</Variables>")!=std::string::npos) break;
      
      if(line.find("<Variable VarIndex")!=std::string::npos){
        size_t varStart = line.find("\"",line.find("Title"))+1;
        size_t varEnd = line.find("\"",varStart);

        TMVAVarNames.push_back(line.substr(varStart,varEnd-varStart));
      }
    }
    
    if(verbose){
      for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
        std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << std::endl;
      }
    }
    
    inFile.close();
  }
}


void HiggsTagger::GetTMVAVars(std::string filePath_, bool verbose){
  
  GetTMVAVarNames(filePath_,verbose);
  
  TMVAVars = new float[50];
  
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
  
  if(verbose){
    for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
      std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
    }
  }  
}


void HiggsTagger::ResetTMVAVars(){
  for(int i=0;i<50;i++){
    TMVAVars[i] = -999;
  }
}


float HiggsTagger::GetHiggsTaggerOutput(const boosted::SubFilterJet& higgsJet, bool verbose){
  
  float failReturn = -1.1;
  if(mode == HiggsTag::SecondCSV){
    failReturn = 0.1;
  }
  
  if(higgsJet.filterjets.size()<2) return failReturn;
  
  std::vector<pat::Jet> filterjets = higgsJet.filterjets;
  filterjets = helper->GetSortedByCSV(filterjets);
  
  switch(mode){
    
    case HiggsTag::SecondCSV:
      {
        return fmax(filterjets[1].bDiscriminator(btagger),-.1);
      }
      break;
      
    case HiggsTag::TMVA:
      {
        ResetTMVAVars();

        double M2 = (filterjets[0].p4()+filterjets[1].p4()).M();
        double M3 = M2;
        if(higgsJet.filterjets.size()>=3)
          M3 = (filterjets[0].p4()+filterjets[1].p4()+filterjets[2].p4()).M();

        for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
          int iVar = itVarName-TMVAVarNames.begin();

          if(*itVarName=="HiggsJet_Pt")                           TMVAVars[iVar] = higgsJet.fatjet.pt();                                                                                                               
          else if(*itVarName=="HiggsJet_M2")                      TMVAVars[iVar] = M2;                                                                                                                                 
          else if(*itVarName=="HiggsJet_M3")                      TMVAVars[iVar] = M3;                                                                                                                                 
          else if(*itVarName=="HiggsJet_CSV1")                    TMVAVars[iVar] = filterjets[0].bDiscriminator(btagger);                                                                                              
          else if(*itVarName=="HiggsJet_CSV2")                    TMVAVars[iVar] = filterjets[1].bDiscriminator(btagger);                                                                                              
          else if(*itVarName=="HiggsJet_NSubjettiness_12_Ratio")  TMVAVars[iVar] = higgsJet.subjettiness2/higgsJet.subjettiness1;                                                                                      
          else if(*itVarName=="HiggsJet_NSubjettiness_23_Ratio")  TMVAVars[iVar] = higgsJet.subjettiness3/higgsJet.subjettiness2;                                                                                      
          else std::cout << "Error! No matching Top Tagger Input Variable found!" << std:: endl;
        }

        if(verbose){
          std::cout << "Higgs Tagger Variables:" << std::endl;
          for(std::vector<string>::const_iterator itVarName=TMVAVarNames.begin();itVarName!=TMVAVarNames.end();++itVarName){
            std::cout << "Variable " << itVarName-TMVAVarNames.begin() << ": " << *itVarName << ": " << TMVAVars[itVarName-TMVAVarNames.begin()] << std::endl;
          }
        }

        float TMVAoutput = TMVAReader->EvaluateMVA("TMVATopTagger");

        if(verbose) std::cout << "TMVAOutput: " << TMVAoutput << std::endl;

        return TMVAoutput;
      }
      break;
    
    default:
    
      std::cout << "Error! No matching Higgs Tagger Mode found!" << std:: endl;
      return -999;
  }
}    


boosted::SubFilterJetCollection HiggsTagger::GetSortedByHiggsTaggerOutput(const boosted::SubFilterJetCollection& higgsJets, bool verbose){
  
  boosted::SubFilterJetCollection result = higgsJets;
  
  HiggsTaggerOutputComparison higgsTagComp(this);
  std::sort(result.begin(), result.end(),higgsTagComp);
  
  return result;
}


float HiggsTagger::GetHiggsCand(boosted::SubFilterJetCollection& higgsJets, boosted::SubFilterJet& higgsCand, bool verbose){
  
  float maxHiggsTag=-9999;
  
  for(boosted::SubFilterJetCollection::iterator itHiggsJet = higgsJets.begin() ; itHiggsJet != higgsJets.end(); ++itHiggsJet){
    
    float higgsTag = GetHiggsTaggerOutput(*itHiggsJet,verbose);
    
    if(higgsTag > maxHiggsTag) {
      maxHiggsTag = higgsTag;
      higgsCand = *itHiggsJet;
    }
  } 
  
  return maxHiggsTag;
}
