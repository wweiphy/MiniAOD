#include "MiniAOD/MiniAODHelper/interface/Systematics.h"
#include "FWCore/Utilities/interface/Exception.h"

sysType::sysType sysType::get(const std::string& name) {
  if( name == "" ) return sysType::NA;

  if( name == "JERup"   ) return sysType::JERup;
  if( name == "JERdown" ) return sysType::JERdown;

  if( name == "JESup"                ) return JESup;
  if( name == "JESPileUpDataMCup"    ) return JESPileUpDataMCup;   
  if( name == "JESPileUpPtRefup"     ) return JESPileUpPtRefup;      
  if( name == "JESPileUpPtBBup"      ) return JESPileUpPtBBup; 	 
  if( name == "JESPileUpPtEC1up"     ) return JESPileUpPtEC1up;    
  if( name == "JESPileUpPtEC2up"     ) return JESPileUpPtEC2up;    
  if( name == "JESPileUpPtHFup"      ) return JESPileUpPtHFup;      	 
  if( name == "JESRelativeJEREC1up"  ) return JESRelativeJEREC1up; 
  if( name == "JESRelativeJEREC2up"  ) return JESRelativeJEREC2up; 
  if( name == "JESRelativeJERHFup"   ) return JESRelativeJERHFup;  
  if( name == "JESRelativeFSRup"     ) return JESRelativeFSRup;    
  if( name == "JESRelativeStatFSRup" ) return JESRelativeStatFSRup;
  if( name == "JESRelativeStatEC2up" ) return JESRelativeStatEC2up;
  if( name == "JESRelativeStatECup"  ) return JESRelativeStatECup; 
  if( name == "JESRelativeStatHFup"  ) return JESRelativeStatHFup; 
  if( name == "JESRelativePtBBup"    ) return JESRelativePtBBup;   
  if( name == "JESRelativePtEC1up"   ) return JESRelativePtEC1up;  
  if( name == "JESRelativePtEC2up"   ) return JESRelativePtEC2up;  
  if( name == "JESRelativePtHFup"    ) return JESRelativePtHFup;   
  if( name == "JESTimeEtaup"         ) return JESTimeEtaup;         	 
  if( name == "JESAbsoluteScaleup"   ) return JESAbsoluteScaleup;  
  if( name == "JESAbsoluteMPFBiasup" ) return JESAbsoluteMPFBiasup;
  if( name == "JESAbsoluteStatup"    ) return JESAbsoluteStatup;   
  if( name == "JESSinglePionECALup"  ) return JESSinglePionECALup; 
  if( name == "JESSinglePionHCALup"  ) return JESSinglePionHCALup; 
  if( name == "JESFragmentationup"   ) return JESFragmentationup;  
  if( name == "JESTimePtup"          ) return JESTimePtup;          	    	 
  if( name == "JESFlavorQCDup"       ) return JESFlavorQCDup;

  if( name == "JESdown"                ) return JESdown;
  if( name == "JESPileUpDataMCdown"    ) return JESPileUpDataMCdown;   
  if( name == "JESPileUpPtRefdown"     ) return JESPileUpPtRefdown;      
  if( name == "JESPileUpPtBBdown"      ) return JESPileUpPtBBdown; 	 
  if( name == "JESPileUpPtEC1down"     ) return JESPileUpPtEC1down;    
  if( name == "JESPileUpPtEC2down"     ) return JESPileUpPtEC2down;    
  if( name == "JESPileUpPtHFdown"      ) return JESPileUpPtHFdown;      	 
  if( name == "JESRelativeJEREC1down"  ) return JESRelativeJEREC1down; 
  if( name == "JESRelativeJEREC2down"  ) return JESRelativeJEREC2down; 
  if( name == "JESRelativeJERHFdown"   ) return JESRelativeJERHFdown;  
  if( name == "JESRelativeFSRdown"     ) return JESRelativeFSRdown;    
  if( name == "JESRelativeStatFSRdown" ) return JESRelativeStatFSRdown;
  if( name == "JESRelativeStatEC2down" ) return JESRelativeStatEC2down;
  if( name == "JESRelativeStatECdown"  ) return JESRelativeStatECdown; 
  if( name == "JESRelativeStatHFdown"  ) return JESRelativeStatHFdown; 
  if( name == "JESRelativePtBBdown"    ) return JESRelativePtBBdown;   
  if( name == "JESRelativePtEC1down"   ) return JESRelativePtEC1down;  
  if( name == "JESRelativePtEC2down"   ) return JESRelativePtEC2down;  
  if( name == "JESRelativePtHFdown"    ) return JESRelativePtHFdown;   
  if( name == "JESTimeEtadown"         ) return JESTimeEtadown;         	 
  if( name == "JESAbsoluteScaledown"   ) return JESAbsoluteScaledown;  
  if( name == "JESAbsoluteMPFBiasdown" ) return JESAbsoluteMPFBiasdown;
  if( name == "JESAbsoluteStatdown"    ) return JESAbsoluteStatdown;   
  if( name == "JESSinglePionECALdown"  ) return JESSinglePionECALdown; 
  if( name == "JESSinglePionHCALdown"  ) return JESSinglePionHCALdown; 
  if( name == "JESFragmentationdown"   ) return JESFragmentationdown;  
  if( name == "JESTimePtdown"          ) return JESTimePtdown;          	    	 
  if( name == "JESFlavorQCDdown"       ) return JESFlavorQCDdown;

  throw cms::Exception("InvalidUncertaintyName") << "No uncertainty with name '" << name << "'";
  return sysType::NA;
}


std::string sysType::toString(const sysType type) {
  if( type == NA ) return "";

  if( type == JERup   ) return "JERup";
  if( type == JERdown ) return "JERdown";

  if( type == JESup                ) return "JESup";
  if( type == JESPileUpDataMCup    ) return "JESPileUpDataMCup";   
  if( type == JESPileUpPtRefup     ) return "JESPileUpPtRefup";      
  if( type == JESPileUpPtBBup      ) return "JESPileUpPtBBup"; 	 
  if( type == JESPileUpPtEC1up     ) return "JESPileUpPtEC1up";    
  if( type == JESPileUpPtEC2up     ) return "JESPileUpPtEC2up";    
  if( type == JESPileUpPtHFup      ) return "JESPileUpPtHFup";      	 
  if( type == JESRelativeJEREC1up  ) return "JESRelativeJEREC1up"; 
  if( type == JESRelativeJEREC2up  ) return "JESRelativeJEREC2up"; 
  if( type == JESRelativeJERHFup   ) return "JESRelativeJERHFup";  
  if( type == JESRelativeFSRup     ) return "JESRelativeFSRup";    
  if( type == JESRelativeStatFSRup ) return "JESRelativeStatFSRup";
  if( type == JESRelativeStatEC2up ) return "JESRelativeStatEC2up";
  if( type == JESRelativeStatECup  ) return "JESRelativeStatECup"; 
  if( type == JESRelativeStatHFup  ) return "JESRelativeStatHFup"; 
  if( type == JESRelativePtBBup    ) return "JESRelativePtBBup";   
  if( type == JESRelativePtEC1up   ) return "JESRelativePtEC1up";  
  if( type == JESRelativePtEC2up   ) return "JESRelativePtEC2up";  
  if( type == JESRelativePtHFup    ) return "JESRelativePtHFup";   
  if( type == JESTimeEtaup         ) return "JESTimeEtaup";         	 
  if( type == JESAbsoluteScaleup   ) return "JESAbsoluteScaleup";  
  if( type == JESAbsoluteMPFBiasup ) return "JESAbsoluteMPFBiasup";
  if( type == JESAbsoluteStatup    ) return "JESAbsoluteStatup";   
  if( type == JESSinglePionECALup  ) return "JESSinglePionECALup"; 
  if( type == JESSinglePionHCALup  ) return "JESSinglePionHCALup"; 
  if( type == JESFragmentationup   ) return "JESFragmentationup";  
  if( type == JESTimePtup          ) return "JESTimePtup";          	    	 
  if( type == JESFlavorQCDup       ) return "JESFlavorQCDup";

  if( type == JESdown                ) return "JESdown";
  if( type == JESPileUpDataMCdown    ) return "JESPileUpDataMCdown";   
  if( type == JESPileUpPtRefdown     ) return "JESPileUpPtRefdown";      
  if( type == JESPileUpPtBBdown      ) return "JESPileUpPtBBdown"; 	 
  if( type == JESPileUpPtEC1down     ) return "JESPileUpPtEC1down";    
  if( type == JESPileUpPtEC2down     ) return "JESPileUpPtEC2down";    
  if( type == JESPileUpPtHFdown      ) return "JESPileUpPtHFdown";      	 
  if( type == JESRelativeJEREC1down  ) return "JESRelativeJEREC1down"; 
  if( type == JESRelativeJEREC2down  ) return "JESRelativeJEREC2down"; 
  if( type == JESRelativeJERHFdown   ) return "JESRelativeJERHFdown";  
  if( type == JESRelativeFSRdown     ) return "JESRelativeFSRdown";    
  if( type == JESRelativeStatFSRdown ) return "JESRelativeStatFSRdown";
  if( type == JESRelativeStatEC2down ) return "JESRelativeStatEC2down";
  if( type == JESRelativeStatECdown  ) return "JESRelativeStatECdown"; 
  if( type == JESRelativeStatHFdown  ) return "JESRelativeStatHFdown"; 
  if( type == JESRelativePtBBdown    ) return "JESRelativePtBBdown";   
  if( type == JESRelativePtEC1down   ) return "JESRelativePtEC1down";  
  if( type == JESRelativePtEC2down   ) return "JESRelativePtEC2down";  
  if( type == JESRelativePtHFdown    ) return "JESRelativePtHFdown";   
  if( type == JESTimeEtadown         ) return "JESTimeEtadown";         	 
  if( type == JESAbsoluteScaledown   ) return "JESAbsoluteScaledown";  
  if( type == JESAbsoluteMPFBiasdown ) return "JESAbsoluteMPFBiasdown";
  if( type == JESAbsoluteStatdown    ) return "JESAbsoluteStatdown";   
  if( type == JESSinglePionECALdown  ) return "JESSinglePionECALdown"; 
  if( type == JESSinglePionHCALdown  ) return "JESSinglePionHCALdown"; 
  if( type == JESFragmentationdown   ) return "JESFragmentationdown";  
  if( type == JESTimePtdown          ) return "JESTimePtdown";          	    	 
  if( type == JESFlavorQCDdown       ) return "JESFlavorQCDdown";

  throw cms::Exception("InvalidUncertaintyType") << "No uncertainty with index '" << type << "'";
  return "";
}


bool sysType::isJECUncertaintyUp(const sysType type) {
  if( type == JESup ) return true;

  if( type == JESPileUpDataMCup    ) return true;		
  if( type == JESPileUpPtRefup     ) return true;
  if( type == JESPileUpPtBBup      ) return true;			
  if( type == JESPileUpPtEC1up     ) return true;
  if( type == JESPileUpPtEC2up     ) return true;
  if( type == JESPileUpPtHFup      ) return true;
  if( type == JESRelativeJEREC1up  ) return true;
  if( type == JESRelativeJEREC2up  ) return true;		
  if( type == JESRelativeJERHFup   ) return true;
  if( type == JESRelativeFSRup     ) return true;
  if( type == JESRelativeStatFSRup ) return true;
  if( type == JESRelativeStatEC2up ) return true;
  if( type == JESRelativeStatECup  ) return true;		
  if( type == JESRelativeStatHFup  ) return true;
  if( type == JESRelativePtBBup    ) return true;
  if( type == JESRelativePtEC1up   ) return true;
  if( type == JESRelativePtEC2up   ) return true;
  if( type == JESRelativePtHFup    ) return true;		
  if( type == JESTimeEtaup         ) return true;
  if( type == JESAbsoluteScaleup   ) return true;
  if( type == JESAbsoluteMPFBiasup ) return true;
  if( type == JESAbsoluteStatup    ) return true;
  if( type == JESSinglePionECALup  ) return true;		
  if( type == JESSinglePionHCALup  ) return true;
  if( type == JESFragmentationup   ) return true;
  if( type == JESTimePtup          ) return true;
  if( type == JESFlavorQCDup       ) return true;			

  else return false;
}

bool sysType::isJECUncertaintyDown(const sysType type) {
  if( type == JESdown ) return true;

  if( type == JESPileUpDataMCdown    ) return true;		
  if( type == JESPileUpPtRefdown     ) return true;
  if( type == JESPileUpPtBBdown      ) return true;			
  if( type == JESPileUpPtEC1down     ) return true;
  if( type == JESPileUpPtEC2down     ) return true;
  if( type == JESPileUpPtHFdown      ) return true;
  if( type == JESRelativeJEREC1down  ) return true;
  if( type == JESRelativeJEREC2down  ) return true;		
  if( type == JESRelativeJERHFdown   ) return true;
  if( type == JESRelativeFSRdown     ) return true;
  if( type == JESRelativeStatFSRdown ) return true;
  if( type == JESRelativeStatEC2down ) return true;
  if( type == JESRelativeStatECdown  ) return true;		
  if( type == JESRelativeStatHFdown  ) return true;
  if( type == JESRelativePtBBdown    ) return true;
  if( type == JESRelativePtEC1down   ) return true;
  if( type == JESRelativePtEC2down   ) return true;
  if( type == JESRelativePtHFdown    ) return true;		
  if( type == JESTimeEtadown         ) return true;
  if( type == JESAbsoluteScaledown   ) return true;
  if( type == JESAbsoluteMPFBiasdown ) return true;
  if( type == JESAbsoluteStatdown    ) return true;
  if( type == JESSinglePionECALdown  ) return true;		
  if( type == JESSinglePionHCALdown  ) return true;
  if( type == JESFragmentationdown   ) return true;
  if( type == JESTimePtdown          ) return true;
  if( type == JESFlavorQCDdown       ) return true;			

  else return false;
}

bool sysType::isJECUncertainty(const sysType type) {
  return isJECUncertaintyUp(type) || isJECUncertaintyDown(type);
}

std::string sysType::GetJECUncertaintyLabel(const sysType type) {
  if( type == JESup                || type == JESdown                ) return "Uncertainty";
  if( type == JESPileUpDataMCup    || type == JESPileUpDataMCdown    ) return "PileUpDataMC";   
  if( type == JESPileUpPtRefup     || type == JESPileUpPtRefdown     ) return "PileUpPtRef";      
  if( type == JESPileUpPtBBup      || type == JESPileUpPtBBdown      ) return "PileUpPtBB";	 
  if( type == JESPileUpPtEC1up     || type == JESPileUpPtEC1down     ) return "PileUpPtEC1";    
  if( type == JESPileUpPtEC2up     || type == JESPileUpPtEC2down     ) return "PileUpPtEC2";    
  if( type == JESPileUpPtHFup      || type == JESPileUpPtHFdown      ) return "PileUpPtHF";	 
  if( type == JESRelativeJEREC1up  || type == JESRelativeJEREC1down  ) return "RelativeJEREC1"; 
  if( type == JESRelativeJEREC2up  || type == JESRelativeJEREC2down  ) return "RelativeJEREC2"; 
  if( type == JESRelativeJERHFup   || type == JESRelativeJERHFdown   ) return "RelativeJERHF";  
  if( type == JESRelativeFSRup     || type == JESRelativeFSRdown     ) return "RelativeFSR";    
  if( type == JESRelativeStatFSRup || type == JESRelativeStatFSRdown ) return "RelativeStatFSR";
  if( type == JESRelativeStatEC2up || type == JESRelativeStatEC2down ) return "RelativeStatEC2";
  if( type == JESRelativeStatECup  || type == JESRelativeStatECdown  ) return "RelativeStatEC"; 
  if( type == JESRelativeStatHFup  || type == JESRelativeStatHFdown  ) return "RelativeStatHF"; 
  if( type == JESRelativePtBBup    || type == JESRelativePtBBdown    ) return "RelativePtBB";   
  if( type == JESRelativePtEC1up   || type == JESRelativePtEC1down   ) return "RelativePtEC1";  
  if( type == JESRelativePtEC2up   || type == JESRelativePtEC2down   ) return "RelativePtEC2";  
  if( type == JESRelativePtHFup    || type == JESRelativePtHFdown    ) return "RelativePtHF";   
  if( type == JESTimeEtaup         || type == JESTimeEtadown         ) return "TimeEta";	 
  if( type == JESAbsoluteScaleup   || type == JESAbsoluteScaledown   ) return "AbsoluteScale";  
  if( type == JESAbsoluteMPFBiasup || type == JESAbsoluteMPFBiasdown ) return "AbsoluteMPFBias";
  if( type == JESAbsoluteStatup    || type == JESAbsoluteStatdown    ) return "AbsoluteStat";   
  if( type == JESSinglePionECALup  || type == JESSinglePionECALdown  ) return "SinglePionECAL"; 
  if( type == JESSinglePionHCALup  || type == JESSinglePionHCALdown  ) return "SinglePionHCAL"; 
  if( type == JESFragmentationup   || type == JESFragmentationdown   ) return "Fragmentation";  
  if( type == JESTimePtup          || type == JESTimePtdown          ) return "TimePt";	    	 
  if( type == JESFlavorQCDup       || type == JESFlavorQCDdown       ) return "FlavorQCD";      

  else return "";
}




    


    












    







