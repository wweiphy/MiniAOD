# Data PU profile 2018
Nominal profile for min-bias xs of 69.2 mb with up/down variation of +/-4.6%:
  - "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasNominal-69200.root";  
  - "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasUp-72383.root";
  - "DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasDown-66017.root";
Profiles compared in `DataPileupHistogram_Run2018_294927-306462_13TeV_PromptReco_Collisions18.pdf`.

Produced following [PileupJSONFileforData TWiki r29](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II):
```
pileupCalc.py -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_MinBiasNominal-69200.root
```
with `pileup_latest.txt` norm tag from 2019-04-01.