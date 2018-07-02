# Data PU profile 2017
Nominal profile for min-bias xs of 69.2 mb with up/down variation of +/-4.6%:
- DataPileupHistogram_Run2017_294927-306462_13TeV_EOY2017ReReco_MinBiasDown-66017.root
- DataPileupHistogram_Run2017_294927-306462_13TeV_EOY2017ReReco_MinBiasNominal-69200.root
- DataPileupHistogram_Run2017_294927-306462_13TeV_EOY2017ReReco_MinBiasUp-72383.root
Profiles compared in `DataPileupHistogram_Run2017_294927-306462_13TeV_EOY2017ReReco.pdf`.

Produced following [PileupJSONFileforData TWiki r29](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II):
```
pileupCalc.py -i Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt  --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 DataPileupHistogram_Run2017_294927-306462_13TeV_EOY2017ReReco_MinBiasNominal-69200.root
```
with `pileup_latest.txt` norm tag from 2018-07-02.