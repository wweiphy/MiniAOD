#!/cvmfs/cms.cern.ch/slc6_amd64_gcc481/cms/cmssw/CMSSW_7_2_3/external/slc6_amd64_gcc481/bin/python

import ROOT

fileDir = "/nfs/dust/cms/user/shwillia/BoostedTrees/Moriond_20160715/"

trueFileDir = fileDir+"ttbar_incl_nominal_Tree.root"
falseFileDir = fileDir+"ttbar_incl_nominal_Tree.root"

trueDef = "(BoostedJet_Top_Pt>0)&&((BoostedJet_Dr_GenB>=0&&BoostedJet_Dr_GenB<1.5)&&(BoostedJet_Dr_GenQ1>=0&&BoostedJet_Dr_GenQ1<1.5)&&(BoostedJet_Dr_GenQ2>=0&&BoostedJet_Dr_GenQ2<1.5))"
falseDef = "(BoostedJet_Top_Pt>0)&&(((BoostedJet_Dr_GenB>=0&&BoostedJet_Dr_GenB>1.5)||(BoostedJet_Dr_GenQ1>=0&&BoostedJet_Dr_GenQ1>1.5)||(BoostedJet_Dr_GenQ2>=0&&BoostedJet_Dr_GenQ2>1.5))&&(BoostedJet_Dr_GenTopHad>=0&&BoostedJet_Dr_GenTopHad>2.0))"

trueFile = ROOT.TFile(trueFileDir)
falseFile = ROOT.TFile(falseFileDir)

trueTree = trueFile.Get("MVATree")
falseTree = falseFile.Get("MVATree")

variableNames = [ ("BoostedJet_Top_M","BoostedJet_Top_M",100,0.,300.),
                  ("BoostedJet_MRatio_23_Top","BoostedJet_M23/BoostedJet_Top_M",80,0.,1.),
                  ("BoostedJet_Atan_1213","atan(BoostedJet_M12/BoostedJet_M13)",80,0.,2.),
                  ("BoostedJet_MRatio_W_Top","BoostedJet_MRatio_W_Top",80,0.,1.),
                  ("BoostedJet_Atan_BW1W2","atan(BoostedJet_BW1_M/BoostedJet_BW2_M)",80,0.,2.),
                  ("BoostedJet_MRatio_Wbtag_Top","BoostedJet_MRatio_Wbtag_Top",80,0.,1.),
                  ("BoostedJet_Atan_BW1W2btag","atan(BoostedJet_BW1btag_M/BoostedJet_BW2btag_M)",80,0.,2.)
                ]

outFile = ROOT.TFile("toplikelihoodtaggerhistos.root","RECREATE")

for name,var,bins,xmin,xmax in variableNames:
  
  print name
  
  trueHist = ROOT.TH1F(name+"_True",name+"_True",bins,xmin,xmax)
  falseHist = ROOT.TH1F(name+"_False",name+"_False",bins,xmin,xmax)
  
  trueTree.Draw(var+">>"+name+"_True",trueDef)
  falseTree.Draw(var+">>"+name+"_False",falseDef)

  trueHist.Scale(1/trueHist.Integral())
  falseHist.Scale(1/falseHist.Integral())
  
  trueHist.Write()
  falseHist.Write()

outFile.Close()

trueFile.Close()
falseFile.Close()
