import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gDirectory.cd('PyROOT:/')
from array import array

infile=ROOT.TFile("muons_HIP_ICHEP.root","read")
g=infile.Get("ratio_eta")

binsX=[]
binsXhigh=[]
binsXlow=[]

binsY=[]
errorYlow=[]
errorYhigh=[]

nP=g.GetN()

for iP in range(nP):
  x=ROOT.Double(0)
  y=ROOT.Double(0)
  
  g.GetPoint(iP,x,y)
  
  xlow=g.GetErrorXlow(iP)
  xhigh=g.GetErrorXhigh(iP)
  ylow=g.GetErrorYlow(iP)
  yhigh=g.GetErrorYhigh(iP)
  
  
  binsX.append(x)
  binsY.append(y)
  binsXhigh.append(x+xhigh)
  binsXlow.append(x-xlow)
  errorYhigh.append(yhigh)
  errorYlow.append(ylow)
  print iP, xlow, x, xhigh, ylow, y, yhigh

l=binsXlow
print l
a=array("d",binsXlow+[binsXhigh[-1]])
print a
#exit(0)
h=ROOT.TH1D("ratio_eta","ratio_eta",len(a)-1,a)
h.Sumw2()
print h.GetNbinsX()

for iP in range(nP):
  h.SetBinContent(iP+1,binsY[iP])
  h.SetBinError(iP+1,max(abs(errorYlow[iP]),abs(errorYhigh[iP])))
  print iP+1,binsY[iP]

h.SaveAs("muon_HIP_ICHEP_HISTO.root")



  
  
  