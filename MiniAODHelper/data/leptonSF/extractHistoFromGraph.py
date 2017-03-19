import ROOT
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gDirectory.cd('PyROOT:/')
from array import array

infilename=sys.argv[1]
outfilename=sys.argv[2]

outfile=ROOT.TFile(outfilename,"RECREATE")
infile=ROOT.TFile(infilename,"read")
keylist=infile.GetListOfKeys()
for key in keylist:
  infile.cd()
  thisname=key.GetName()
  print "doing ", thisname
  g=infile.Get(thisname)

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
  print binsXhigh
  a=array("d",binsXlow+[binsXhigh[-1]])
  print a
  #exit(0)
  h=ROOT.TH1D(thisname,thisname,len(a)-1,a)
  h.Sumw2()
  print h.GetNbinsX()

  for iP in range(nP):
    h.SetBinContent(iP+1,binsY[iP])
    h.SetBinError(iP+1,max(abs(errorYlow[iP]),abs(errorYhigh[iP])))
    print iP+1,binsY[iP]
  outfile.cd()
  h.Write()

outfile.Close()



  
  
  