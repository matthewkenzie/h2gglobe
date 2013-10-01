#!/usr/bin/env python

import sys

logfile = open(sys.argv[1])

import array
import ROOT as r

hcut=0.
jcut=0.

atRelevantPart=False

cuts = [50,75,100,125,150,200,-1]

hcutHist = r.TH1F('hcut','',len(cuts),0,len(cuts))
jcutHist = r.TH1F('jcut','',len(cuts),0,len(cuts))
hcutHist.GetXaxis().SetLabelSize(0.045)
jcutHist.GetXaxis().SetLabelSize(0.045)
hcutHist.GetXaxis().SetTitleOffset(1.1)
jcutHist.GetXaxis().SetTitleOffset(1.1)
hcutHist.GetYaxis().SetTitleOffset(1.1)
jcutHist.GetYaxis().SetTitleOffset(1.1)

for p,cut in enumerate(cuts):
	name = '<%d'%cut
	if cut<0: name='None'
	hcutHist.GetXaxis().SetBinLabel(p+1,name)
	jcutHist.GetXaxis().SetBinLabel(p+1,name)

for line in logfile.readlines():
	if line.startswith('[\'qqbarComb'):
		els = line.split()
		el = els[1].split('/qmu')[0]
		hcut = int(el.split('_jcut')[0].split('hcut')[1])
		jcut = int(el.split('_jcut')[1])
		print hcut, jcut
	if line.startswith('Plotting separation for fqq =  0.0'):
		atRelevantPart=True
	if line.startswith('Separation from histograms'):
		atRelevantPart=False
	if line.startswith('Median point prob ALT:') and atRelevantPart:
		els = line.split()
		#print els
		prob = float(els[4])
		sigma = r.RooStats.PValueToSignificance(prob)
		expcls = 1.-prob/0.5
		print prob, sigma, expcls
		if hcut<0.:
			jcutHist.SetBinContent(cuts.index(jcut)+1,prob)
			jcutHist.SetBinError(cuts.index(jcut)+1,r.TMath.Sqrt((1.-prob)*prob*5000.)/5000.)
		if jcut<0.:
			hcutHist.SetBinContent(cuts.index(hcut)+1,prob)
			hcutHist.SetBinError(cuts.index(hcut)+1,r.TMath.Sqrt((1.-prob)*prob*5000.)/5000.)

jcutHist.SetLineWidth(3)
jcutHist.SetLineColor(r.kRed)
jcutHist.GetXaxis().SetTitle('Cut value in p_{T}')
jcutHist.GetYaxis().SetTitle('Tail probability')
jcutHist.SetMinimum(0.)

hcutHist.SetLineWidth(3)
hcutHist.SetLineColor(r.kBlue)
hcutHist.GetXaxis().SetTitle('Cut value in p_{T}')
hcutHist.GetYaxis().SetTitle('Tail probability')
hcutHist.SetMinimum(0.)

r.gStyle.SetOptStat(0)
r.gROOT.SetBatch()

leg = r.TLegend(0.2,0.2,0.6,0.4)
leg.SetFillColor(0)
leg.AddEntry(hcutHist,"Higgs p_{T} cut","LEP")
leg.AddEntry(jcutHist,"Lead jet p_{T} cut","LEP")

canv = r.TCanvas()

sigmaLine = r.TLine()
sigmaLine.SetLineWidth(2)
sigmaLine.SetLineStyle(r.kDashed)

jcutHist.Draw('LEP')
for sigma in [0.5,0.75,1.,1.25,1.50,2.0]:
	sigmaLine.DrawLine(0.,r.RooStats.SignificanceToPValue(sigma),len(cuts),r.RooStats.SignificanceToPValue(sigma))
	text = r.TLatex()
	text.DrawLatex(len(cuts)+0.1,r.RooStats.SignificanceToPValue(sigma),'%4.2f#sigma'%sigma)
jcutHist.Draw('LEPsame')
hcutHist.Draw('LEPsame')
leg.Draw("same")
canv.Print('ptCutSep.pdf')


	

