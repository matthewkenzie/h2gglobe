import ROOT as r

import sys

tf = r.TFile(sys.argv[1])
tf_truth = r.TFile(sys.argv[2])

r.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
r.gSystem.Load("../../../libLoopAll.so")

lumi = 19600
normalizer = r.Normalization_8TeV()
normalizer.Init(8)

cuts = [50,75,100,125,150,200,-1]

smTrees = [tf.Get('spin_trees/ggh_m125_8TeV'),tf.Get('spin_trees/vbf_m125_8TeV')] #,tf.Get('spin_trees/wzh_m125_8TeV'),tf.Get('spin_trees/tth_m125_8TeV')]
gravTrees = [tf.Get('spin_trees/gg_grav_m125_8TeV')]
truthTree = tf_truth.Get('CosThetaTree')

r.TH1.SetDefaultSumw2()

smEvs = {}
gravEvs = {}

for cut in cuts:
	smEvs[cut] = r.TH1F('sm%d'%cut,'',20,0,1)
	gravEvs[cut] = r.TH1F('grav%d'%cut,'',20,0,1)
smTruth = r.TH1F('smTruth','',20,0,1)
gravTruth = r.TH1F('gravTruth','',20,0,1)

for tree in smTrees:
	for event in range(tree.GetEntries()):
		tree.GetEntry(event)
		for cut in cuts:
			if cut<0.:
				smEvs[cut].Fill(r.TMath.Abs(tree.costheta_cs),tree.evweight)
			else:
				higgs_pt = r.TLorentzVector(tree.higgs_px,tree.higgs_py,tree.higgs_pz,tree.higgs_E).Pt()
				if higgs_pt<=cut:
					smEvs[cut].Fill(r.TMath.Abs(tree.costheta_cs),tree.evweight)


for tree in gravTrees:
	for event in range(tree.GetEntries()):
		tree.GetEntry(event)
		for cut in cuts:
			if cut<0.:
				gravEvs[cut].Fill(r.TMath.Abs(tree.costheta_cs),tree.evweight)
			else:
				higgs_pt = r.TLorentzVector(tree.higgs_px,tree.higgs_py,tree.higgs_pz,tree.higgs_E).Pt()
				if higgs_pt<=cut:
					gravEvs[cut].Fill(r.TMath.Abs(tree.costheta_cs),tree.evweight)

truthTree.Draw("cosThetaCS>>smTruth","type==2 || type==4","goff")
truthTree.Draw("cosThetaCS>>gravTruth","type==6","goff")

r.gStyle.SetOptStat(0)

smTruth.Scale(1./smTruth.Integral())
gravTruth.Scale(1./gravTruth.Integral())
	
r.gROOT.SetBatch()
canv = r.TCanvas()

colors = [r.kMagenta+2,r.kMagenta-4,r.kBlue,r.kBlue-7,r.kGreen+1,r.kOrange-3,r.kRed]

th2dummy = r.TH2F('h','',1,0.,1.,1,0.35,1.5)
th2dummy.GetXaxis().SetTitle("|cos(#theta*)|")
th2dummy.GetYaxis().SetTitle("eff*acc ratio to SM")
canv.cd()
th2dummy.Draw()
l = r.TLine()
l.SetLineWidth(2)
l.SetLineStyle(r.kDashed)
l.DrawLine(0.,1.,1.,1.)

leg = r.TLegend(0.15,0.65,0.35,0.99)
leg2 = r.TLegend(0.35,0.65,0.55,0.99)
#leg.SetHeader('Cut (GeV)')
leg.SetFillColor(0)
leg.SetLineColor(0)
leg2.SetFillColor(0)
leg2.SetLineColor(0)
leg.SetTextSize(0.04)
leg2.SetTextSize(0.04)

smEffAcc={}
gravEffAcc={}
ratio={}

for i, cut in enumerate(reversed(cuts)):
	# normalize
	smEvs[cut].Scale(1./smEvs[cut].Integral())
	gravEvs[cut].Scale(1./gravEvs[cut].Integral())

	"""
	smEvs.SetLineColor(r.kRed)
	smEvs.SetMarkerColor(r.kRed)
	smEvs.SetLineWidth(2)
	gravEvs.SetLineColor(r.kBlue)
	gravEvs.SetMarkerColor(r.kBlue)
	gravEvs.SetLineWidth(2)
	smTruth.SetLineColor(r.kRed)
	smTruth.SetMarkerColor(r.kRed)
	smTruth.SetLineWidth(2)
	gravTruth.SetLineColor(r.kBlue)
	gravTruth.SetMarkerColor(r.kBlue)
	gravTruth.SetLineWidth(2)


	smEvs.Draw("LEP")
	gravEvs.Draw("LEPsame")
	canv.Print('ctdist.pdf')

	gravTruth.Draw("LEP")
	smTruth.Draw("LEPsame")
	canv.Print('truthdist.pdf')
	"""

	smEffAcc[cut] = smEvs[cut].Clone('%s_ea'%smEvs[cut].GetName())
	gravEffAcc[cut] = gravEvs[cut].Clone('%s_ea'%gravEvs[cut].GetName())

	smEffAcc[cut].Divide(smTruth)
	gravEffAcc[cut].Divide(gravTruth)

	gravEffAcc[cut].Scale(1./gravEffAcc[cut].Integral())
	smEffAcc[cut].Scale(1./smEffAcc[cut].Integral())

	#gravEffAcc[cut].Draw("LEP")
	#smEffAcc[cut].Draw("LEPsame")
	#canv.Print('eadist.pdf')

	ratio[cut] = gravEffAcc[cut].Clone('ratio%d'%cut)
	ratio[cut].Divide(smEffAcc[cut])
	ratio[cut].SetLineColor(colors[i])
	ratio[cut].SetLineWidth(3)
	name = '<%d'%cut
	if cut<0.: name='None'
	if i<4:
		leg.AddEntry(ratio[cut],name,"LEP")
	else:
		leg2.AddEntry(ratio[cut],name,"LEP")
	ratio[cut].Draw("LEPsame")

leg.Draw("same")
leg2.Draw("same")
canv.Print('ratio.pdf')

