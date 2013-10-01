import sys

import ROOT as r

r.gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so')
r.gSystem.Load('../../libLoopAll.so')
r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

ncats=20
samples = ['data','ggh','vbf','wzh','tth','gg_grav','qq_grav']
cuts = [50,75,100,125,150,200,-1]
th1fs_hcut = {}
th1fs_jcut = {}
files_hcut = {}
files_jcut = {}

for cut in cuts:
	files_hcut[cut] = r.TFile('CMS-HGG_hcut%d_jcut-1.root'%cut)
	files_jcut[cut] = r.TFile('CMS-HGG_hcut-1_jcut%d.root'%cut)
	name = str(cut)
	if cut<0.: name = 'None'
	th1fs_hcut[cut] = r.TH1F('hcut_%s'%(name),'',7,0,7)
	th1fs_jcut[cut] = r.TH1F('jcut_%s'%(name),'',7,0,7)


def makeCutFlowStylePlot(cuts,files,th1fs,name,title):

	loosestVal = {}
	for samp in samples: loosestVal[samp]=0.

	for cut in cuts:
		ws = files[cut].Get('cms_hgg_workspace')
		for j,samp in enumerate(samples):
			dsTot=0.
			for cat in range(ncats):
				if samp=='data':
					ds = ws.data('%s_mass_cat%d'%(samp,cat))
				else:
					ds = ws.data('sig_%s_mass_m125_cat%d'%(samp,cat))
				dsTot += ds.sumEntries()
				if dsTot>loosestVal[samp]: loosestVal[samp]=dsTot
			th1fs[cut].SetBinContent(j+1,dsTot)
			th1fs[cut].GetXaxis().SetBinLabel(j+1,samp)

	for cut in cuts:
		th1fs[cut].SetBarWidth(0.9)
		th1fs[cut].SetBarOffset(0.05)
		for j, samp in enumerate(samples):
			bC = th1fs[cut].GetBinContent(j+1)
			th1fs[cut].SetBinContent(j+1,bC/loosestVal[samp])

	canv = r.TCanvas()
	colors = [r.kMagenta+2,r.kMagenta-4,r.kBlue,r.kBlue-7,r.kGreen+1,r.kOrange-3,r.kRed]
	th1fs[-1].SetTitle(title)
	th1fs[-1].GetYaxis().SetRangeUser(0.,1.1)
	th1fs[-1].GetYaxis().SetTitle('Fraction of events lost')
	th1fs[-1].GetXaxis().SetLabelSize(0.05)
	leg = r.TLegend(0.87,0.5,0.99,0.99)
	leg.SetHeader('Cut (GeV)')
	leg.SetFillColor(0)
	leg.SetTextSize(0.035)
	for i,cut in enumerate(reversed(cuts)):
		th1fs[cut].SetLineColor(colors[i])
		th1fs[cut].SetFillColor(colors[i])
		legname = '<'+str(cut)
		if cut<0.: legname='None'
		leg.AddEntry(th1fs[cut],legname,'F')
		if i==0:
			th1fs[cut].Draw("BAR")
		else:
			th1fs[cut].Draw("BARsame")

	leg.Draw("same")
	#raw_input('?')
	canv.Print('%s.pdf'%name)

# __main__

makeCutFlowStylePlot(cuts,files_hcut,th1fs_hcut,'hcut','Higgs p_{T} cut')
makeCutFlowStylePlot(cuts,files_jcut,th1fs_jcut,'jcut','Lead jet p_{T} cut')


"""
for file in files:
	
	tf = r.TFile(file)
	ws = tf.Get('cms_hgg_workspace')

	dataTot=0.
	for cat in range(20):
		ds = ws.data('data_mass_cat%d'%cat)
		#print ds.GetName(), '--', ds.sumEntries()
		dataTot+=ds.sumEntries()
	print tf.GetName(), '--', dataTot

	for cat in range(20):
		ds = ws.data('sig_ggh_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()

	for cat in range(20):
		ds = ws.data('sig_vbf_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()

	for cat in range(20):
		ds = ws.data('sig_wzh_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()

	for cat in range(20):
		ds = ws.data('sig_tth_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()

	for cat in range(20):
		ds = ws.data('sig_gg_grav_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()

	for cat in range(20):
		ds = ws.data('sig_qq_grav_mass_m125_cat%d'%cat)
		print ds.GetName(), '--', ds.sumEntries()
	"""



