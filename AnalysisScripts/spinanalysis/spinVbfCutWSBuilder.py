import os
import sys

if len(sys.argv)<5:
  print 'usage: python spinVbfCutWSBuilder.py <tree_file> <ws_file> <higgs_pt_cut> <lead_j_pt_cut>'
  sys.exit()

higgs_pt_cut = float(sys.argv[3])
lead_j_pt_cut = float(sys.argv[4])

import ROOT as r
r.gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so')
r.gSystem.Load('../../libLoopAll.so')

tf = r.TFile(sys.argv[1])
wstf = r.TFile(sys.argv[2])
outf = r.TFile('CMS-HGG_hcut%.0f_jcut%.0f.root'%(higgs_pt_cut,lead_j_pt_cut),'RECREATE')
outws = r.RooWorkspace('cms_hgg_workspace','cms_hgg_workspace')
bkgws = wstf.Get('cms_hgg_workspace')

treeToDset = {'data_mass': tf.Get('spin_trees/Data')}
treeToDset['sig_ggh_mass_m125'] = tf.Get('spin_trees/ggh_m125_8TeV') 
treeToDset['sig_vbf_mass_m125'] = tf.Get('spin_trees/vbf_m125_8TeV')
treeToDset['sig_wzh_mass_m125'] = tf.Get('spin_trees/wzh_m125_8TeV')
treeToDset['sig_tth_mass_m125'] = tf.Get('spin_trees/tth_m125_8TeV')
treeToDset['sig_gg_grav_mass_m125'] = tf.Get('spin_trees/gg_grav_m125_8TeV')
treeToDset['sig_qq_grav_mass_m125'] = tf.Get('spin_trees/qq_grav_m125_8TeV')

mgg = r.RooRealVar('CMS_hgg_mass','CMS_hgg_mass',100,180)
mgg.setBins(160)
wt = r.RooRealVar('weight','weight',0.,100.)

for name,tree in treeToDset.items():
	dataSet=[]
	th1f=[]
	for c in range(20):
		dataSet.append(r.RooDataSet('%s_cat%d'%(name,c),'',r.RooArgSet(mgg,wt),'weight'))
		th1f.append(r.TH1F('th1f_%s'%dataSet[c].GetName(),'th1f_%s'%dataSet[c].GetName(),160,100,180))
	for ev in range(tree.GetEntries()):
		if ev%1000==0: print ev, '/', tree.GetEntries()
		tree.GetEntry(ev)
		cat = tree.category
		weight = tree.evweight
		#if 'gg_grav' in name:
		#	weight *= 1.064
		#if 'qq_grav' in name:
		#	weight *= 1.19
		higgs_p4 = r.TLorentzVector(tree.higgs_px,tree.higgs_py,tree.higgs_pz,tree.higgs_E)
		higgs_pt = higgs_p4.Pt()
		higgs_mass = higgs_p4.M()
		lead_j_pt = tree.myVBFLeadJPt
		if (higgs_pt_cut<0. or higgs_pt<higgs_pt_cut) and (lead_j_pt_cut<0. or lead_j_pt<lead_j_pt_cut):
			mgg.setVal(higgs_mass)
			dataSet[cat].add(r.RooArgSet(mgg),weight)
			th1f[cat].Fill(higgs_mass,weight)
	for c in range(20):
		if 'data' in name:
			bkgMod = bkgws.pdf('data_pol_model_8TeV_cat%d'%c)
			bkgMod.fitTo(dataSet[c])
			getattr(outws,'import')(bkgMod)
		getattr(outws,'import')(dataSet[c])
		dataBinned = dataSet[c].binnedClone('roohist_%s'%dataSet[c].GetName())
		getattr(outws,'import')(dataBinned)
		th1f[c].Write()

outf.cd()
outws.Write()
outf.Close()
tf.Close()
