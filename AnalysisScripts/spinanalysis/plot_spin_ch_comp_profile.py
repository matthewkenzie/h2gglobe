#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--dir")
parser.add_option("-s","--spinCats",type="int",default=5)
parser.add_option("-y","--yaxis",default="-1.5,3.5")
parser.add_option("-u","--unblind",default=False,action="store_true")
parser.add_option("-b","--batch",default=False,action="store_true")
parser.add_option("-o","--outfile",default="chcomp")
parser.add_option("-S","--sqrtS",default="comb")
parser.add_option("-L","--legPos",default="TL")
parser.add_option("-H","--drawAsHist",default=False,action="store_true")
(options,args) = parser.parse_args()

import ROOT as r
import sys,os,array
r.gROOT.ProcessLine(".x hggPaperStyle.C")

if options.legPos!="TL" and options.legPos!="TR" and options.legPos!="BL" and options.legPos!="BR":
	sys.exit("Invalid legend position"+options.legPos)

def findQuantile(pts,cl):

	#gr is a list of r,nll
	# start by walking along the variable and check if crosses a CL point
	crossbound = [ pt[1]<=cl for pt in pts ]
	rcrossbound = crossbound[:]
	rcrossbound.reverse()

	minci = 0
	maxci = len(crossbound)-1
	min = pts[0][0]
	max = pts[maxci][0]

	for c_i,c in enumerate(crossbound): 
		if c : 
			minci=c_i
			break
	
	for c_i,c in enumerate(rcrossbound): 
		if c : 
			maxci=len(rcrossbound)-c_i-1
			break

	if minci>0: 
		y0,x0 = pts[minci-1][0],pts[minci-1][1]
		y1,x1 = pts[minci][0],pts[minci][1]
		min = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)
		
	if maxci<len(crossbound)-1: 
		y0,x0 = pts[maxci][0],pts[maxci][1]
		y1,x1 = pts[maxci+1][0],pts[maxci+1][1]
		max = y0+((cl-x0)*y1 - (cl-x0)*y0)/(x1-x0)

	return min,max

def plot1DNLL(file,poiname,name=""):

	x = poiname

	clean_graphs=[]

	tf = r.TFile(file)
	tree = tf.Get('limit')
	gr = r.TGraph()
	gr.SetName('gr_%s_%s'%(x,name))

	res=[]
	for i in range(tree.GetEntries()):
		tree.GetEntry(i)
		xv = getattr(tree,x)
		#if tree.quantileExpected==1: continue
		if tree.deltaNLL<0: continue
		if xv in [re[0] for re in res]: continue
		if tree.deltaNLL>100 or tree.deltaNLL<-100: continue
		if tree.deltaNLL != tree.deltaNLL: continue
		if xv<-1 and tree.deltaNLL==0: continue
		res.append([xv,2.*tree.deltaNLL])
	res.sort()
	minNLL = min([re[1] for re in res])
	for re in res: 
		#if options.correctNLL and re[1]==0.: re[1]=-1
		re[1]-=minNLL
	
	p=0
	for re, nll in res: 
		if nll>=0.:
			gr.SetPoint(p,re,nll)
			p+=1

	m,m1 = findQuantile(res,0);
	#m = 1.
	#m1 = 1.
	l,h  = findQuantile(res,1);
	l2,h2  = findQuantile(res,4);
	
	clean_graphs.append(gr)

	xmin = m
	eplus = h-m
	eminus = m-l
	eplus2 = h2-m
	eminus2 = m-l2

	print "%15s - %4.2f +%4.2g -%4.2g" % ( name, xmin, eplus , eminus )
	return [xmin,eplus,eminus,eplus2,eminus2]

if options.batch:
	r.gROOT.SetBatch()

fnames={}
fnames['Data'] = 'SMFitToDataSpinCat'
#fnames['Data'] = 'None'
#fnames['DataNoErr'] = 'SMFitToDataMultiSignal' 
fnames['SM'] = '_sm_fit_to_sm_asimov_'
fnames['GravGG'] = '_sm_fit_to_gravgg_asimov_' 
fnames['Grav50GG50QQ'] = '_sm_fit_to_grav50gg50qq_asimov_' 
fnames['GravQQ'] = '_sm_fit_to_gravqq_asimov_'

files={}
for key in fnames.keys(): files[key] = []

for root,dir,fs in os.walk(options.dir):
	if os.path.basename(root) != os.path.basename(options.dir): continue
	for f in fs:
		for key, fname in fnames.items():
			if 'higgsCombine%s'%fname in f and 'Job' not in f:
				files[key].append(f)

print files

key_order = ['Data','SM','GravGG','Grav50GG50QQ','GravQQ']
opts={}
opts['Data'] = 						[r.kBlack,		r.kFullCircle,			 1.2, "Observed"]
#opts['DataNoErr'] = 			[r.kBlack,		r.kFullCircle,			 1.2, "Observed"]
opts['SM'] = 							[r.kRed,			r.kFullSquare,			 1.2, "Expected J^{P} = 0^{+} (SM)"]
opts['GravGG'] = 					[r.kBlue,			r.kFullTriangleUp,	 1.6, "Expected J^{P} = 2^{+}_{m} (f_{q#bar{q}}=0.0)"]
opts['Grav50GG50QQ'] = 		[r.kMagenta,	r.kFullDiamond,			 1.9, "Expected J^{P} = 2^{+}_{m} (f_{q#bar{q}}=0.5)"]
opts['GravQQ'] = 					[r.kGreen+1,	r.kFullTriangleDown, 1.6, "Expected J^{P} = 2^{+}_{m} (f_{q#bar{q}}=1.0)"]
dummyHist = r.TH2F('dummy','',1,0.,1.,1,float(options.yaxis.split(',')[0]),float(options.yaxis.split(',')[1]))
dummyHist.GetXaxis().SetTitle('|cos#theta*|')
dummyHist.GetYaxis().SetTitle('#mu')#sigma/#sigma_{SM}')

x=[0.1,0.2875,0.4625,0.65,0.875]
boundaries = array.array('f',[0.0,0.2,0.375,0.55,0.75,1.0])
boundHist = r.TH1F('bound','',len(boundaries)-1,boundaries)
canv = r.TCanvas()
#canv.SetGrid()

dummyHist.Draw("AXISG")
dummyHist.SetStats(0)
if options.legPos=="TL":
	leg = r.TLegend(0.15,0.62,0.6,0.91)
elif options.legPos=="BL":
	leg = r.TLegend(0.11,0.11,0.4,0.35)
elif options.legPos=="TR":
	leg = r.TLegend(0.6,0.65,0.89,0.89)
else: # BR
	leg = r.TLegend(0.6,0.11,0.89,0.35)
leg.SetFillColor(0)
leg.SetLineColor(0)

fkeys = files.keys()
fkeys.sort()
print fkeys

graphs = {}
hists = {}
for key in key_order:
	graphs[key] = r.TGraphAsymmErrors()
	graphs[key].SetName(key)
	graphs[key].SetLineWidth(2)
	graphs[key].SetLineColor(opts[key][0])
	graphs[key].SetMarkerColor(opts[key][0])
	graphs[key].SetMarkerStyle(opts[key][1])
	graphs[key].SetMarkerSize(opts[key][2])
	hists[key] = r.TH1F('hist_%s'%key,'',len(boundaries)-1,boundaries)
	hists[key].SetLineColor(opts[key][0])
	hists[key].SetLineWidth(4)

print graphs
print hists

for key in key_order:
	if key=='Data': 
		leg.AddEntry(graphs[key],opts[key][3],"EP")
		print 'oiy\n', files[key]
		sfiles = files[key].sort()
		print sfiles
		for p,file in enumerate(files[key]):
			print file
			res = plot1DNLL(options.dir+'/'+file,"r_spinCat%d"%p,"spinCat%d"%p)
			graphs[key].SetPoint(p,x[p],res[0])
			graphs[key].SetPointError(p,x[p]-boundHist.GetBinLowEdge(p+1),boundHist.GetBinLowEdge(p+2)-x[p],res[2],res[1])
	else: 
		if key!='DataNoErr': 
			if options.drawAsHist: leg.AddEntry(hists[key],opts[key][3],"L")
			else: leg.AddEntry(graphs[key],opts[key][3],"LP")
		tf = r.TFile(options.dir+'/'+files[key][0])
		tree = tf.Get('limit')
		tree.GetEntry(0)
		for cat in range(options.spinCats):
			graphs[key].SetPoint(cat,x[cat],getattr(tree,'r_spinCat%d'%cat))
			hists[key].SetBinContent(cat+1,getattr(tree,'r_spinCat%d'%cat))
	
	if key=='Data': 
		if options.unblind: 
			graphs[key].Draw("PEsame")
	elif key=='DataNoErr':
		if options.unblind:
			graphs[key].Draw("Psame")
	else: 
		if options.drawAsHist:
			hists[key].Draw("HISTsame")
		else:
			graphs[key].Draw("LPsame")

line = r.TLine(0.,0.,1.,0)
line.SetLineWidth(2)
line.SetLineStyle(r.kDashed)
line.Draw("same")
leg.Draw("same")
if options.unblind: 
	if 'DataNoErr' in graphs.keys(): graphs['DataNoErr'].Draw("Psame")
	if 'Data' in graphs.keys(): graphs['Data'].Draw("PEsame")

lat = r.TLatex()
lat.SetNDC()

if options.sqrtS=="comb":
	text = '#sqrt{s}=7TeV L=5.1fb^{-1}, #sqrt{s}=8TeV L=19.7fb^{-1}'
elif options.sqrtS=="7" or options.sqrtS=="7TeV":
	text = '#sqrt{s}=7TeV L=5.1fb^{-1}'
elif options.sqrtS=="8" or options.sqrtS=="8TeV":
	text = '#sqrt{s}=8TeV L=19.7fb^{-1}'
else:
	sys.exit("Invalid sqrtS"+options.sqrtS)

lat.DrawLatex(0.129,0.93,"CMS H#rightarrow#gamma#gamma")
lat.DrawLatex(0.438,0.93,text)
canv.Update()
canv.Modified()
if not options.batch:
	raw_input('Good?\n')
canv.Print('%s.pdf'%options.outfile)
canv.Print('%s.png'%options.outfile)
canv.Print('%s.root'%options.outfile)

