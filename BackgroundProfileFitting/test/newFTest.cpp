#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TIterator.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLatex.h"

#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooDataHist.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "../interface/PdfModelBuilder.h"

using namespace std;
using namespace RooFit;
using namespace boost;
namespace po = program_options;

string infilename_;
string outfilename_;
string outdatfilename_;
string sigfilename_;
bool doSBversion_=false;
bool is2011_=false;
int sqrtS_=8;
int singleCat_;
int startingCat_=0;
int ncats_;
bool blind_=true;
bool verbose_=false;

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext=""){
  
  if (type=="Bernstein") return pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order); 
  else if (type=="Chebychev") return pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order); 
  else if (type=="Exponential") return pdfsModel.getExponentialSingle(Form("%s_exp%d",ext,order),order,false); 
  else if (type=="PowerLaw") return pdfsModel.getPowerLawSingle(Form("%s_pow%d",ext,order),order,false); 
  else if (type=="Laurent") return pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order); 
  else {
    cerr << "ERROR -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

RooAbsPdf* getSBPdf(RooAbsPdf *bkgPdf, RooWorkspace *sigWS, int cat, RooRealVar *MH, RooRealVar *mu, double dataEvs) {
	
	RooRealVar *bkgNEvs = new RooRealVar("bkgNEvs","",dataEvs,dataEvs-5.*TMath::Sqrt(dataEvs),dataEvs+5.*TMath::Sqrt(dataEvs));
	RooAbsPdf *sigPdf = (RooAbsPdf*)sigWS->pdf(Form("sigpdfrelcat%d_allProcs",cat));
	MH->setConstant(false);
	MH->setVal(125);
	MH->setRange(124,126);
	// have to set nuisances to constant
	RooArgSet vars = sigWS->allVars();
	TIterator *iter = vars.createIterator();
	RooAbsArg *parg;
	while ((parg=(RooAbsArg*)iter->Next())) {
		if (TString(parg->GetName()).Contains("nuisance") or TString(parg->GetName()).Contains("IntLumi")) {
			RooRealVar *var = (RooRealVar*)sigWS->var(parg->GetName());
			var->setConstant(true);
		}
	}
	// back to where we were
	RooAddPdf *sbPdf = new RooAddPdf(Form("sb%s",bkgPdf->GetName()),"",RooArgList(*bkgPdf,*sigPdf),RooArgList(*bkgNEvs,*mu));
	return sbPdf;
}

void calcChi2WithPlot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, int binning=80) {

	TCanvas *canv = new TCanvas();
	RooPlot *plot = mass->frame();
	data->plotOn(plot,Binning(binning),Name("data"));
	pdf->plotOn(plot,Name("model"));
	plot->Draw();
	canv->Print("temp.pdf");
	double chi2_rf = plot->chiSquare("model","data",4);

	RooCurve *curve = (RooCurve*)plot->findObject("model");
	RooHist *hist = (RooHist*)plot->findObject("data");

	hist->Draw("EP");
	curve->Draw("Lsame");
	canv->Print("temp2.pdf");

	double chi2_me = 0.;
	for (int i=0; i<hist->GetN(); i++){
		double x,y,err,fy;
		hist->GetPoint(i,x,y);
		fy = curve->Eval(x);
		if (fy >= y) {
			err = hist->GetErrorYhigh(i);
		}
		else {
			err = hist->GetErrorYlow(i);
		}
		//cout << i << " " << x << " " << y << " " << fy << endl;
		//cout << "\t" << err << endl;
		double pull = (y-fy)/err;
		chi2_me += pull*pull;
	}
	
	int nbins = hist->GetN();
	int ndof = pdf->getVariables()->getSize()-1;

	double chi2_me_pdof = chi2_me/(nbins-ndof);
	double prob_rf = TMath::Prob(chi2_rf*(nbins-ndof),nbins-ndof);
	double prob_me = TMath::Prob(chi2_me,nbins-ndof);
	cout << nbins << " " << ndof << " " << nbins-ndof << endl;
	cout << chi2_rf << " " << chi2_me << " " << chi2_me_pdof << " " << prob_rf << " " << prob_me << endl;

}

double calcChi2Full(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, int ndof, int binning=80, bool perdof=false, bool chi2prob=false){
	
	RooPlot *plot = mass->frame();
	data->plotOn(plot,Binning(binning),Name("data"));
	pdf->plotOn(plot,Name("model"));
	plot->Draw();
	
	RooCurve *curve = (RooCurve*)plot->findObject("model");
	RooHist *hist = (RooHist*)plot->findObject("data");
	
	double chi2_me = 0.;
	for (int i=0; i<hist->GetN(); i++){
		double x,y,err,fy;
		hist->GetPoint(i,x,y);
		fy = curve->Eval(x);
		if (fy >= y) {
			err = hist->GetErrorYhigh(i);
		}
		else {
			err = hist->GetErrorYlow(i);
		}
		//cout << i << " " << x << " " << y << " " << fy << endl;
		//cout << "\t" << err << endl;
		double pull = (y-fy)/err;
		chi2_me += pull*pull;
	}
	int nbins = hist->GetN();
	if (perdof) return chi2_me/(nbins-ndof);
	double prob_me = TMath::Prob(chi2_me,nbins-ndof);
	if (chi2prob) return prob_me; 
	return chi2_me;
}

double calcChi2(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, int ndof, int binning=80){
	return calcChi2Full(mass,pdf,data,ndof,binning,false,false);
}

double calcChi2pdof(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, int ndof, int binning=80){
	return calcChi2Full(mass,pdf,data,ndof,binning,true,false);
}

double calcChi2prob(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, int ndof, int binning=80){
	return calcChi2Full(mass,pdf,data,ndof,binning,false,true);
}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, RooArgSet *params, string name, double chi2, double gofProb, RooAbsPdf *bkgPdf=NULL){

	TCanvas *canv = new TCanvas();
	RooPlot *plot = mass->frame();
	plot->SetTitle("");
	plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	if (blind_) {
		data->plotOn(plot,Binning(80),Invisible());
		mass->setRange("cutRangeLow",100,110);
		mass->setRange("cutRangeHigh",150,180);
		data->plotOn(plot,Binning(80),CutRange("cutRangeLow,cutRangeHigh"));
	}
	else {
		data->plotOn(plot,Binning(80));
	}
	if (doSBversion_) {
		pdf->plotOn(plot,Components(*bkgPdf));
	}
	else {
		pdf->plotOn(plot);
	}
	pdf->paramOn(plot,Parameters(*params));
	plot->Draw();
  
	TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("#chi^{2} = %.3f, Prob = %.2f ",chi2,gofProb));
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  delete canv;
  delete lat;
}

void plot(RooRealVar *mass, RooMultiPdf *mpdf, RooAbsData *data, RooCategory *category, int bestFitIndex, int cat, string name){
	
  int colors[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
	TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
	RooPlot *plot = mass->frame();
	plot->SetTitle(Form("Category %d",cat));
	plot->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
	
	data->plotOn(plot,Binning(80));
	TObject *datLeg = plot->getObject(int(plot->numItems()-1));
	leg->AddEntry(datLeg,Form("Data - cat%d",cat),"LEP");
	int colorInd=0;
	int style=1;
	for (int ind=0; ind<category->numTypes(); ind++){
		category->setIndex(ind);
		if (colorInd<6) colorInd++;
		else {
			colorInd=0;
			style++;
		}
		mpdf->getCurrentPdf()->plotOn(plot,LineColor(colors[colorInd]),LineStyle(style));
		TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
		string ext="";
		if (bestFitIndex==ind) ext= " (Best Fit Pdf)";
		leg->AddEntry(pdfLeg,Form("%s%s",mpdf->getCurrentPdf()->GetTitle(),ext.c_str()),"L");
	}
	plot->Draw();
	leg->Draw("same");
	canv->SaveAs(Form("%s.pdf",name.c_str()));
	delete canv;
}

int getBestFitFunction(RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, bool silent=false){

	RooRealVar nBackground("bkgshape_norm","nbkg",data->sumEntries(),0,10E8);
	RooExtendPdf extPdf("internal","testextend",*bkg,nBackground);

	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters(*data);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		std::cout << "CLEAN SET OF PARAMETERS" << std::endl;
		params->Print("V");
		std::cout << "-----------------------" << std::endl;
	}
	
	//bkg->setDirtyInhibit(1);
	RooAbsReal *nllm = extPdf.createNLL(*data,RooFit::Extended());
	RooMinimizer minim(*nllm);
	minim.setStrategy(1);
	
	for (int id=0;id<number_of_indeces;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		if (!silent) {
			std::cout << "BEFORE FITTING" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;
		}
		
		minim.minimize("Minuit2","minimize");
		double minNll = (nllm->getVal())+bkg->getCorrection();
		if (!silent) {
			std::cout << "After Minimization ------------------  " <<std::endl;
			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			std::cout << " ------------------------------------  " << std::endl;
	
			std::cout << "AFTER FITTING" << std::endl;
			params->Print("V");
			std::cout << " Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << " Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << " NLL + c = " <<std::setprecision(10) <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}
	params->assignValueOnly(snap);
    	cat->setIndex(best_index);
	
	if (!silent) {
		std::cout << "Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}

void OptionParser(int argc, char* argv[]){

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
		("infilename,i", po::value<string>(&infilename_),																						"In file with data")
		("outfilename,o", po::value<string>(&outfilename_)->default_value("mpdfws.root"),						"Save multipdf workspace here")
		("outdatfilename,d", po::value<string>(&outdatfilename_)->default_value("dat/fTest.dat"),		"Save bias study dat file here")
		("sigfilename,s", po::value<string>(&sigfilename_),																					"Add a file signal file for SB fits")
		("is2011",																																									"Is 2011 run")
		("singleCat", po::value<int>(&singleCat_)->default_value(-1),																"Run single category")
		("ncats,n", po::value<int>(&ncats_)->default_value(14),																			"Number of categories")
		("unblind", 																																								"Unblind")
    ("verbose,v",                                                                               "Run with more output")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);
	if (vm.count("help")) { cout << desc << endl; exit(1); }
	if (vm.count("is2011")) { is2011_=true; sqrtS_=7; }
	if (vm.count("sigfilename")) doSBversion_=true;
	if (vm.count("unblind")) blind_=false;
	if (vm.count("verbose")) verbose_=true;
	if (singleCat_<0) {
		startingCat_=0;
	}
	else {
		startingCat_ = singleCat_;
		ncats_ = singleCat_+1;
	}
  if (!verbose_) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
  }
}

int main(int argc, char* argv[]){

	OptionParser(argc,argv);

	gROOT->SetBatch();
	TFile *inFile = TFile::Open(infilename_.c_str());
	RooWorkspace *inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
	RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hgg_mass");
	PdfModelBuilder pdfsModel;
	pdfsModel.setObsVar(mass);

	TFile *outFile = new TFile(outfilename_.c_str(),"RECREATE");
	RooWorkspace *outWS = new RooWorkspace("multipdf","multipdf");

	// signal specific
	TFile *sigFile = NULL;
	RooWorkspace *sigWS = NULL;
	RooRealVar *MH = NULL;
	RooRealVar *mu = NULL;
	if (doSBversion_) {
		sigFile = TFile::Open(sigfilename_.c_str());
		sigWS = (RooWorkspace*)sigFile->Get(Form("wsig_%dTeV",sqrtS_));
		MH = (RooRealVar*)sigWS->var("MH");
		mu = new RooRealVar("mu","mu",0.01,-5.,10.);
	}

	vector<string> pdfTypes;
	pdfTypes.push_back("Bernstein");
	pdfTypes.push_back("Exponential");
	pdfTypes.push_back("PowerLaw");
	pdfTypes.push_back("Laurent");

	double fTestThreshold = 0.05;
	double envThreshold = 0.1;

	for (int cat=startingCat_; cat<ncats_; cat++) {
		RooAbsData *data = (RooDataSet*)inWS->data(Form("data_mass_cat%d",cat));
		vector<pair<string,int> > fTestResult;
		RooArgList envelopeResult("store");
		RooArgList envelopeResultBonly("bonlystore");
		for (unsigned int iPdf=0; iPdf<pdfTypes.size(); iPdf++){
			bool fTestReady=false;
			bool envelopeReady=false;
			int order=1;
			int prev_order=0;
			int truth_order=-1;
			double thisNll=0.;
			double prevNll=0.;
			while (!fTestReady && !envelopeReady) {
				RooAbsPdf *bkgPdf = getPdf(pdfsModel,pdfTypes[iPdf],order,Form("pdf_cat%d_%dTeV",cat,sqrtS_));
				if (!bkgPdf) {
					order++;
				}
				else {
					RooAbsPdf *thePdf;
					if (doSBversion_) thePdf = getSBPdf(bkgPdf,sigWS,cat,MH,mu,data->sumEntries());
					else thePdf = bkgPdf;
					assert(data);
					assert(thePdf);
					//thePdf->Print();
					//data->Print();
					bkgPdf->SetTitle(Form("%s %d",pdfTypes[iPdf].c_str(),order));
					thePdf->SetTitle(Form("%s %d",pdfTypes[iPdf].c_str(),order));
					RooFitResult *fitRes = thePdf->fitTo(*data,Save(true));
					
					// fTest for truth
					thisNll = 2*fitRes->minNll(); 
					double chi2_nlldiff = prevNll-thisNll;
					if (chi2_nlldiff<0. && order>1) chi2_nlldiff=0.;
					double prob = TMath::Prob(chi2_nlldiff,order-prev_order);

					cout << "FTEST: " << thePdf->GetName() << " prevNll: " << prevNll << " thisNll: " << thisNll << " chi2diff: " << chi2_nlldiff << " chi2prob: " << prob << endl;

					if (prob >= fTestThreshold) {
						truth_order = prev_order;
						fTestResult.push_back(make_pair(pdfTypes[iPdf],truth_order));
						fTestReady=true;
					}
					
					// gof and fTest for envelope
					int nFreeParams = bkgPdf->getVariables()->getSize()-1;
					if (doSBversion_) nFreeParams+=3;
					double chi2 = calcChi2(mass,thePdf,data,nFreeParams,80);
					double gofProb = calcChi2prob(mass,thePdf,data,nFreeParams,80);
					
					cout << "ENVELOPE: " << thePdf->GetName() << " gofProb: " << gofProb << endl;

					if ( (gofProb > 0.01 && prob < envThreshold) || order == truth_order) {
						envelopeResult.add(*thePdf);
						envelopeResultBonly.add(*bkgPdf);
					}
					if (prob >= envThreshold) envelopeReady=true;

					// plot it
					if (doSBversion_) plot(mass,thePdf,data,bkgPdf->getParameters(*data),Form("plots/fTest/%s%d_cat%d",pdfTypes[iPdf].c_str(),order,cat),chi2,gofProb,bkgPdf);
					else plot(mass,thePdf,data,bkgPdf->getParameters(*data),Form("plots/fTest/%s%d_cat%d",pdfTypes[iPdf].c_str(),order,cat),chi2,gofProb);

					prev_order = order;
					prevNll = thisNll;
					order++;
				}
			}
		} // end pdf loop
		envelopeResult.Print("v");

		cout << "fTest Result: " << endl;
		for (unsigned int i=0; i<fTestResult.size(); i++){
			cout << "\t" << fTestResult[i].first << " : " << fTestResult[i].second << endl;
		}

		RooCategory catIndex(Form("pdfindex_%d_%dTeV",cat,sqrtS_),"c");
		RooCategory catIndexSB(Form("sbpdfindex_%d_%dTeV",cat,sqrtS_),"c");
		RooMultiPdf *bmpdf = new RooMultiPdf(Form("CMS_hgg_cat%d_%dTeV_bkgshape",cat,sqrtS_),"all pdfs",catIndex,envelopeResultBonly);
		RooMultiPdf *mpdf = new RooMultiPdf(Form("CMS_hgg_cat%d_%dTeV_bkgshape_forBFSetting",cat,sqrtS_),"all pdfs",catIndexSB,envelopeResult);
		RooRealVar nBackground(Form("CMS_hgg_cat%d_%dTeV_bkgshape_norm",cat,sqrtS_),"nbkg",data->sumEntries(),data->sumEntries()-5.*TMath::Sqrt(data->sumEntries()),data->sumEntries()+5.*TMath::Sqrt(data->sumEntries()));

		int bestFitIndex = getBestFitFunction(bmpdf,data,&catIndex,true);
		catIndex.setIndex(bestFitIndex);
		
		/*
		if (doSBversion_) {
			bestFitIndex = getBestFitFunction(mpdf,data,&catIndexSB,true);
			catIndexSB.setIndex(bestFitIndex);
			catIndex.setIndex(bestFitIndex);
		}
		*/
		
		//if (doSBversion_) plot(mass,mpdf,data,&catIndexSB,bestFitIndex,cat,Form("plots/fTest/multipdf_cat%d",cat));

		plot(mass,bmpdf,data,&catIndex,bestFitIndex,cat,Form("plots/fTest/multipdf_cat%d",cat));

		mass->setBins(320);
		RooDataHist dataBinned(Form("roohist_data_mass_cat%d",cat),"data",*mass,*data);

		outWS->import(*data);
		outWS->import(dataBinned);
		outWS->import(nBackground);
		outWS->import(catIndex);
		outWS->import(*bmpdf,RecycleConflictNodes());
		
	}
	
	outFile->cd();
	outWS->Write();
	outFile->Close();

	return 0;
}
