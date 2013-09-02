#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMacro.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"

#include <iostream>

using namespace RooFit;
using namespace std;

/*
void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname);
*/

void doBandsFit(TGraphAsymmErrors *onesigma, TGraphAsymmErrors *twosigma, 
		RooRealVar * hmass,
		RooAbsPdf *cpdf, 
		RooCurve *nomcurve,  RooAbsData *datanorm,
		RooPlot *plot, 
		TString & catname)
{
	RooRealVar *nlim = new RooRealVar(TString::Format("nlim%s",catname.Data()),"",0.0,0.0,1e+5);

	for (int i=1; i<(plot->GetXaxis()->GetNbins()+1); ++i) {

		double lowedge = plot->GetXaxis()->GetBinLowEdge(i);
		double upedge = plot->GetXaxis()->GetBinUpEdge(i);
		double center = plot->GetXaxis()->GetBinCenter(i);
        
		double nombkg = nomcurve->interpolate(center);

		nlim->setVal(nombkg);
		hmass->setRange("errRange",lowedge,upedge);
		RooAbsPdf *epdf = 0;
		epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
        
		RooAbsReal *nll = epdf->createNLL(*datanorm,Extended(),NumCPU(4));
		RooMinimizer minim(*nll);
		minim.setStrategy(0);
		minim.setPrintLevel(-1);
		//double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
		double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
        
		minim.migrad();
		minim.minos(*nlim);
		
		printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
		
		onesigma->SetPoint(i-1,center,nombkg);
		onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
        
		minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
		// eventually if cl = 0.95 this is the usual 1.92!      
        	minim.migrad();
		minim.minos(*nlim);
		
		twosigma->SetPoint(i-1,center,nombkg);
		twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());      
        		
		delete nll;
		delete epdf;
	}
	//onesigma->Print("V");
}

void getConfigFromFile(TFile *inFile, bool &is2011, bool &baseline){
  
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      if (macro->GetName()==TString("statanalysis") || macro->GetName()==TString("spinanalysis")) baseline=true;
      if (macro->GetName()==TString("massfactorizedmvaanalysis")) baseline=false;
      TList *list = (TList*)macro->GetListOfLines();
      for (int l=0; l<list->GetSize(); l++){
        TObjString *obStr = (TObjString*)list->At(l);
        TString line(obStr->GetName());
        TString varName(line(0,line.First("=")));
        TString varVal(line(line.First("=")+1,line.Length()));
        if (varName=="dataIs2011") is2011 = varVal.Atoi();
      }
    }
  }
}

void getNCats(RooWorkspace *ws, int mh, int &ncats){
  ncats=0;
  RooDataSet *data = (RooDataSet*)ws->data(Form("sig_mass_m%d_cat0",mh));
  while (1) {
    data = (RooDataSet*)ws->data(Form("sig_mass_m%d_cat%d",mh,ncats));
    if (!data) break;
    else ncats++;
  }
}

void makeBkgPlotsGeneric(string filebkg, string filesig="", bool blind=true, bool doBands=true, bool useBinnedData=false, int ncats=9, bool is2011=false, bool baseline=false, bool spin=false) {
	
  // Globals
	gROOT->SetStyle("Plain");
	gROOT->SetBatch(1);
	gStyle->SetOptStat(0);
  
  TFile *outf = TFile::Open("BkgPlotCanvs.root","RECREATE");
	
  RooMsgService::instance().setGlobalKillBelow(RooFit::MsgLevel(RooFit::FATAL));
  RooMsgService::instance().setSilentMode(true);
	
  TFile *fb = TFile::Open(filebkg.c_str());
	TFile *fs = ( filesig.empty() ? fb : TFile::Open(filesig.c_str()) );

  getConfigFromFile(fb,is2011,baseline);
	
	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
  RooWorkspace *w_sig = (RooWorkspace*)fs->Get("wsig_8TeV");

  bool doSignal=false;
  if (w_sig) doSignal=true;

	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
  RooRealVar *MH = 0;
	double lumi = intL->getVal()/1000.;

  //getNCats(w_bkg,125,ncats);

	TLatex *latex = new TLatex();	
	latex->SetTextSize(0.03);
	latex->SetNDC();
	TLatex *cmslatex = new TLatex();
	cmslatex->SetTextSize(0.03);
	cmslatex->SetNDC();

	for (int cat=0;cat<ncats;cat++){
		
		TCanvas *can = new TCanvas("c","",800,800);
		TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);

		// Get Dataset ->
		RooAbsData *data;
		if (useBinnedData) data =  (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
		else  data =  (RooDataHist*)w_bkg->data(Form("roohist_data_mass_cat%d",cat));
		data->Print();

		// Background Pdf ->
		/// RooExtendPdf *bkg =  (RooExtendPdf*)w_bkg->pdf(Form("data_pol_model_8TeV_cat%d",cat));
		RooAbsPdf *bkg =  (RooAbsPdf*)w_bkg->pdf(Form("pdf_data_pol_model_8TeV_cat%d",cat));
		bkg->Print();
		bkg->fitTo(*data);
		RooFitResult *r = bkg->fitTo(*data,RooFit::Save(1));
    r->Print();

    // Get signal mc 
    RooDataSet *sigMC = 0; 
    TH1F *sigMCHist = 0; 
    if (doSignal) {
      MH = (RooRealVar*)w_sig->var("MH");
      sigMC = (RooDataSet*)w_sig->data(Form("sig_mass_m125_cat%d",cat));
      sigMCHist = new TH1F("sigHist","sigHist",80,100,180);
      sigMC->fillHistogram(sigMCHist,RooArgList(*x));
      sigMCHist->SetLineColor(4);
      sigMCHist->SetFillColor(38);
      sigMCHist->SetFillStyle(3001);
      sigMCHist->SetLineWidth(2);
    }

    // Get signal pdf
    RooAbsPdf *sigPDF = 0;
    if (doSignal) sigPDF = (RooAbsPdf*)w_sig->pdf(Form("sigpdfrelcat%d_allProcs",cat));

    // Plotting
    RooPlot *frame = x->frame();
    data->plotOn(frame,Binning(80),Invisible());
    TObject *dataLeg = (TObject*)frame->getObject(frame->numItems()-1);
    leg->AddEntry(dataLeg,"Data","LEP");
    bkg->plotOn(frame,LineColor(kRed));
    TObject *bkgLeg = (TObject*)frame->getObject(frame->numItems()-1);
    leg->AddEntry(bkgLeg,"Bkg model","L");
	  
    // Get error bands
    TGraphAsymmErrors *onesigma = 0, *twosigma = 0;
		if( doBands ) {
			onesigma = new TGraphAsymmErrors();
 			onesigma->SetLineColor(kYellow);
 			onesigma->SetFillColor(kYellow);
 			onesigma->SetMarkerColor(kYellow);
			twosigma = new TGraphAsymmErrors();
 			twosigma->SetLineColor(kGreen);
 			twosigma->SetFillColor(kGreen);
 			twosigma->SetMarkerColor(kGreen);
			TString name=Form("cat%d",cat);
      doBandsFit(onesigma, twosigma, x, bkg, dynamic_cast<RooCurve*>(frame->getObject(frame->numItems()-1)),data, frame, name );
      leg->AddEntry(onesigma,"#pm 1#sigma Bkg model","LF");
      leg->AddEntry(twosigma,"#pm 2#sigma Bkg model","LF");
		}
   
    // do plotting
    if (blind) {
			x->setRange("unblind_up",150,180);
			data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_up"));
			x->setRange("unblind_down",100,110);
			data->plotOn(frame,RooFit::Binning(80),RooFit::CutRange("unblind_down"));
    }
    else {
      data->plotOn(frame,Binning(80));
    }

    if (doSignal) {
      MH->setVal(125);
      sigPDF->plotOn(frame,Normalization(sigMCHist->Integral(),RooAbsReal::NumEvent),LineColor(kBlue),LineWidth(3));
      TObject *sigLeg = (TObject*)frame->getObject(frame->numItems()-1);
      leg->AddEntry(sigLeg,"Sig model m_{H}=125GeV","L");
    }

    frame->SetTitle("");
    frame->SetMinimum(0.0001);
		frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
		frame->GetXaxis()->SetNdivisions(5,5,0);
		frame->GetYaxis()->SetTitle("Events / (1 GeV)");
		frame->GetYaxis()->SetTitleOffset(1.2);
    
    frame->Draw();
 		if( doBands ) {
 			twosigma->Draw("L3 SAME");     
 			onesigma->Draw("L3 SAME");     
 			frame->Draw("same");
 		}
    //if (doSignal) sigMCHist->Draw("FHISTsame");
		leg->Draw();
		cmslatex->DrawLatex(0.2,0.85,Form("#splitline{CMS Preliminary}{#sqrt{s} = 8TeV L = %2.1ffb^{-1}}",lumi));
		//latex->DrawLatex(0.2,0.75,labels[cat].c_str());

    can->SaveAs(Form("newplotcat%d.pdf",cat));
  
  }
  outf->Close();
  fb->Close();
		
}


