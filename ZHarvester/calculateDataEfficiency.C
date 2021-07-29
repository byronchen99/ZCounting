#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TTree.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TEfficiency.h>            // class to handle efficiency calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TPad.h>
#include <TLatex.h>

#include "Utils/CPlot.hh"             // helper class for plots
#include "Utils/MitStyleRemix.hh"         // style settings for drawing
#include "Utils/CEffUser1D.hh"            // class for handling efficiency graphss
#include "Utils/CEffUser2D.hh"            // class for handling efficiency tables

#include "Utils/ZMMSignals.hh"
#include "Utils/ZMMBackgrounds.hh"
#include "Utils/RooGaussDoubleSidedExp.h"
#endif

// RooFit headers
#include "RooWorkspace.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"

// constants
const Float_t etaCutTag   = 2.4;
const Float_t etaCutProbe = 2.4;
const Float_t etaBound    = 0.9;

// Set up pile-up bounds available get MC template for fitting
const Int_t minPU = 1;
const Int_t maxPU = 60;


// Global settings
//--------------------------------------------------------------------------------------------------
Float_t massLo  = 56.;
Float_t massHi  = 116.;
UInt_t  massBin = (UInt_t)(massHi-massLo);

Float_t ptCutTag = 30.;
Float_t ptCutProbe = 30.;

TString outputDir="";

char lumitext[100] = "";

// Global setters
//--------------------------------------------------------------------------------------------------
void set_massRange(Float_t massLo_, Float_t massHi_, UInt_t nBins=0){
    massLo = massLo_;
    massHi = massHi_;
    if(nBins == 0){
        massBin = (UInt_t)(massHi-massLo);
    }
    else{
        massBin = nBins;
    }

    std::cout<<"Set mass range to ["<<massLo_<<","<<massHi_<<"] with "<<massBin<<" bins."<<std::endl;
}

void set_ptCut(Float_t pt_){
    ptCutTag = pt_;
    ptCutProbe = pt_;
    std::cout<<"Set pT cut of tag and probe muons to "<<pt_<<" GeV"<<std::endl;
}

void set_output(TString dir_){
    outputDir = dir_;
    std::cout<<"Set output directory to "<<outputDir<<std::endl;
}

void set_luminosity(Float_t lumi_){
    sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi_);
}

//--------------------------------------------------------------------------------------------------
// generate template for extraction of muon efficiency in barrel or endcap region
void generateTemplate(
    const TString mcfilename,
    TH1D *hPV=0
);

//--------------------------------------------------------------------------------------------------
// generate template for extraction of Z yield in {BB, BE, EE} region
void generateTemplate_ZYield(
	const TString mcfilename,
	TH1D          *hPV
);

//--------------------------------------------------------------------------------------------------
// count events in pass and fail region to measure the efficiency
void performCount(
    Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
    const TString effType, const Bool_t etaRegion, const Int_t iBin, const TString format
);

//--------------------------------------------------------------------------------------------------
// perform fit in 1 region to extract the number of signal events
void performFit(
	Double_t &resNsig, Double_t &resErrl, Double_t &resErrh,
    Double_t &resChi2, Double_t &resPurity, Double_t &resPurityErrl, Double_t &resPurityErrh,
    TH1D *h, const Int_t sigMod, const Int_t bkgMod,
	const TString etaRegion, const Bool_t passRegion,
    const Int_t iBin, const TString format, const Int_t fitStrategy=0
);

//--------------------------------------------------------------------------------------------------
// perform fit in 2 region (pass and fail) to extract the efficiency
void performFit(
    Double_t &resEff, Double_t &resErrl, Double_t &resErrh,
    Double_t &resChi2Pass, Double_t &resChi2Fail,
    TH1D *passHist, TH1D *failHist,
    const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
    const TString bkgQCDTemplate, const TString bkgTTTemplate,
    const TString effType, const Bool_t etaRegion, const Int_t iBin, const TString format, const Int_t fitStrategy=0
);

//--------------------------------------------------------------------------------------------------
// perform fit in all regions together to extract the HLT efficiency and the number of signal events
// BB: both muons are in barrel
// BE: one muon is in barrel, one in endcap
// EE: both muons are in endcap
// passHist{BB,BE,EE}:  dimuon events where both muons pass the HLT
// failHist{BB,BE,EE}:  dimuon events where only one muon passes the HLT
// c{BB,BE,EE}:         correlation coefficient for the second muon passing the HLT, if the first has already passed the HLT
void performCombinedFit(
    Double_t &resEffB, Double_t &resErrBl, Double_t &resErrBh,
    Double_t &resEffE, Double_t &resErrEl, Double_t &resErrEh,
    Double_t &resYieldBB, Double_t &resYieldBBl, Double_t &resYieldBBh,
    Double_t &resYieldEE, Double_t &resYieldEEl, Double_t &resYieldEEh,
    Double_t &resYieldBE, Double_t &resYieldBEl, Double_t &resYieldBEh,
    Double_t &resChi2PassBB, Double_t &resChi2FailBB,
    Double_t &resChi2PassEE, Double_t &resChi2FailEE,
    Double_t &resChi2PassBE, Double_t &resChi2FailBE,
    TH1D *passHistBB, TH1D *failHistBB,
    TH1D *passHistEE, TH1D *failHistEE,
    TH1D *passHistBE, TH1D *failHistBE,
    Double_t cBB, Double_t cBE, Double_t cEE,
    const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
    const TString bkgQCDTemplate, const TString bkgTTTemplate,
    const Int_t iBin,
    const TString format,
    const Int_t fitStrategy=0
);

std::vector<double> preFit(TH1D *failHist);

//--------------------------------------------------------------------------------------------------
// Signal Model:
// 	    1: Breit-Wigner convolved with Crystal Ball function
// 	    2: MC template convolved with Gaussian
Int_t set_signal_model(
    Int_t model_type,
    CSignalModel *&model,
    RooRealVar &param_mass,
    const Bool_t pass,
    Int_t ibin=0,
    TH1D *hist=0
){

    switch(model_type) {
        case 1:
            model = new CBreitWignerConvCrystalBall(param_mass, pass, ibin);
            return 4;
        case 2:
            model = new CMCTemplateConvGaussian(param_mass, hist, pass, ibin);
            return 2;
    }
    return 0;
}

//--------------------------------------------------------------------------------------------------
// Background Model:
//      1: exponential model
//      2: quadratic model
//      3: exponential + quadratic model
//      4: Das (wide gaussian with exponential tails) model
//      5: Das + decaying exponential model
//      6: bkg QCD template
//      (7: bkg QCD + ttbar template)
Int_t set_background_model(
    Int_t model_type,
    CBackgroundModel *&model,
    RooRealVar &param_mass,
    const Bool_t pass,
    Int_t ibin=0,
    TH1D *hist=0,
    std::vector<double> *params=0
){
    switch(model_type) {
        case 1:
            model = new CExponential(param_mass, pass, ibin);
            return 1;
        case 2:
            if(params != 0)
                model = new CQuadratic(param_mass, pass, ibin,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadratic(param_mass, pass, 0, 0.,0.,0.,0.,0.,0.);
            return 3;
        case 3:
            if(params != 0)
                model = new CQuadPlusExp(param_mass, pass, ibin,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadPlusExp(param_mass, pass, ibin, 0.,0.,0.,0.,0.,0.);
                return 4;
        case 4:
            model = new CDas(param_mass, pass, ibin);
            return 4;
        case 5:
            model = new CDasPlusExp(param_mass, pass, ibin);
            return 6;
        case 6:
            model = new CRooCMSShape(param_mass, pass, ibin);
            return 3;
        case 7:
            model = new CQCD(param_mass, hist, pass, ibin);
            return 1;
    }
    return 0;
}

//--------------------------------------------------------------------------------------------------
void generateTemplate(
	const TString mcfilename,
	TH1D          *hPV
){
  cout << "Creating histogram templates... "; cout.flush();

  TFile *infile    = new TFile(mcfilename);
  TTree *eventTree = (TTree*)infile->Get("tree");
  TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

  if(hPV)
    hPV->Divide(hPVtemplate);

  Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
  Double_t wgt;
  UInt_t npv;
  Bool_t pass;

  eventTree->SetBranchAddress("mass",       &mass);
  eventTree->SetBranchAddress("ptTag",      &ptTag);
  eventTree->SetBranchAddress("ptProbe",    &ptProbe);
  eventTree->SetBranchAddress("etaTag",     &etaTag);
  eventTree->SetBranchAddress("etaProbe",   &etaProbe);
  eventTree->SetBranchAddress("nPV",        &npv);
  eventTree->SetBranchAddress("pass",       &pass);

  TH1D *h_mass_pass_central = new TH1D("h_mass_pass_central", "", massBin, massLo, massHi);
  TH1D *h_mass_fail_central = new TH1D("h_mass_fail_central", "", massBin, massLo, massHi);
  TH1D *h_mass_pass_forward = new TH1D("h_mass_pass_forward", "", massBin, massLo, massHi);
  TH1D *h_mass_fail_forward = new TH1D("h_mass_fail_forward", "", massBin, massLo, massHi);

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    eventTree->GetEntry(ientry);

    if(mass < massLo)  continue;
    if(mass > massHi)  continue;
    if(ptTag   < ptCutTag)   continue;
    if(ptProbe   < ptCutProbe)   continue;
    if(fabs(etaTag) > etaCutTag) continue;
    if(fabs(etaProbe) > etaCutProbe) continue;

    wgt = 1.0;
    if(hPV)
      wgt *= hPV->GetBinContent(hPV->FindBin(npv));

    if(fabs(etaProbe) < etaBound){
      if(pass) h_mass_pass_central->Fill(mass, wgt);
      else     h_mass_fail_central->Fill(mass, wgt);
    }else{
      if(pass) h_mass_pass_forward->Fill(mass, wgt);
      else     h_mass_fail_forward->Fill(mass, wgt);
    }
  }

  TFile outfile = TFile(outputDir+"/histTemplates.root", "RECREATE");
  h_mass_pass_central->Write();
  h_mass_fail_central->Write();
  h_mass_pass_forward->Write();
  h_mass_fail_forward->Write();
  outfile.Write();
  outfile.Close();

  infile->Close();
  delete infile;

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
void generateTemplate_ZYield(
	const TString mcfilename,
	TH1D          *hPV
){
  cout << "Creating histogram templates... "; cout.flush();

  TFile *infile    = new TFile(mcfilename);
  TTree *eventTree = (TTree*)infile->Get("tree");
  TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

  if(hPV)
    hPV->Divide(hPVtemplate);

  Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
  Double_t wgt;
  UInt_t npv;
  Bool_t pass;

  eventTree->SetBranchAddress("mass",       &mass);
  eventTree->SetBranchAddress("ptTag",      &ptTag);
  eventTree->SetBranchAddress("ptProbe",    &ptProbe);
  eventTree->SetBranchAddress("etaTag",     &etaTag);
  eventTree->SetBranchAddress("etaProbe",   &etaProbe);
  eventTree->SetBranchAddress("nPV",        &npv);
  eventTree->SetBranchAddress("pass",       &pass);

  TH1D *h_mass_2HLT_BB = new TH1D("h_mass_2HLT_BB", "", massBin, massLo, massHi);
  TH1D *h_mass_2HLT_BE = new TH1D("h_mass_2HLT_BE", "", massBin, massLo, massHi);
  TH1D *h_mass_2HLT_EE = new TH1D("h_mass_2HLT_EE", "", massBin, massLo, massHi);
  TH1D *h_mass_1HLT_BB = new TH1D("h_mass_1HLT_BB", "", massBin, massLo, massHi);
  TH1D *h_mass_1HLT_BE = new TH1D("h_mass_1HLT_BE", "", massBin, massLo, massHi);
  TH1D *h_mass_1HLT_EE = new TH1D("h_mass_1HLT_EE", "", massBin, massLo, massHi);

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    eventTree->GetEntry(ientry);

    if(mass < massLo)  continue;
    if(mass > massHi)  continue;
    if(ptTag   < ptCutTag)   continue;
    if(ptProbe   < ptCutProbe)   continue;
    if(fabs(etaTag) > etaCutTag) continue;
    if(fabs(etaProbe) > etaCutProbe) continue;

    wgt = 1.0;
    if(hPV)
      wgt *= hPV->GetBinContent(hPV->FindBin(npv));

    if(fabs(etaProbe) < etaBound && fabs(etaTag) < etaBound){
        if(pass) h_mass_2HLT_BB->Fill(mass, wgt);
        else     h_mass_1HLT_BB->Fill(mass, wgt);
    }
    else if(fabs(etaProbe) >= etaBound && fabs(etaTag) >= etaBound){
        if(pass) h_mass_2HLT_EE->Fill(mass, wgt);
        else     h_mass_1HLT_EE->Fill(mass, wgt);
    }
    else{
        if(pass) h_mass_2HLT_BE->Fill(mass, wgt);
        else     h_mass_1HLT_BE->Fill(mass, wgt);
    }

  }

  TFile outfile = TFile(outputDir+"/histTemplates.root", "RECREATE");
  h_mass_2HLT_BB->Write();
  h_mass_2HLT_BE->Write();
  h_mass_2HLT_EE->Write();
  h_mass_1HLT_BB->Write();
  h_mass_1HLT_BE->Write();
  h_mass_1HLT_EE->Write();
  outfile.Write();
  outfile.Close();

  infile->Close();
  delete infile;

  cout << "Done!" << endl;
}


//--------------------------------------------------------------------------------------------------
template<typename T>
Double_t make_plot(
    const Int_t   nfl,
    const RooRealVar &param_mass,
    RooAbsData *data,
    const RooAbsPdf &modelPdf,
    CSignalModel *sigModel,
    CBackgroundModel *bkgModel,
    RooFitResult *fitResult,
    T &Nsig,
    RooRealVar &Nbkg,
    const char effstr[],
    const Int_t nEntries,
    const Int_t iBin,
    const TString effType="yield",
    const Bool_t etaRegion=kTRUE,
    const Bool_t passRegion=kTRUE,
    const Double_t yMax=0.,
    const TString format="png"
){

    char pname[50];
    char ctitle[100];
    char binlabelx[100];
    char binlabely[100];
    char ylabel[50];
    char yield[50];
    char nsigstr[100];
    char nbkgstr[100];
    char chi2str[100];

    if(effType == "BB" || effType == "BE" || effType == "EE"){
        sprintf(pname,"yield_%s_%s_%d", effType.Data(), passRegion ? "2HLT" : "1HLT", iBin);
        sprintf(ctitle,"%s %s %d ", passRegion ? "2HLT" : "1HLT", effType.Data(), iBin);
        if(effType == "BB")
            sprintf(binlabelx, "|#eta| < %.1f",etaBound);
        if(effType == "BE")
            sprintf(binlabelx, "|#eta| < 2.4");
        if(effType == "EE")
            sprintf(binlabelx, "%.1f < |#eta| < 2.4",etaBound);
    }
    else{
        sprintf(pname,"%s_%s_%s_%d", effType.Data(), etaRegion ? "forward" : "central", passRegion ? "pass" : "fail", iBin);
        sprintf(ctitle,"%s %s (%d) ", effType.Data(), passRegion ? "pass" : "fail", iBin);

        if(!etaRegion) sprintf(binlabelx, "|#eta| < 0.9");
        else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");
    }

    const double margin_left = 0.15;
    const double margin_right = 0.03;

    TCanvas *canvas = MakeCanvas(pname, ctitle,800,800);

    TPad *pad1 = new TPad("pad1", "pad1", 0., 0.35, 1, 1.0);

    canvas->SetTicks();
    pad1->SetLeftMargin(margin_left);
    pad1->SetRightMargin(margin_right);
    pad1->SetTopMargin(0.01);
    pad1->SetBottomMargin(0.025);
    pad1->SetTickx();
    pad1->SetTicky();
    pad1->Draw();
    pad1->cd();

    const double textsize1 = 28./(pad1->GetWh()*pad1->GetAbsHNDC());

    TLegend *legend = new TLegend(1-margin_right-0.6, 0.63, 1-margin_right-0.05, 0.79);
    legend->SetNColumns(2);

    TLegendEntry *entry1 = new TLegendEntry();
    entry1->SetTextSize(textsize1);
    entry1->SetLabel("Data");
    entry1->SetOption("PE");
    entry1->SetMarkerStyle(kFullCircle);
    entry1->SetMarkerColor(kBlack);
    entry1->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry1);

    TLegendEntry *entry4 = new TLegendEntry();
    entry4->SetTextSize(textsize1);
    entry4->SetLabel("Full model");
    entry4->SetOption("l");
    entry4->SetLineColor(2);
    entry4->SetLineWidth(2);
    entry4->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry4);

    TLegendEntry *entry2 = new TLegendEntry();
    entry2->SetTextSize(textsize1);
    entry2->SetLabel(bkgModel->model->GetTitle());
    entry2->SetOption("l");
    entry2->SetLineColor(8);
    entry2->SetLineWidth(2);
    entry2->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry2);

    TLegendEntry *entry3 = new TLegendEntry();
    entry3->SetTextSize(textsize1);
    entry3->SetLabel(sigModel->model->GetTitle());
    entry3->SetOption("l");
    entry3->SetLineColor(9);
    entry3->SetLineWidth(2);
    entry3->SetTextAlign(12);
    legend->GetListOfPrimitives()->Add(entry3);

    RooPlot *mframe = param_mass.frame(Bins(massBin));
    data->plotOn(mframe, MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelPdf.plotOn(mframe, Components(*(bkgModel->model)), LineColor(8));
    modelPdf.plotOn(mframe, Components(*(sigModel->model)), LineColor(9));
    modelPdf.plotOn(mframe, LineColor(kRed));

    // Construct a histogram with the pulls of the data w.r.t the curve
    RooHist *hpull = mframe->pullHist();

    double a = Nsig.getVal(), aErr = Nsig.getPropagatedError(*fitResult);
    double b = Nbkg.getVal(), bErr = Nbkg.getPropagatedError(*fitResult);
    double resChi2 = mframe->chiSquare(nfl);

    sprintf(binlabely, "p_{T} > %i GeV",(Int_t)ptCutTag);
    sprintf(ylabel,"Events / 1 GeV");
    sprintf(yield,"%u Events",nEntries);
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",a,aErr);
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",resChi2);
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",b,bErr);

    mframe->SetTitle("");

    mframe->GetXaxis()->SetTitleSize(0.);
    mframe->GetXaxis()->SetLabelSize(0.);

    mframe->GetYaxis()->SetTitle(ylabel);
    mframe->GetYaxis()->SetLabelSize(textsize1);
    mframe->GetYaxis()->SetTitleFont(42);
    mframe->GetYaxis()->SetTitleSize(textsize1*1.2);
    mframe->GetYaxis()->SetTitleOffset(1.1);

    mframe->SetMinimum(0.);
    mframe->SetMaximum(mframe->GetMaximum()*1.2);

    mframe->Draw();

    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(11);
    latex->SetTextFont(62);
    latex->SetTextSize(textsize1*1.2);
    latex->DrawLatex(margin_left+0.04, 0.91, "CMS");
    latex->SetTextFont(52);
    latex->SetTextSize(textsize1);
    latex->DrawLatex(margin_left+0.14, 0.91, "Work in progress");
    latex->SetTextFont(42);
    latex->DrawLatex(margin_left+0.04, 0.84, lumitext);

    latex->DrawLatex(1-margin_right-0.6, 0.54, ctitle);
    latex->DrawLatex(1-margin_right-0.6, 0.47, yield);
    latex->DrawLatex(1-margin_right-0.35, 0.54, nsigstr);
    latex->DrawLatex(1-margin_right-0.35, 0.47, nbkgstr);
    latex->DrawLatex(1-margin_right-0.35, 0.40, chi2str);
    latex->DrawLatex(1-margin_right-0.35, 0.33, effstr);

    latex->SetTextAlign(31);
    latex->DrawLatex(1-margin_right-0.04, 0.91, binlabelx);
    latex->DrawLatex(1-margin_right-0.04, 0.84, binlabely);

    legend->Draw("same");

    canvas->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.35);
    pad2->SetLeftMargin(margin_left);
    pad2->SetRightMargin(margin_right);
    pad2->SetTopMargin(0.025);
    pad2->SetBottomMargin(0.3);
    pad2->SetTickx();
    pad2->SetTicky();
    pad2->Draw("ALPF");
    pad2->cd();

    const double textsize2 = 28./(pad2->GetWh()*pad2->GetAbsHNDC());

    RooPlot *rframe = param_mass.frame(Bins(massBin));

    rframe->SetTitle("");
    rframe->addPlotable(hpull, "PX");

    rframe->GetYaxis()->SetTitle("Pulls");
    rframe->GetYaxis()->SetTitleSize(textsize2*1.2);
    rframe->GetYaxis()->SetTitleOffset(0.6);
    rframe->GetYaxis()->SetLabelSize(textsize2);
    rframe->GetYaxis()->SetTitleFont(42);
    rframe->GetYaxis()->CenterTitle(true);
    rframe->GetYaxis()->SetNdivisions(405);

    rframe->GetXaxis()->SetTitle("tag-and-probe mass [GeV]");
    rframe->GetXaxis()->SetTitleSize(textsize2*1.2);
    rframe->GetXaxis()->SetTitleOffset(1.);
    rframe->GetXaxis()->SetLabelSize(textsize2);
    rframe->GetXaxis()->SetTitleFont(42);
    rframe->Draw();

    TLine *line0 = new TLine(massLo, 0., massHi, 0.);
    line0->SetLineStyle(1);
    line0->Draw("same");

    canvas->SaveAs(outputDir+"/"+pname+".png");
    canvas->SaveAs(outputDir+"/"+pname+".eps");
    canvas->Close();

    return resChi2;
}

//--------------------------------------------------------------------------------------------------
std::vector<float> getZyield(
                TH1D *h_yield,
                const Int_t    iBin,         // Label of measurement in currect run
                const Int_t    sigMod,
                const Int_t    bkgMod,
                const TString  etaRegion="",    // {BB, BE or EE}
                const Bool_t   passRegion=true,
                const TString  mcfilename="",
                TH1D*          _hQCD=0,
                TH1D*          hPV=0,
                const TString format="png"          // plot format
){

  CPlot::sOutDir = outputDir;

  if(sigMod == 0){      // perform count - expect 1% fakes
    const Double_t NsigMax = h_yield->Integral();
    std::vector<float> resultEff = {};
    resultEff.push_back(NsigMax*0.99);
    resultEff.push_back(std::sqrt(NsigMax)*0.99);
    resultEff.push_back(std::sqrt(NsigMax)*0.99);
    resultEff.push_back(0.);
    resultEff.push_back(0.01);
    return resultEff;
  }


  TH1D *h=0;
  if(sigMod==2) {
    TFile *histfile = 0;
    generateTemplate_ZYield(mcfilename, hPV);
  }

  Double_t resNsig;
  Double_t resErrl;
  Double_t resErrh;
  Double_t resChi2;
  Double_t resPurity;
  Double_t resPurityErrl;
  Double_t resPurityErrh;

  int i = 0;

  Double_t best_Chi2 = -1.;
  int best_strategy = 0;

  do{
      performFit(resNsig, resErrl, resErrh, resChi2, resPurity, resPurityErrl, resPurityErrh, h_yield,
        sigMod, bkgMod, etaRegion, passRegion, iBin, format, i);

      std::cout<<"---------------------------------------"<<std::endl;
      std::cout<<"------ Fit ("<<i<<") ------------------"<<std::endl;
      std::cout<<"------ chi2 = " << resChi2 <<std::endl;
      std::cout<<"---------------------------------------"<<std::endl;

      if(best_Chi2 < 0 || best_Chi2 > resChi2){
        best_Chi2 = resChi2;
        best_strategy = i;
      }

      i++;

  } while(i < 3 && (best_Chi2 > 5. || best_Chi2 < 0.));

  if(best_strategy != i-1){
      // make sure the strategy that gives the best chi2 is taken
      performFit(resNsig, resErrl, resErrh, resChi2, resPurity, resPurityErrl, resPurityErrh, h_yield,
        sigMod, bkgMod, etaRegion, passRegion, iBin, format, best_strategy);
  }

  std::vector<float> resultEff = {};

  resultEff.push_back(resNsig);
  resultEff.push_back(resErrl);
  resultEff.push_back(resErrh);
  resultEff.push_back(resChi2);
  resultEff.push_back(resPurity);
  resultEff.push_back(resPurityErrl);
  resultEff.push_back(resPurityErrh);

  return resultEff;
}

//--------------------------------------------------------------------------------------------------
std::vector<float> calculateDataEfficiency(
        TH1D          *h_mass_pass,
        TH1D          *h_mass_fail,
		const Int_t   iBin,                 // Label of measurement in currect run
		const TString effType,              // "HLT" or "SIT" or "Glo" or "Sta" or "Trk"
        const Bool_t  etaRegion,
		const Int_t   sigModPass,           // signal extraction method for PASS sample
		const Int_t   bkgModPass,           // background model for PASS sample
		const Int_t   sigModFail,           // signal extraction method for FAIL sample
		const Int_t   bkgModFail,           // background model for FAIL sample
        TH1D          *hPV=0,
        const TString mcfilename="",        // ROOT file containing MC events to generate templates from
        const TString bkgQCDFilename="",    // ROOT file containing bkg template
        const TString bkgTTFilename="",
		const TString format="png"          // plot format
){

  // CPlot::sOutDir = outputDir;

  if(etaRegion)
    std::cout<<">>> Do fit in B for "<<effType<<std::endl;
  else
    std::cout<<">>> Do fit in E for "<<effType<<std::endl;

  Double_t eff  = 0.;
  Double_t errl = 0.;
  Double_t errh = 0.;
  Double_t chi2pass = 999.;
  Double_t chi2fail = 999.;

  if(sigModPass == 0){
    performCount(eff, errl, errh, h_mass_pass, h_mass_fail,
      effType, etaRegion, iBin, format);
  }else{
    // Generate histogram templates from MC if necessary
    if(sigModPass==2 || sigModFail==2) {
        generateTemplate(mcfilename, hPV);
    }
    int i = 0;
    int best_strategy = 0;
    double best_Chi2 = -1;
    double resChi2 = -1;
    do {
        performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail,
            sigModPass, bkgModPass, sigModFail, bkgModFail,
            bkgQCDFilename, bkgTTFilename,
            effType, etaRegion, iBin, format, i);

        std::cout<<"---------------------------------------"<<std::endl;
        std::cout<<"------ Fit ("<<i<<") ------------------"<<std::endl;
        std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
        std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
        std::cout<<"---------------------------------------"<<std::endl;

        resChi2 = std::max(chi2pass, chi2fail);
        if(best_Chi2 < 0 || (chi2pass > 0 && chi2fail > 0 && best_Chi2 > resChi2)){
          best_Chi2 = resChi2;
          best_strategy = i;
        }

        i++;

    } while(i < 3 && (resChi2 > 5. || chi2fail < 0. || chi2pass < 0.));

    if(best_strategy != i-1){
        // make sure the strategy that gives the best chi2 is taken
        performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail,
            sigModPass, bkgModPass, sigModFail, bkgModFail,
            bkgQCDFilename, bkgTTFilename,
            effType, etaRegion, iBin, format, best_strategy);
    }
  }



  std::vector<float> resultEff = {};

  resultEff.push_back(eff);
  resultEff.push_back(errl);
  resultEff.push_back(errh);
  resultEff.push_back(chi2pass);
  resultEff.push_back(chi2fail);

  return resultEff;
}

//--------------------------------------------------------------------------------------------------
std::vector<float> calculateHLTEfficiencyAndYield(
    TH1D          *h_mass_BB_pass,
    TH1D          *h_mass_BB_fail,
    TH1D          *h_mass_EE_pass,
    TH1D          *h_mass_EE_fail,
    TH1D          *h_mass_BE_pass,
    TH1D          *h_mass_BE_fail,
	const Int_t   iBin,                 // Label of measurement in currect run
	const Int_t   sigModPass,           // signal extraction method for PASS sample
	const Int_t   bkgModPass,           // background model for PASS sample
	const Int_t   sigModFail,           // signal extraction method for FAIL sample
	const Int_t   bkgModFail,           // background model for FAIL sample
    const Double_t cBB = 1.,            // correlation factors for second muon
    const Double_t cBE = 1.,            // correlation factors for second muon
    const Double_t cEE = 1.,            // correlation factors for second muon
    TH1D          *hPV=0,
    const TString mcfilename="",        // ROOT file containing MC events to generate templates from
    const TString bkgQCDFilename="",    // ROOT file containing bkg template
    const TString bkgTTFilename="",
	const TString format="png"          // plot format
){


    Double_t effB  = 0.;
    Double_t errBl = 0.;
    Double_t errBh = 0.;
    Double_t effE  = 0.;
    Double_t errEl = 0.;
    Double_t errEh = 0.;
    Double_t yieldBB  = 0.;
    Double_t yieldBBl = 0.;
    Double_t yieldBBh = 0.;
    Double_t yieldEE  = 0.;
    Double_t yieldEEl = 0.;
    Double_t yieldEEh = 0.;
    Double_t yieldBE  = 0.;
    Double_t yieldBEl = 0.;
    Double_t yieldBEh = 0.;
    Double_t chi2passBB = 999.;
    Double_t chi2failBB = 999.;
    Double_t chi2passEE = 999.;
    Double_t chi2failEE = 999.;
    Double_t chi2passBE = 999.;
    Double_t chi2failBE = 999.;

    // Generate histogram templates from MC if necessary
    if(sigModPass==2 || sigModFail==2) {
        generateTemplate_ZYield(mcfilename, hPV);
    }
    int i = 0;
    int best_strategy = 0;
    double best_Chi2 = -1;
    double resChi2 = -1;

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ Fit ("<<i<<") ------------------"<<std::endl;
    std::cout<<"------ c(BB) = " << cBB <<" ---------"<<std::endl;
    std::cout<<"------ c(BE) = " << cBE <<" ---------"<<std::endl;
    std::cout<<"------ c(EE) = " << cEE <<" ---------"<<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    // do {
        performCombinedFit(
            effB, errBl, errBh,
            effE, errEl, errEh,
            yieldBB, yieldBBl, yieldBBh,
            yieldEE, yieldEEl, yieldEEh,
            yieldBE, yieldBEl, yieldBEh,
            chi2passBB, chi2failBB,
            chi2passEE, chi2failEE,
            chi2passBE, chi2failBE,
            h_mass_BB_pass, h_mass_BB_fail,
            h_mass_EE_pass, h_mass_EE_fail,
            h_mass_BE_pass, h_mass_BE_fail,
            cBB, cBE, cEE,
            sigModPass, bkgModPass, sigModFail, bkgModFail,
            bkgQCDFilename, bkgTTFilename,
            iBin,
            format,
            i                                   // strategy
        );

        std::cout<<"---------------------------------------"<<std::endl;
        std::cout<<"------ Fit ("<<i<<") ------------------"<<std::endl;
        std::cout<<"------ eff(B) = "   << effB   << " -"<<errBl   <<" +"<< errBh   << std::endl;
        std::cout<<"------ eff(E) = "   << effE   << " -"<<errEl   <<" +"<< errEh   << std::endl;
        std::cout<<"------ yield(BB) = " << yieldBB << " -"<<yieldBBl <<" +"<< yieldBBh << std::endl;
        std::cout<<"------ yield(EE) = " << yieldEE << " -"<<yieldEEl <<" +"<< yieldEEh << std::endl;
        std::cout<<"------ yield(BE) = " << yieldBE << " -"<<yieldBEl <<" +"<< yieldBEh << std::endl;
        std::cout<<"------ chi2pass(BB) = " << chi2passBB <<std::endl;
        std::cout<<"------ chi2fail(BB) = " << chi2failBB <<std::endl;
        std::cout<<"------ chi2pass(EE) = " << chi2passEE <<std::endl;
        std::cout<<"------ chi2fail(EE) = " << chi2failEE <<std::endl;
        std::cout<<"------ chi2pass(BE) = " << chi2passBE <<std::endl;
        std::cout<<"------ chi2fail(BE) = " << chi2failBE <<std::endl;
        std::cout<<"---------------------------------------"<<std::endl;

    //     resChi2 = std::max(chi2pass, chi2fail);
    //     if(best_Chi2 < 0 || (chi2pass > 0 && chi2fail > 0 && best_Chi2 > resChi2)){
    //       best_Chi2 = resChi2;
    //       best_strategy = i;
    //     }
    //
    //     i++;
    //
    // } while(i < 3 && (resChi2 > 5. || chi2fail < 0. || chi2pass < 0.));

    // if(best_strategy != i-1){
    //     // make sure the strategy that gives the best chi2 is taken
    //     performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail,
    //         sigModPass, bkgModPass, sigModFail, bkgModFail,
    //         bkgQCDFilename, bkgTTFilename,
    //         effType, etaRegion, iBin, format, best_strategy);
    // }
  // }


  std::vector<float> resultEff = {};

  resultEff.push_back(effB);
  resultEff.push_back(errBl);
  resultEff.push_back(errBh);
  resultEff.push_back(effE);
  resultEff.push_back(errEl);
  resultEff.push_back(errEh);
  resultEff.push_back(yieldBB);
  resultEff.push_back(yieldBBl);
  resultEff.push_back(yieldBBh);
  resultEff.push_back(yieldEE);
  resultEff.push_back(yieldEEl);
  resultEff.push_back(yieldEEh);
  resultEff.push_back(yieldBE);
  resultEff.push_back(yieldBEl);
  resultEff.push_back(yieldBEh);
  resultEff.push_back(chi2passBB);
  resultEff.push_back(chi2failBB);
  resultEff.push_back(chi2passEE);
  resultEff.push_back(chi2failEE);
  resultEff.push_back(chi2passBE);
  resultEff.push_back(chi2failBE);
  return resultEff;
}

//--------------------------------------------------------------------------------------------------
// perform count for tag and probe efficiency
void performCount(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
	const TString effType, const Bool_t etaRegion, const Int_t iBin, const TString format
){
  Double_t npass=0, ntotal=0;
  npass  = passHist->Integral();
  ntotal = failHist->Integral() + npass;
  resEff  = (ntotal>0) ? npass/ntotal : 0;

  // Calculate the boundaries for the frequentist Clopper-Pearson interval
  resErrl = resEff - TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kFALSE);
  resErrh = TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kTRUE) - resEff;

  // Plot tag and passing- and failing- probe mass distribution
  char binlabelx[100];
  char binlabely[100];
  char effstr[100];
  char ylabel[50];
  char pname[50];
  char yield[50];

  if(!etaRegion) sprintf(binlabelx, "0.0 < |#eta| < 0.9");
  else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");

  sprintf(binlabely, "p_{T} > %i GeV/c",(Int_t)ptCutProbe);
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);
  sprintf(ylabel,"Events / 1 GeV/c^{2}");

  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);

  sprintf(pname,"%s_%s_pass_%d", effType.Data(), etaRegion ? "forward" : "central", iBin);
  sprintf(yield,"%i Events",(UInt_t)npass);
  CPlot plotPass(pname,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.AddHist1D(passHist,"E");
  plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
  plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
  plotPass.AddTextBox(effstr,0.70,0.85,0.95,0.90,0,kBlack,-1);
  plotPass.AddTextBox("#bf{CMS} Work in progress",0.19,0.83,0.54,0.89,0);
  plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  plotPass.Draw(cpass,kTRUE,format);

  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540);
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,0);

  sprintf(pname,"%s_%s_fail_%d", effType.Data(), etaRegion ? "forward" : "central", iBin);
  sprintf(yield,"%i Events",(UInt_t)(ntotal-npass));
  CPlot plotFail(pname,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotFail.AddHist1D(failHist,"E");
  plotFail.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  plotFail.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
  plotFail.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
  plotFail.AddTextBox(effstr,0.70,0.85,0.95,0.90,0,kBlack,-1);
  plotFail.AddTextBox("#bf{CMS} Work in progress",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);

  plotFail.Draw(cfail,kTRUE,format);
}

//--------------------------------------------------------------------------------------------------
// perform fit on one histogram for yield extraxtion
void performFit(
	Double_t &resNsig, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2,
    Double_t &resPurity, Double_t &resPurityErrl, Double_t &resPurityErrh,
    TH1D *h_yield, const Int_t sigMod, const Int_t bkgMod,
	const TString etaRegion, const Bool_t passRegion,
    const Int_t iBin, const TString format, const Int_t fitStrategy
){

    RooRealVar m("m", "mass", massLo, massHi);
    m.setBins(10000);

    CSignalModel     *sigModel = 0;
    CBackgroundModel *bkgModel = 0;

    Int_t nfl=0;

    TFile *histfile = 0;
    TH1D *h=0;
    if(sigMod==2) {
        histfile = new TFile(outputDir+"/histTemplates.root");
        if(passRegion)
            h = (TH1D*)histfile->Get("h_mass_2HLT_"+etaRegion);
        else
            h = (TH1D*)histfile->Get("h_mass_1HLT_"+etaRegion);

      assert(h);
    }

    nfl += set_signal_model(sigMod, sigModel, m, kTRUE, 0, h);
    nfl += set_background_model(bkgMod, bkgModel, m, kTRUE, 0);

    RooAbsData *data = 0;
    data = new RooDataHist("ZReco","ZReco",RooArgList(m),h_yield);

    const Double_t NsigMax = h_yield->Integral();
    RooRealVar Nsig("Nsig","sigYield",NsigMax,0.,1.5*NsigMax);
    RooRealVar Nbkg("Nbkg","bkgYield",0.01*NsigMax,0.,NsigMax);
    RooAddPdf modelPdf("modelTot","Z sig+bkg",RooArgList(*(sigModel->model),*(bkgModel->model)),RooArgList(Nsig,Nbkg));

    RooFitResult *fitResult=0;

    if(fitStrategy == 1){
        // fit bkg shape to sideband region only
        m.setRange("backgroundLow", massLo, 76);
        m.setRange("backgroundHigh", 106, massHi);

        fitResult = bkgModel->model->fitTo(*data,
            RooFit::Range("backgroundLow,backgroundHigh"),
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Strategy(1), // MINOS STRATEGY
            RooFit::Save());
    }
    if(fitStrategy == 2){
        // fit bkg shape to sideband region only
        m.setRange("signalCenter", 81, 101);

        fitResult = sigModel->model->fitTo(*data,
            RooFit::Range("signalCenter"),
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Strategy(1), // MINOS STRATEGY
            RooFit::Save());
    }

    fitResult = modelPdf.fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(1), // MINOS STRATEGY
        //RooFit::Minos(RooArgSet()),
        RooFit::Save());

    // bkgModel->Print();

    fitResult = modelPdf.fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(1), // MINOS STRATEGY
        //RooFit::Minos(RooArgSet()),
        RooFit::Save());

    bkgModel->Print();

    const Int_t nEntries = (Int_t)h_yield->Integral();

    resNsig  = Nsig.getVal();
    resErrl = fabs(Nsig.getErrorLo());
    resErrh = Nsig.getErrorHi();
    resChi2 = 0.;
    resPurity = resNsig/nEntries;
    resPurityErrl = resErrh/resNsig;
    resPurityErrh = resErrl/resNsig;

    char pname[50];
    char strPurity[100];
    sprintf(strPurity,"purity = %.4f_{ -%.4f}^{ +%.4f}",resPurity,resPurityErrl,resPurityErrh);
    resChi2 = make_plot(nfl, m, data, modelPdf, sigModel, bkgModel, fitResult,
        Nsig, Nbkg, strPurity, nEntries, iBin, etaRegion, 0, passRegion);

    delete sigModel;
    delete bkgModel;
    delete data;
}

//--------------------------------------------------------------------------------------------------
// perform fit on two histograms for tag and probe efficiency
void performFit(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2Pass, Double_t &resChi2Fail, TH1D *passHist, TH1D *failHist,
	const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
	const TString bkgQCDTemplate, const TString bkgTTTemplate,
	const TString effType, const Bool_t etaRegion, const Int_t iBin, const TString format, const Int_t fitStrategy
){

  RooRealVar m("m","mass",massLo,massHi);
  m.setBins(10000);

  RooCategory sample("sample","");
  sample.defineType("Pass",1);
  sample.defineType("Fail",2);

  RooAbsData *dataPass=0;
  RooAbsData *dataFail=0;
  RooAbsData *dataCombined=0;

  dataPass     = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
  dataFail     = new RooDataHist("dataFail","dataFail",RooArgSet(m),failHist);
  dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
                                   RooFit::Index(sample),
                                   RooFit::Import("Pass",*((RooDataHist*)dataPass)),
                                   RooFit::Import("Fail",*((RooDataHist*)dataFail)));



  TFile *histfile = 0;
  if(sigpass==2 || sigfail==2) {
    histfile = new TFile(outputDir+"/histTemplates.root");
    assert(histfile);
  }
  std::vector<double> vBkgPars;
  if(bkgfail== 2 or bkgfail==3){
    vBkgPars = preFit(failHist);
  }
  TFile *histbkgQCDfile = 0;
  // TFile *histbkgTTfile = 0;
  if(bkgfail== 7 or bkgpass==7){
    histbkgQCDfile = new TFile(bkgQCDTemplate);
    assert(histbkgQCDfile);
    // if(bkgfail== 7 or bkgpass==7){
    //   histbkgTTfile =  new TFile(bkgTTTemplate);
    //   assert(histbkgTTfile);
    // }
  }
  TH1D *hbkgQCDPass = 0;
  TH1D *hbkgQCDFail = 0;
  // TH1D *hbkgTTPass = 0;
  // TH1D *hbkgTTFail = 0;
  if(bkgpass == 7){
    TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(Form("h_mass_%s_pass_%s", effType.Data(), etaRegion ? "forward" : "central"));
    hbkgQCDPass = hbkg->ProjectionY();
    // smoothening - careful: can introduce bias
    // hbkgQCDPass->GetXaxis()->SetRangeUser(massLo,massHi);
    // hbkgQCDPass->Smooth(10000,"R");
    assert(hbkgQCDPass);
    // if(bkgpass == 7){
    //   hbkgTTPass = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Pass"));
    //   assert(hbkgTTPass);
    // }
  }
  if(bkgfail == 7){
    TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(Form("h_mass_%s_fail_%s", effType.Data(), etaRegion ? "forward" : "central"));
    hbkgQCDFail = hbkg->ProjectionY();
    // smoothening - careful: can introduce bias
    // hbkgQCDFail->GetXaxis()->SetRangeUser(massLo,massHi);
    // hbkgQCDFail->Smooth(10000,"R");
    assert(hbkgQCDFail);
    // if(bkgfail == 7){
    //   hbkgTTFail = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Fail"));
    //   assert(hbkgTTFail);
    // }
  }

  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;

  Int_t nflpass=0, nflfail=0;

  TH1D *h=0;
  if(sigpass==2) {
    h = (TH1D*)histfile->Get(Form("h_mass_pass_%s", etaRegion ? "forward" : "central"));
    h->SetDirectory(0);
  }

  nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 0, h);
  nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 0, hbkgQCDPass);

  if(sigfail==2) {
    TH1D *h = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
    h->SetDirectory(0);
  }

  nflfail += set_signal_model(sigfail, sigFail, m, kFALSE, 0, h);
  nflfail += set_background_model(bkgfail, bkgFail, m, kFALSE, 0, hbkgQCDFail, &vBkgPars);


  Double_t NsigMax     = passHist->Integral()+failHist->Integral();
  Double_t NbkgFailMax = failHist->Integral();
  Double_t NbkgPassMax = passHist->Integral();
  RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0,1.5*NsigMax);
  RooRealVar eff("eff","Efficiency",0.98,0.,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);
  if(bkgpass==0) NbkgPass.setVal(0);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax, 0.0, NbkgFailMax);

  RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
  RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));

  RooAddPdf modelPass("modelPass","Model for PASS sample",
                              (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
                              (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
  RooAddPdf modelFail("modelFail","Model for FAIL sample",RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(modelPass,"Pass");
  totalPdf.addPdf(modelFail,"Fail");

  RooFitResult *fitResult=0;
  Int_t strategy = 2;
  if(effType == "Sel" or effType == "HLT") {
    Nsig.setRange(0,2.0*NsigMax);
    strategy=1;
  }

  RooMsgService::instance().setSilentMode(kTRUE);
  // fit bkg shape in failing probes to sideband region only

  // different fit strategies:
  // Strategy 0: just fit fail and pass region simultaneously in full range
  // Strategy 1: first fit fail pdf in sideband range,
  //    then in complete fail range, then fit fail and pass region simultaneously in full range
  // Strategy 2: Same as 1 but with wider sideband regions

  std::cout<<">>> Fit with strategy "<<fitStrategy<<std::endl;

    if(fitStrategy > 0) {

        if(fitStrategy==1) {
            m.setRange("rangeLow", massLo, 81);
            m.setRange("rangeHigh", 101, massHi);
        }
        else{
            m.setRange("rangeLow", massLo, 76);
            m.setRange("rangeHigh", 106, massHi);
        }

        std::cout<<">>> Fit sideband regions in fail"<<std::endl;

        fitResult = bkgFail->model->fitTo(*dataFail,
                                 RooFit::Range("rangeLow,rangeHigh"),
                                 RooFit::PrintEvalErrors(-1),
                                 RooFit::PrintLevel(-1),
                                 RooFit::Warnings(0),
                                 RooFit::Strategy(strategy), // MINUIT STRATEGY
                                 RooFit::Save());

        std::cout<<">>> Fit full range pdf in fail "<<std::endl;
        // fit total pdf in Fail
        fitResult = modelFail.fitTo(*dataFail,
                                RooFit::PrintEvalErrors(-1),
                                RooFit::PrintLevel(-1),
                                RooFit::Warnings(0),
                                RooFit::Extended(1),
                                RooFit::Strategy(strategy), // MINUIT STRATEGY
                                RooFit::Save());
    }

    fitResult = totalPdf.fitTo(*dataCombined,
                           RooFit::PrintEvalErrors(-1),
                           RooFit::PrintLevel(-1),
                           RooFit::Warnings(0),
                           RooFit::Extended(1),
                           RooFit::Strategy(strategy), // MINUIT STRATEGY
                           RooFit::Minos(RooArgSet(eff)),
                           RooFit::Save());

    fitResult = totalPdf.fitTo(*dataCombined,
                           RooFit::PrintEvalErrors(-1),
                           RooFit::PrintLevel(-1),
                           RooFit::Warnings(0),
                           RooFit::Extended(1),
                           RooFit::Strategy(strategy), // MINUIT STRATEGY
                           RooFit::Minos(RooArgSet(eff)),
                           RooFit::Save());

  bkgFail->Print();
  bkgPass->Print();

  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();

  char effstr[100];
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);

  resChi2Pass = make_plot(nflpass, m, dataPass, modelPass,
      sigPass, bkgPass, fitResult, NsigPass, NbkgPass, effstr, passHist->Integral(), iBin, effType.Data(), etaRegion, kTRUE);

  Double_t yMax = 0.;
  if(effType == "Glo")
    yMax = 1.5*failHist->GetMaximum();

  resChi2Fail = make_plot(nflfail, m, dataFail, modelFail,
      sigFail, bkgFail, fitResult, NsigFail, NbkgFail, effstr, failHist->Integral(), iBin, effType.Data(), etaRegion, kFALSE, yMax);


  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigPass;
  delete bkgPass;
  delete sigFail;
  delete bkgFail;
  delete histfile;
  delete histbkgQCDfile;
  // delete histbkgTTfile;
}

//--------------------------------------------------------------------------------------------------
// perform fit on multiple histograms for simultaneous HLT efficiency and signal yield extraction
// The extraction is based on the following:
//   The number of Z->mumu events where 1 HLT has fired is:
//     N_1 = 2 * epsilon_HLT * (1 - epsilon_HLT) * N_Z + N_Bkg_1
//   The number of Z->mumu events where 2 HLT have fired is:
//     N_2 = epsilon_HLT^2 * N_Z + N_Bkg_1
void performCombinedFit(
	Double_t &resEffB, Double_t &resErrBl, Double_t &resErrBh,
    Double_t &resEffE, Double_t &resErrEl, Double_t &resErrEh,
    Double_t &resYieldBB, Double_t &resYieldBBl, Double_t &resYieldBBh,
    Double_t &resYieldEE, Double_t &resYieldEEl, Double_t &resYieldEEh,
    Double_t &resYieldBE, Double_t &resYieldBEl, Double_t &resYieldBEh,
    Double_t &resChi2PassBB, Double_t &resChi2FailBB,
    Double_t &resChi2PassEE, Double_t &resChi2FailEE,
    Double_t &resChi2PassBE, Double_t &resChi2FailBE,
    TH1D *passHistBB, TH1D *failHistBB,
    TH1D *passHistEE, TH1D *failHistEE,
    TH1D *passHistBE, TH1D *failHistBE,
    Double_t cBB_, Double_t cBE_, Double_t cEE_,
	const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
	const TString bkgQCDTemplate, const TString bkgTTTemplate,
    const Int_t iBin,
    const TString format,
    const Int_t fitStrategy
){

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("PassBB",1);
    sample.defineType("FailBB",2);
    sample.defineType("PassEE",3);
    sample.defineType("FailEE",4);
    sample.defineType("PassBE",5);
    sample.defineType("FailBE",6);

    RooAbsData *dataPassBB = new RooDataHist("dataPassBB","dataPassBB",RooArgSet(m),passHistBB);
    RooAbsData *dataFailBB = new RooDataHist("dataFailBB","dataFailBB",RooArgSet(m),failHistBB);
    RooAbsData *dataPassEE = new RooDataHist("dataPassEE","dataPassEE",RooArgSet(m),passHistEE);
    RooAbsData *dataFailEE = new RooDataHist("dataFailEE","dataFailEE",RooArgSet(m),failHistEE);
    RooAbsData *dataPassBE = new RooDataHist("dataPassBE","dataPassBE",RooArgSet(m),passHistBE);
    RooAbsData *dataFailBE = new RooDataHist("dataFailBE","dataFailBE",RooArgSet(m),failHistBE);

    RooAbsData *dataCombinedBB = new RooDataHist(
        "dataCombinedBB","dataCombinedBB",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("PassBB",*((RooDataHist*)dataPassBB)),
        RooFit::Import("FailBB",*((RooDataHist*)dataFailBB)));

    RooAbsData *dataCombinedEE = new RooDataHist(
        "dataCombinedEE","dataCombinedEE",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("PassEE",*((RooDataHist*)dataPassEE)),
        RooFit::Import("FailEE",*((RooDataHist*)dataFailEE)));

    RooAbsData *dataCombinedBE = new RooDataHist(
        "dataCombinedBE","dataCombinedBE",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("PassBE",*((RooDataHist*)dataPassBE)),
        RooFit::Import("FailBE",*((RooDataHist*)dataFailBE)));

    RooAbsData *dataCombined = new RooDataHist(
        "dataCombined", "dataCombined", RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("PassBB", *((RooDataHist*)dataPassBB)),
        RooFit::Import("FailBB", *((RooDataHist*)dataFailBB)),
        RooFit::Import("PassEE", *((RooDataHist*)dataPassEE)),
        RooFit::Import("FailEE", *((RooDataHist*)dataFailEE)),
        RooFit::Import("PassBE", *((RooDataHist*)dataPassBE)),
        RooFit::Import("FailBE", *((RooDataHist*)dataFailBE))
    );

    TFile *histfile = 0;
    if(sigpass==2 || sigfail==2) {
        histfile = new TFile(outputDir+"/histTemplates.root");
        assert(histfile);
    }
    // std::vector<double> vBkgPars;
    // if(bkgfail== 2 or bkgfail==3){
    //   vBkgPars = preFit(failHist);
    // }
    // TFile *histbkgQCDfile = 0;
    // if(bkgfail== 7 or bkgpass==7){
    //   histbkgQCDfile = new TFile(bkgQCDTemplate);
    //   assert(histbkgQCDfile);
    // }
    TH1D *hbkgQCDPass = 0;
    TH1D *hbkgQCDFail = 0;
    // if(bkgpass == 7){
    //   TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(Form("h_mass_%s_pass_%s", effType.Data(), etaRegion ? "forward" : "central"));
    //   hbkgQCDPass = hbkg->ProjectionY();
    //   assert(hbkgQCDPass);
    // }
    // if(bkgfail == 7){
    //   TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(Form("h_mass_%s_fail_%s", effType.Data(), etaRegion ? "forward" : "central"));
    //   hbkgQCDFail = hbkg->ProjectionY();
    //   assert(hbkgQCDFail);
    // }

    CSignalModel     *sigPassBB = 0;
    CBackgroundModel *bkgPassBB = 0;
    CSignalModel     *sigFailBB = 0;
    CBackgroundModel *bkgFailBB = 0;

    CSignalModel     *sigPassEE = 0;
    CBackgroundModel *bkgPassEE = 0;
    CSignalModel     *sigFailEE = 0;
    CBackgroundModel *bkgFailEE = 0;

    CSignalModel     *sigPassBE = 0;
    CBackgroundModel *bkgPassBE = 0;
    CSignalModel     *sigFailBE = 0;
    CBackgroundModel *bkgFailBE = 0;

    Int_t nflpassBB=0, nflfailBB=0;
    Int_t nflpassEE=0, nflfailEE=0;
    Int_t nflpassBE=0, nflfailBE=0;

    TH1D *hPassBB=0;
    TH1D *hPassEE=0;
    TH1D *hPassBE=0;
    if(sigpass==2) {
        hPassBB = (TH1D*)histfile->Get("h_mass_pass_BB");
        hPassEE = (TH1D*)histfile->Get("h_mass_pass_EE");
        hPassBE = (TH1D*)histfile->Get("h_mass_pass_BE");
        hPassBB->SetDirectory(0);
        hPassEE->SetDirectory(0);
        hPassBE->SetDirectory(0);
    }

    nflpassBB += set_signal_model(sigpass, sigPassBB, m, kTRUE, 0, hPassBB);
    nflpassBB += set_background_model(bkgpass, bkgPassBB, m, kTRUE, 0, hbkgQCDPass);
    nflpassEE += set_signal_model(sigpass, sigPassEE, m, kTRUE, 1, hPassEE);
    nflpassEE += set_background_model(bkgpass, bkgPassEE, m, kTRUE, 1, hbkgQCDPass);
    nflpassBE += set_signal_model(sigpass, sigPassBE, m, kTRUE, 2, hPassBE);
    nflpassBE += set_background_model(bkgpass, bkgPassBE, m, kTRUE, 2, hbkgQCDPass);

    TH1D *hFailBB=0;
    TH1D *hFailEE=0;
    TH1D *hFailBE=0;
    if(sigfail==2) {
        hFailBB = (TH1D*)histfile->Get("h_mass_fail_BB");
        hFailEE = (TH1D*)histfile->Get("h_mass_fail_EE");
        hFailBE = (TH1D*)histfile->Get("h_mass_fail_BE");
        hFailBB->SetDirectory(0);
        hFailEE->SetDirectory(0);
        hFailBE->SetDirectory(0);
    }

    nflfailBB += set_signal_model(sigfail, sigFailBB, m, kFALSE, 0, hFailBB);
    nflfailBB += set_background_model(bkgfail, bkgFailBB, m, kFALSE, 0, hbkgQCDFail);
    nflfailEE += set_signal_model(sigfail, sigFailEE, m, kFALSE, 1, hFailEE);
    nflfailEE += set_background_model(bkgfail, bkgFailEE, m, kFALSE, 1, hbkgQCDFail);
    nflfailBE += set_signal_model(sigfail, sigFailBE, m, kFALSE, 2, hFailBE);
    nflfailBE += set_background_model(bkgfail, bkgFailBE, m, kFALSE, 2, hbkgQCDFail);

    Double_t NsigMaxBB     = passHistBB->Integral() + failHistBB->Integral();
    Double_t NbkgFailMaxBB = failHistBB->Integral();
    Double_t NbkgPassMaxBB = passHistBB->Integral();
    Double_t NsigMaxEE     = passHistEE->Integral() + failHistEE->Integral();
    Double_t NbkgFailMaxEE = failHistEE->Integral();
    Double_t NbkgPassMaxEE = passHistEE->Integral();
    Double_t NsigMaxBE     = passHistBE->Integral() + failHistBE->Integral();
    Double_t NbkgFailMaxBE = failHistBE->Integral();
    Double_t NbkgPassMaxBE = passHistBE->Integral();

    RooRealVar effB("effB","Efficiency barrel",0.95,0.,1.0);
    RooRealVar effE("effE","Efficiency endcap",0.95,0.,1.0);

    RooRealVar NsigBB("NsigBB","Signal Yield BB",NsigMaxBB,0,1.5*NsigMaxBB);
    RooRealVar NsigEE("NsigEE","Signal Yield EE",NsigMaxEE,0,1.5*NsigMaxEE);
    RooRealVar NsigBE("NsigBE","Signal Yield BE",NsigMaxBE,0,1.5*NsigMaxBE);

    RooRealVar cBB("cBB", "Correlation factor BB", cBB_);
    RooRealVar cBE("cBE", "Correlation factor BE", cBE_);
    RooRealVar cEE("cEE", "Correlation factor EE", cEE_);

    RooFormulaVar NsigPassBB("NsigPassBB","effB*effB*NsigBB*cBB",RooArgList(effB,NsigBB,cBB));
    RooFormulaVar NsigFailBB("NsigFailBB","2*effB*(1.0-cBB*effB)*NsigBB",RooArgList(effB,NsigBB,cBB));
    RooFormulaVar NsigPassEE("NsigPassEE","effE*effE*NsigEE*cEE",RooArgList(effE,NsigEE,cEE));
    RooFormulaVar NsigFailEE("NsigFailEE","2*effE*(1.0-cEE*effE)*NsigEE",RooArgList(effE,NsigEE,cEE));
    RooFormulaVar NsigPassBE("NsigPassBE","effB*effE*NsigBE*cBE",RooArgList(effB,effE,NsigBE,cBE));
    RooFormulaVar NsigFailBE("NsigFailBE","( effB + effE -2*effB*effE*cBE )*NsigBE",RooArgList(effB,effE,NsigBE,cBE));

    RooRealVar NbkgPassBB("NbkgPassBB","Background count in BB PASS sample",0.01*NbkgPassMaxBB, 0.0, NbkgPassMaxBB);
    RooRealVar NbkgFailBB("NbkgFailBB","Background count in BB FAIL sample",0.05*NbkgFailMaxBB, 0.0, NbkgFailMaxBB);
    RooRealVar NbkgPassEE("NbkgPassEE","Background count in EE PASS sample",0.01*NbkgPassMaxEE, 0.0, NbkgPassMaxEE);
    RooRealVar NbkgFailEE("NbkgFailEE","Background count in EE FAIL sample",0.05*NbkgFailMaxEE, 0.0, NbkgFailMaxEE);
    RooRealVar NbkgPassBE("NbkgPassBE","Background count in BE PASS sample",0.01*NbkgPassMaxBE, 0.0, NbkgPassMaxBE);
    RooRealVar NbkgFailBE("NbkgFailBE","Background count in BE FAIL sample",0.05*NbkgFailMaxBE, 0.0, NbkgFailMaxBE);

    RooAddPdf modelPassBB("modelPassBB ","Model for BB PASS sample",
        RooArgList(*(sigPassBB->model),*(bkgPassBB->model)),
        RooArgList(NsigPassBB, NbkgPassBB));
    RooAddPdf modelFailBB("modelFailBB ","Model for BB FAIL sample",
        RooArgList(*(sigFailBB->model),*(bkgFailBB->model)),
        RooArgList(NsigFailBB, NbkgFailBB));

    RooAddPdf modelPassEE("modelPassEE ","Model for EE PASS sample",
        RooArgList(*(sigPassEE->model),*(bkgPassEE->model)),
        RooArgList(NsigPassEE, NbkgPassEE));
    RooAddPdf modelFailEE("modelFailEE ","Model for EE FAIL sample",
        RooArgList(*(sigFailEE->model),*(bkgFailEE->model)),
        RooArgList(NsigFailEE, NbkgFailEE));

    RooAddPdf modelPassBE("modelPassBE ","Model for BE PASS sample",
        RooArgList(*(sigPassBE->model),*(bkgPassBE->model)),
        RooArgList(NsigPassBE, NbkgPassBE));
    RooAddPdf modelFailBE("modelFailBE ","Model for BE FAIL sample",
        RooArgList(*(sigFailBE->model),*(bkgFailBE->model)),
        RooArgList(NsigFailBE, NbkgFailBE));

    RooSimultaneous pdfBB("pdfBB","pdfBB",sample);
    pdfBB.addPdf(modelPassBB,"PassBB");
    pdfBB.addPdf(modelFailBB,"FailBB");

    RooSimultaneous pdfEE("pdfEE","pdfEE",sample);
    pdfEE.addPdf(modelPassEE,"PassEE");
    pdfEE.addPdf(modelFailEE,"FailEE");

    RooSimultaneous pdfBE("pdfBE","pdfBE",sample);
    pdfBE.addPdf(modelPassBE,"PassBE");
    pdfBE.addPdf(modelFailBE,"FailBE");

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPassBB,"PassBB");
    totalPdf.addPdf(modelFailBB,"FailBB");
    totalPdf.addPdf(modelPassEE,"PassEE");
    totalPdf.addPdf(modelFailEE,"FailEE");
    totalPdf.addPdf(modelPassBE,"PassBE");
    totalPdf.addPdf(modelFailBE,"FailBE");

    RooFitResult *fitResult=0, *fitResultBB=0, *fitResultEE=0, *fitResultBE=0;
    Int_t strategy = 2;

    RooMsgService::instance().setSilentMode(kTRUE);

    m.setRange("rangeLow", massLo, 76);
    m.setRange("rangeHigh", 106, massHi);

    // fit BB region only
    std::cout<<"--- fit BB region only -- "<<std::endl;
    // fitResult = bkgFailBB->model->fitTo(*dataFailBB,
    //                          RooFit::Range("rangeLow,rangeHigh"),
    //                          RooFit::PrintEvalErrors(-1),
    //                          RooFit::PrintLevel(-1),
    //                          RooFit::Warnings(0),
    //                          RooFit::Strategy(strategy), // MINUIT STRATEGY
    //                          RooFit::Save());

    fitResultBB = pdfBB.fitTo(*dataCombinedBB,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        RooFit::Minos(RooArgSet(effB)),
        RooFit::Save());

    // fit EE region only
    std::cout<<"--- fit EE region only -- "<<std::endl;
    // fitResult = bkgFailEE->model->fitTo(*dataFailEE,
    //                          RooFit::Range("rangeLow,rangeHigh"),
    //                          RooFit::PrintEvalErrors(-1),
    //                          RooFit::PrintLevel(-1),
    //                          RooFit::Warnings(0),
    //                          RooFit::Strategy(strategy), // MINUIT STRATEGY
    //                          RooFit::Save());

    fitResultEE = pdfEE.fitTo(*dataCombinedEE,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        RooFit::Minos(RooArgSet(effE)),
        RooFit::Save());

    // fit BE region only
    std::cout<<"--- fit BE region only -- "<<std::endl;
    // fitResult = bkgFailBE->model->fitTo(*dataFailBE,
    //                          RooFit::Range("rangeLow,rangeHigh"),
    //                          RooFit::PrintEvalErrors(-1),
    //                          RooFit::PrintLevel(-1),
    //                          RooFit::Warnings(0),
    //                          RooFit::Strategy(strategy), // MINUIT STRATEGY
    //                          RooFit::Save());

    fitResultBE = pdfBE.fitTo(*dataCombinedBE,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        RooFit::Minos(RooArgSet(effE, effB)),
        RooFit::Save());

    // fit all regions together
    std::cout<<"--- fit all regions together -- "<<std::endl;
    fitResult = totalPdf.fitTo(*dataCombined,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        // RooFit::Minos(RooArgSet(effB, effE)),
        RooFit::Save());

    fitResult = totalPdf.fitTo(*dataCombined,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        RooFit::Minos(RooArgSet(NsigBB, NsigEE, NsigBE)),
        RooFit::Save());

    bkgFailBB->Print();
    bkgPassBB->Print();
    bkgFailEE->Print();
    bkgPassEE->Print();
    bkgFailBE->Print();
    bkgPassBE->Print();

    resEffB  = effB.getVal();
    resErrBl = fabs(effB.getErrorLo());
    resErrBh = effB.getErrorHi();
    resEffE  = effE.getVal();
    resErrEl = fabs(effE.getErrorLo());
    resErrEh = effE.getErrorHi();

    resYieldBB  = NsigBB.getVal();
    resYieldBBl = fabs(NsigBB.getErrorLo());
    resYieldBBh = NsigBB.getErrorHi();
    resYieldEE  = NsigEE.getVal();
    resYieldEEl = fabs(NsigEE.getErrorLo());
    resYieldEEh = NsigEE.getErrorHi();
    resYieldBE  = NsigBE.getVal();
    resYieldBEl = fabs(NsigBE.getErrorLo());
    resYieldBEh = NsigBE.getErrorHi();

    char effstr[100];
    sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEffB,resErrBl,resErrBh);

    resChi2PassBB = make_plot(nflpassBB, m, dataPassBB, modelPassBB,
        sigPassBB, bkgPassBB, fitResult, NsigPassBB, NbkgPassBB, effstr, passHistBB->Integral(), iBin, "BB", 0, kTRUE);

    resChi2FailBB = make_plot(nflfailBB, m, dataFailBB, modelFailBB,
        sigFailBB, bkgFailBB, fitResult, NsigFailBB, NbkgFailBB, effstr, failHistBB->Integral(), iBin, "BB", 0, kFALSE);

    sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEffE,resErrEl,resErrEh);

    resChi2PassEE = make_plot(nflpassEE, m, dataPassEE, modelPassEE,
        sigPassEE, bkgPassEE, fitResult, NsigPassEE, NbkgPassEE, effstr, passHistEE->Integral(), iBin, "EE", 1, kTRUE);

    resChi2FailEE = make_plot(nflfailEE, m, dataFailEE, modelFailEE,
        sigFailEE, bkgFailEE, fitResult, NsigFailEE, NbkgFailEE, effstr, failHistEE->Integral(), iBin, "EE", 1, kFALSE);

    sprintf(effstr,"#varepsilon_{B} = %.3f_{ -%.3f}^{ +%.3f} | #varepsilon_{E} = %.3f_{ -%.3f}^{ +%.3f}",
        resEffB,resErrBl,resErrBh,resEffE,resErrEl,resErrEh);

    resChi2PassBE = make_plot(nflpassBE, m, dataPassBE, modelPassBE,
        sigPassBE, bkgPassBE, fitResult, NsigPassBE, NbkgPassBE, effstr, passHistBE->Integral(), iBin, "BE", 1, kTRUE);

    resChi2FailBE = make_plot(nflfailBE, m, dataFailBE, modelFailBE,
        sigFailBE, bkgFailBE, fitResult, NsigFailBE, NbkgFailBE, effstr, failHistBE->Integral(), iBin, "BE", 1, kFALSE);


    delete dataCombined;
    delete dataPassBB;
    delete dataFailBB;
    delete dataPassEE;
    delete dataFailEE;
    delete dataPassBE;
    delete dataFailBE;
    delete sigPassBB;
    delete bkgPassBB;
    delete sigFailBB;
    delete bkgFailBB;
    delete sigPassEE;
    delete bkgPassEE;
    delete sigFailEE;
    delete bkgFailEE;
    delete sigPassBE;
    delete bkgPassBE;
    delete sigFailBE;
    delete bkgFailBE;
    delete histfile;
    // delete histbkgQCDfile;
}

//--------------------------------------------------------------------------------------------------
std::vector<double> preFit(TH1D* failHist){
    std::cout<<"PREFIT"<<std::endl;
    // do not fit between 81 to 101 GeV

//  std::vector<float> v = {1.,0.,1.,0.,1.,0.};return v;
  TH1D *h = new TH1D("h", "", massBin, massLo, massHi);
  TF1 *fq = new TF1("fq", "[0]+[1]*x+[2]*x*x", massLo, massHi);
  TF1 *fl = new TF1("fl", "[0]+[1]*x", massLo, massHi);

  const int binwidth = (massHi-massLo)/massBin;

  for(int i = 0; i < (int)(81-massLo)/binwidth; i++){
    h->SetBinContent(i+1, failHist->GetBinContent(i+1));
    h->SetBinError(i+1, failHist->GetBinError(i+1));
  }
  for(int i = (int)(massHi-101)/binwidth; i < massBin; i++){
    h->SetBinContent(i+1, failHist->GetBinContent(i+1));
    h->SetBinError(i+1, failHist->GetBinError(i+1));
  }

  h->Fit("fq");
  std::cout<<fq->GetParameter(0)<<","<<fq->GetParError(0)<<" ,"<<fq->GetParameter(1)<<","<<fq->GetParError(1)<<","<<fq->GetParameter(2)<<","<<fq->GetParError(2)<<std::endl;
  std::vector<double> v;

  if(fq->GetParameter(2) > 0.){
    h->Fit("fl");
    std::cout<<"SWITCH TO LINEAR!"<<std::endl;
    std::cout<<fl->GetParameter(0)<<","<<fl->GetParError(0)<<" ,"<<fl->GetParameter(1)<<","<<fl->GetParError(1)<<std::endl;
    v = {fl->GetParameter(0), fl->GetParError(0), fl->GetParameter(1), fl->GetParError(1), 0.,0.};
  }else{
    v = {fq->GetParameter(0), fq->GetParError(0), fq->GetParameter(1), fq->GetParError(1), fq->GetParameter(2), fq->GetParError(2)};
  }

  return v;
}
