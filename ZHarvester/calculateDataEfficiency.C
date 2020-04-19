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


// Set up fiducial region
const Float_t massLo  = 56.;
const Float_t massHi  = 116.;
const UInt_t  massBin = (UInt_t)(massHi-massLo);

const Float_t etaCutTag   = 2.4;
const Float_t etaCutProbe = 2.4;
const Float_t etaBound    = 0.9;

// Set up pile-up bounds available get MC template for fitting
const Int_t minPU = 1;
const Int_t maxPU = 60;

void generateTemplate(
        const TString mcfilename,
        const Float_t ptCutTag, const Float_t ptCutProbe, TH1D *hPV=0, const TString outputDir="");

void generateTemplate_ZYield(
	const TString mcfilename,
	const Float_t ptCutTag,
	const Float_t ptCutProbe,
	TH1D          *hPV,
    const TString outputDir
);

void performCount(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
        const Float_t ptCutTag, const Float_t ptCutProbe,
        const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

void performFit(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2Pass, Double_t &resChi2Fail, TH1D *passHist, TH1D *failHist,
        const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
        const Float_t ptCutTag, const Float_t ptCutProbe, const TString outputDir, const TString bkgQCDTemplate, const TString bkgTTTemplate,
        const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

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
    TH1D *hist=0
){

    switch(model_type) {
        case 1:
            model = new CBreitWignerConvCrystalBall(param_mass, pass, 0);
            return 4;
        case 2:
            model = new CMCTemplateConvGaussian(param_mass, hist, pass, 0);
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
    TH1D *hist=0,
    std::vector<double> *params=0
){
    switch(model_type) {
        case 1:
            model = new CExponential(param_mass, pass, 0);
            return 1;
        case 2:
            if(params != 0)
                model = new CQuadratic(param_mass, pass, 0,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadratic(param_mass, pass, 0, 0.,0.,0.,0.,0.,0.);
            return 3;
        case 3:
            if(params != 0)
                model = new CQuadPlusExp(param_mass, pass, 0,
                    params->at(0), params->at(1), params->at(2), params->at(3), params->at(4), params->at(5));
            else
                model = new CQuadPlusExp(param_mass, pass, 0, 0.,0.,0.,0.,0.,0.);
                return 4;
        case 4:
            model = new CDas(param_mass, pass, 0);
            return 4;
        case 5:
            model = new CDasPlusExp(param_mass, pass, 0);
            return 6;
        case 6:
            model = new CQCD(param_mass, hist, pass, 0);
            return 1;
    }
    return 0;
}

//--------------------------------------------------------------------------------------------------
void generateTemplate(
	const TString mcfilename,
	const Float_t ptCutTag,
	const Float_t ptCutProbe,
	TH1D          *hPV,
    const TString outputDir
){
  cout << "Creating histogram templates... "; cout.flush();

  TFile *infile    = new TFile(mcfilename);
  TTree *eventTree = (TTree*)infile->Get("tree");
  TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

  if(hPV)
    hPV->Divide(hPVtemplate);

  Float_t mass, ptTag, etaTag, ptProbe, etaProbe;
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
	const Float_t ptCutTag,
	const Float_t ptCutProbe,
	TH1D          *hPV,
    const TString outputDir
){
  cout << "Creating histogram templates... "; cout.flush();

  TFile *infile    = new TFile(mcfilename);
  TTree *eventTree = (TTree*)infile->Get("tree");
  TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

  if(hPV)
    hPV->Divide(hPVtemplate);

  Float_t mass, ptTag, etaTag, ptProbe, etaProbe;
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

  TH1D *h_mass_zyield = new TH1D("h_mass_zyield", "", massBin, massLo, massHi);


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

    h_mass_zyield->Fill(mass, wgt);

  }

  TFile outfile = TFile(outputDir+"/histTemplates.root", "RECREATE");
  h_mass_zyield->Write();
  outfile.Write();
  outfile.Close();

  infile->Close();
  delete infile;

  cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
template<typename T>
Double_t make_plot(
    const Float_t ptCutTag,
    const Float_t lumi,
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
    const TString format="png"
){

    char pname[50];
    char ctitle[100];
    char binlabelx[100];
    char binlabely[100];
    char lumitext[100];
    char ylabel[50];
    char yield[50];
    char nsigstr[100];
    char nbkgstr[100];
    char chi2str[100];

    if(effType == "yield"){
        sprintf(pname,"%s_%d", effType.Data(), iBin);
        sprintf(ctitle,"%s %d ", effType.Data(), iBin);
        sprintf(binlabelx, "0.0 < |#eta| < 2.4");
    }
    else{
        sprintf(pname,"%s_%s_%s_%d", effType.Data(), etaRegion ? "forward" : "central", passRegion ? "pass" : "fail", iBin);
        sprintf(ctitle,"%s %s %s %d ", effType.Data(), etaRegion ? "forward" : "central", passRegion ? "pass" : "fail", iBin);

        if(!etaRegion) sprintf(binlabelx, "0.0 < |#eta| < 0.9");
        else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");
    }

    TCanvas *cyield = MakeCanvas(pname, ctitle,720,540);
    cyield->SetWindowPosition(cyield->GetWindowTopX()+cyield->GetBorderSize()+800,0);

    sprintf(binlabely, "%i GeV/c < p_{T} < 13000 GeV/c",(Int_t)ptCutTag);
    sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
    sprintf(ylabel,"Events / 1 GeV/c^{2}");

    RooPlot *mframe = param_mass.frame(Bins(massBin));
    data->plotOn(mframe,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelPdf.plotOn(mframe, Components(*(bkgModel->model)), LineColor(8));
    modelPdf.plotOn(mframe, Components(*(sigModel->model)), LineColor(9));
    modelPdf.plotOn(mframe, LineColor(kRed));

    double a = Nsig.getVal(), aErr = Nsig.getPropagatedError(*fitResult);
    double b = Nbkg.getVal(), bErr = Nbkg.getPropagatedError(*fitResult);

    sprintf(yield,"%u Events",nEntries);
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",a,aErr);
    Double_t resChi2 = mframe->chiSquare(nfl);
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",resChi2);
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",b,bErr);

    CPlot plotPass(pname, mframe, ctitle, "tag-probe mass [GeV/c^{2}]",ylabel);
    plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
    plotPass.AddTextBox(sigModel->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotPass.AddTextBox(bkgModel->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotPass.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);

    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
    plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);

    plotPass.Draw(cyield,kTRUE, format);
    delete cyield;

    return resChi2;
}

//--------------------------------------------------------------------------------------------------
std::vector<float> getZyield(
                TH1D *h_yield,
                const TString outputDir,    // output directory
                const Int_t   iBin,         // Label of measurement in currect run
                const Int_t   sigMod,
                const Int_t   bkgMod,
                const Float_t ptCutTag,
                const Float_t ptCutProbe,
                const Float_t lumi=10.,     // luminosity for plot label
                const TString  mcfilename="",
                TH1D* _hQCD = 0,
                TH1D* hPV = 0
){

  CPlot::sOutDir = outputDir;

  Double_t NsigMax = h_yield->GetEntries();
  if(sigMod == 0){      // perform count - expect 1% fakes
    std::vector<float> resultEff = {};
    resultEff.push_back(NsigMax*0.99);
    resultEff.push_back(std::sqrt(NsigMax)*0.99);
    resultEff.push_back(std::sqrt(NsigMax)*0.99);
    resultEff.push_back(0.);
    resultEff.push_back(0.01);
    return resultEff;
  }

  RooRealVar m("m","mass",massLo,massHi);
  m.setBins(10000);

  CSignalModel     *sigModel = 0;
  CBackgroundModel *bkgModel = 0;

  Int_t nfl=0;

  TH1D *h=0;
  if(sigMod==2) {
    TFile *histfile = 0;
    generateTemplate_ZYield(mcfilename, ptCutTag, ptCutProbe, hPV, outputDir);
    histfile = new TFile(outputDir+"/histTemplates.root");
    h = (TH1D*)histfile->Get("h_mass_zyield");
    h->SetDirectory(0);
  }

  nfl += set_signal_model(sigMod, sigModel, m, kTRUE, h);
  nfl += set_background_model(bkgMod, bkgModel, m, kTRUE, _hQCD);

  RooAbsData *data = 0;
  data = new RooDataHist("ZReco","ZReco",RooArgList(m),h_yield);

  RooRealVar Nsig("Nsig","sigYield",NsigMax,0.,1.5*NsigMax);
  RooRealVar Nbkg("Nbkg","bkgYield",0.01*NsigMax,0.,NsigMax);
  RooAddPdf modelPdf("model","Z sig+bkg",RooArgList(*(sigModel->model),*(bkgModel->model)),RooArgList(Nsig,Nbkg));

  // fit bkg shape in failing probes to sideband region only
  m.setRange("backgroundLow", massLo, 76);
  m.setRange("backgroundHigh", 106, massHi);

  RooFitResult *fitResult=0;
  fitResult = bkgModel->model->fitTo(*data,
                            RooFit::Range("backgroundLow,backgroundHigh"),
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Strategy(1), // MINOS STRATEGY
                             RooFit::Save());

  fitResult = modelPdf.fitTo(*data,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(1), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet()),
                             RooFit::Save());

  const Int_t nEntries = (Int_t)h_yield->GetEntries();

  Double_t resNsig  = Nsig.getVal();
  Double_t resErrl = fabs(Nsig.getErrorLo());
  Double_t resErrh = Nsig.getErrorHi();
  Double_t resChi2 = 0.;
  Double_t resPurity = resNsig/nEntries;
  Double_t resPurityErrl = resErrh/resNsig;
  Double_t resPurityErrh = resErrl/resNsig;

  char pname[50];
  char strPurity[100];
  sprintf(strPurity,"purity = %.4f_{ -%.4f}^{ +%.4f}",resPurity,resPurityErrl,resPurityErrh);
  resChi2 = make_plot(ptCutTag, lumi, nfl, m, data, modelPdf, sigModel, bkgModel, fitResult, Nsig, Nbkg, strPurity, nEntries, iBin);

  delete sigModel;
  delete bkgModel;
  delete data;

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
		const TString outputDir,            // output directory
		const Int_t   iBin,                 // Label of measurement in currect run
		const TString effType,              // "HLT" or "SIT" or "Glo" or "Sta" or "Trk"
        const Bool_t  etaRegion,
		const Int_t   sigModPass,           // signal extraction method for PASS sample
		const Int_t   bkgModPass,           // background model for PASS sample
		const Int_t   sigModFail,           // signal extraction method for FAIL sample
		const Int_t   bkgModFail,           // background model for FAIL sample
		const Float_t ptCutTag,             // for generating template and printing cuts on plots
		const Float_t ptCutProbe,           //
        TH1D          *hPV=0,
		const Float_t lumi=10.,             // luminosity for plot label
        const TString mcfilename="",        // ROOT file containing MC events to generate templates from
        const TString bkgQCDFilename="",    // ROOT file containing bkg template
        const TString bkgTTFilename="",
		const TString format="png"          // plot format
){

  CPlot::sOutDir = outputDir;

  // Generate histogram templates from MC if necessary
  if(sigModPass==2 || sigModFail==2) {
    generateTemplate(mcfilename, ptCutTag, ptCutProbe, hPV, outputDir);
  }

  Double_t eff  = 0.;
  Double_t errl = 0.;
  Double_t errh = 0.;
  Double_t chi2pass = 999.;
  Double_t chi2fail = 999.;

  if(sigModPass == 0){
    performCount(eff, errl, errh, h_mass_pass, h_mass_fail, ptCutTag, ptCutProbe,
      effType, etaRegion, iBin, lumi, format);
  }else{
    performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail,
      sigModPass, bkgModPass, sigModFail, bkgModFail,
      ptCutTag, ptCutProbe, outputDir, bkgQCDFilename, bkgTTFilename,
      effType, etaRegion, iBin, lumi, format);
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
void performCount(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
	const Float_t ptCutTag, const Float_t ptCutProbe,
	const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format
){
  Double_t npass=0, ntotal=0;
  npass  = passHist->GetEntries();
  ntotal = failHist->GetEntries() + npass;
  resEff  = (ntotal>0) ? npass/ntotal : 0;

  // Calculate the boundaries for the frequentist Clopper-Pearson interval
  resErrl = resEff - TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kFALSE);
  resErrh = TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kTRUE) - resEff;

  // Plot tag and passing- and failing- probe mass distribution
  char binlabelx[100];
  char binlabely[100];
  char effstr[100];
  char lumitext[100];
  char ylabel[50];
  char pname[50];
  char yield[50];

  if(!etaRegion) sprintf(binlabelx, "0.0 < |#eta| < 0.9");
  else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");

  sprintf(binlabely, "%i GeV/c < p_{T} < 13000 GeV/c",(Int_t)ptCutProbe);
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
//  sprintf(lumitext,"  %.0f lumi-sections at  #sqrt{s} = 13 TeV",lumi);
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
  plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
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
  plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  plotFail.Draw(cfail,kTRUE,format);
}

//--------------------------------------------------------------------------------------------------
void performFit(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2Pass, Double_t &resChi2Fail, TH1D *passHist, TH1D *failHist,
	const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
	const Float_t ptCutTag, const Float_t ptCutProbe, const TString outputDir, const TString bkgQCDTemplate, const TString bkgTTTemplate,
	const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format
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
  TFile *histbkgTTfile = 0;
  if(bkgfail== 6 or bkgpass==6 or bkgfail== 7 or bkgpass==7){
    histbkgQCDfile =  new TFile(bkgQCDTemplate);
    assert(histbkgQCDfile);
    if(bkgfail== 7 or bkgpass==7){
      histbkgTTfile =  new TFile(bkgTTTemplate);
      assert(histbkgTTfile);
    }
  }
  TH1D *hbkgQCDPass = 0;
  TH1D *hbkgQCDFail = 0;
  TH1D *hbkgTTPass = 0;
  TH1D *hbkgTTFail = 0;
  if(bkgpass == 6 or bkgpass == 7){
    hbkgQCDPass = (TH1D*)histbkgQCDfile->Get(Form("bkg_template_%s%s", effType.Data(), "Pass"));
    assert(hbkgQCDPass);
    if(bkgpass == 7){
      hbkgTTPass = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Pass"));
      assert(hbkgTTPass);
    }
  }
  if(bkgfail == 6 or bkgfail == 7){
    hbkgQCDFail = (TH1D*)histbkgQCDfile->Get(Form("bkg_template_%s%s", effType.Data(), "Fail"));
    assert(hbkgQCDFail);
    if(bkgfail == 7){
      hbkgTTFail = (TH1D*)histbkgTTfile->Get(Form("bkg_template_%s%s", effType.Data(), "Fail"));
      assert(hbkgTTFail);
    }
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

  nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, h);
  nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE);

  if(sigfail==2) {
    TH1D *h = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
    h->SetDirectory(0);
  }

  nflfail += set_signal_model(sigfail, sigFail, m, kFALSE, h);
  nflfail += set_background_model(bkgfail, bkgFail, m, kFALSE, 0, &vBkgPars);


  Double_t NsigMax     = passHist->Integral()+failHist->Integral();
  Double_t NbkgFailMax = failHist->Integral();
  Double_t NbkgPassMax = passHist->Integral();
  RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0,1.5*NsigMax);
  RooRealVar eff("eff","Efficiency",1.0,0.0,1.0);
  RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.1*NbkgPassMax,0.01,NbkgPassMax);
  if(bkgpass==0) NbkgPass.setVal(0);
  RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax,0.01,NbkgFailMax);

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
  m.setRange("backgroundLow", massLo, 76);
  m.setRange("backgroundHigh", 106, massHi);

  fitResult = bkgFail->model->fitTo(*dataFail,
                            RooFit::Range("backgroundLow,backgroundHigh"),
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             RooFit::Save());

  // fit total pdf
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(1),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(1),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(1),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();


  char effstr[100];
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);

  resChi2Pass = make_plot(ptCutTag, lumi, nflpass, m, dataPass, modelPass,
      sigPass, bkgPass, fitResult, NsigPass, NbkgPass, effstr, passHist->GetEntries(), iBin, effType.Data(), etaRegion, kTRUE);

  resChi2Fail = make_plot(ptCutTag, lumi, nflfail, m, dataFail, modelFail,
      sigFail, bkgFail, fitResult, NsigFail, NbkgFail, effstr, failHist->GetEntries(), iBin, effType.Data(), etaRegion, kFALSE);


  delete dataCombined;
  delete dataPass;
  delete dataFail;
  delete sigPass;
  delete bkgPass;
  delete sigFail;
  delete bkgFail;
  delete histfile;
  delete histbkgQCDfile;
  delete histbkgTTfile;
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

//--------------------------------------------------------------------------------------------------
std::vector<float> calculateZYield(
        TH1D           *hYieldOSLowPU,        // histogram with reconstructed z candidates
        TH1D           *hYieldSSLowPU,        // histogram with reconstructed z candidates in same charge control region
        TH1D           *hYieldOSHighPU,       // histogram with reconstructed z candidates in high PU
        TH1D           *hYieldSSHighPU,       // histogram with reconstructed z candidates in high PU in same charge control region
        Int_t          sigOS=1,
        Int_t          bkgSS=2,
        Int_t          fitStrategy=1,
        Double_t       ptCutTag=27.,
        Double_t       ptCutProbe=27.,
        TH1D           *hPV=0,
        const TString  mcfilename="",
        const TString  outputDir="./",
        Int_t          iMeasurement = 0
){
    CPlot::sOutDir = outputDir;

    Int_t ptCut = 30;

    // fitStrategy:
    //  1: fit OS low PU region only
    //  2: fit bkg on SS low PU -> freeze bkg model -> fit bkg(scale) + signal on OS Low PU
    //  3: fit bkg + signal on OS and bkg on SS simultanious
    //  4: git bkg model in SS high PU -> freeze bkg model -> fit bkg(scale) + signal on OS Low PU


    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("OSLowPU",1);
    sample.defineType("SSLowPU",2);
    sample.defineType("OSHighPU",3);
    sample.defineType("SSHighPU",4);

    RooAbsData *dataOSLowPU  = new RooDataHist("dataOSLowPU","dataOSLowPU",RooArgSet(m),hYieldOSLowPU);
    RooAbsData *dataSSLowPU  = new RooDataHist("dataSSLowPU","dataSSLowPU",RooArgSet(m),hYieldSSLowPU);
    RooAbsData *dataOSHighPU = new RooDataHist("dataOSHighPU","dataOSHighPU",RooArgSet(m),hYieldOSHighPU);
    RooAbsData *dataSSHighPU = new RooDataHist("dataSSHighPU","dataSSHighPU",RooArgSet(m),hYieldSSHighPU);
    RooAbsData *dataSS = new RooDataHist("dataSS","dataSS",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("SSLowPU",*((RooDataHist*)dataSSLowPU)),
        RooFit::Import("SSHighPU",*((RooDataHist*)dataSSHighPU)));
    RooAbsData *dataOS = new RooDataHist("dataOS","dataOS",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("OSLowPU",*((RooDataHist*)dataOSLowPU)),
        RooFit::Import("OSHighPU",*((RooDataHist*)dataOSHighPU)));
    RooAbsData *dataLowPU = new RooDataHist("dataLowPU","dataLowPU",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("OSLowPU",*((RooDataHist*)dataOSLowPU)),
        RooFit::Import("SSLowPU",*((RooDataHist*)dataSSLowPU)));
    RooAbsData *dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("OSLowPU",*((RooDataHist*)dataOSLowPU)),
        RooFit::Import("SSLowPU",*((RooDataHist*)dataSSLowPU)),
        RooFit::Import("OSHighPU",*((RooDataHist*)dataOSHighPU)),
        RooFit::Import("SSHighPU",*((RooDataHist*)dataSSHighPU)));

    TFile *histfile = 0;
    if(sigOS==2) {
        generateTemplate_ZYield(mcfilename, ptCutTag, ptCutProbe, hPV, outputDir);
        histfile = new TFile(outputDir+"/histTemplates.root");
        assert(histfile);
    }

    CSignalModel     *modelSigOS = 0;
    CBackgroundModel *modelBkgSS = 0;

    Int_t nflOS=0, nflSS=0;

    switch(sigOS) {
      case 1:
        modelSigOS = new CBreitWignerConvCrystalBall(m, kFALSE, 0);
        nflOS += 4; break;
      case 2:
        TH1D *h = (TH1D*)histfile->Get("h_mass_zyield");
        assert(h);
        modelSigOS = new CMCTemplateConvGaussian(m,h,kFALSE,0);
        nflOS += 2; break;
    }
    switch(bkgSS) {
      case 1:
        modelBkgSS = new CExponential(m, kFALSE, 0);
        nflSS += 1; break;
      case 2:
        modelBkgSS = new CQuadratic(m, kFALSE, 0);
        nflSS += 3; break;
      case 3:
        modelBkgSS = new CQuadPlusExp(m, kFALSE, 0);
        nflSS += 4; break;
      case 4:
        modelBkgSS = new CDas(m, kFALSE, 0);
        nflSS += 4; break;
      case 5:
        modelBkgSS = new CDasPlusExp(m, kFALSE, 0);
        nflSS += 6; break;
    }

    Double_t NsigLowPUMax = hYieldOSLowPU->Integral();
    Double_t NbkgLowPUMax = hYieldSSLowPU->Integral();
    Double_t NsigHighPUMax = hYieldOSHighPU->Integral();
    Double_t NbkgHighPUMax = hYieldSSHighPU->Integral();

    RooRealVar nsigOSLowPU("nsigOS","signalfraction in OS", 0.99*NsigLowPUMax, 0., 1.5*NsigLowPUMax);
    RooRealVar nbkgOSLowPU("nbkgOS","backgroundfraction in OS", NbkgLowPUMax, 0., NsigLowPUMax);
    RooRealVar nbkgSSLowPU("nbkgSS","backgroundfraction in SS", NbkgLowPUMax, 0., 1.5*NbkgLowPUMax);

    RooRealVar nsigOSHighPU("nsigOSHighPU","signalfraction in OS high PU", 0.99*NsigHighPUMax, 0., 1.5*NsigHighPUMax);
    RooRealVar nbkgOSHighPU("nbkgOSHighPU","backgroundfraction in OS high PU", NbkgHighPUMax, 0., NsigHighPUMax);
    RooRealVar nbkgSSHighPU("nbkgSSHighPU","backgroundfraction in SS high PU", NbkgHighPUMax, 0., 1.5*NbkgHighPUMax);

    RooFormulaVar fr_lowPU("fr_lowPU","@0/(@0 + @1)",RooArgList(nbkgOSLowPU,nsigOSLowPU));
    RooFormulaVar fr_highPU("fr_highPU","@0/(@0 + @1)",RooArgList(nbkgOSHighPU,nsigOSHighPU));

    RooFormulaVar tf_lowPU("tf_lowPU","@0/@1",RooArgList(nbkgOSLowPU,nbkgSSLowPU));
    RooFormulaVar tf_highPU("tf_highPU","@0/@1",RooArgList(nbkgOSHighPU,nbkgSSHighPU));

    RooAddPdf *modelOSLowPU = new RooAddPdf("modelOSLowPU","Model for OS sample Low PU", RooArgList(*(modelSigOS->model),*(modelBkgSS->model)), RooArgList(nsigOSHighPU, nbkgOSLowPU));
    RooAddPdf *modelSSLowPU = new RooAddPdf("modelSSLowPU","Model for SS sample Low PU", RooArgList(*(modelBkgSS->model)), RooArgList(nbkgSSLowPU));
    RooAddPdf *modelOSHighPU = new RooAddPdf("modelOSHighPU","Model for OS high PU sample", RooArgList(*(modelSigOS->model),*(modelBkgSS->model)), RooArgList(nsigOSHighPU, nbkgOSHighPU));
    RooAddPdf *modelSSHighPU = new RooAddPdf("modelSSHighPU","Model for SS high PU sample", RooArgList(*(modelBkgSS->model)), RooArgList(nbkgSSHighPU));

    RooSimultaneous modelOS("modelOS","modelOS",sample);
    modelOS.addPdf(*modelOSLowPU,"OSLowPU");
    modelOS.addPdf(*modelOSHighPU,"OSHighPU");

    RooSimultaneous modelSS("modelSS","modelSS",sample);
    modelSS.addPdf(*modelSSLowPU,"SSLowPU");
    modelSS.addPdf(*modelSSHighPU,"SSHighPU");

    RooSimultaneous modelLowPU("modelLowPU","modelLowPU",sample);
    modelLowPU.addPdf(*modelSSLowPU,"SSLowPU");
    modelLowPU.addPdf(*modelOSLowPU,"OSLowPU");

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(*modelOSLowPU,"OSLowPU");
    totalPdf.addPdf(*modelOSHighPU,"OSHighPU");
    totalPdf.addPdf(*modelSSLowPU,"SSLowPU");
    totalPdf.addPdf(*modelSSHighPU,"SSHighPU");



    Int_t strategy = 1;
    RooFitResult *fitResult=0;

    if(fitStrategy == 5){
        // --- fit bkg on SS region low PU and high PU simultaniously
        fitResult = modelSS.fitTo(*dataSS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(1),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = modelSS.fitTo(*dataSS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(2),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());
    }


    if(fitStrategy==2 || fitStrategy==3 || fitStrategy==4 || fitStrategy==5){
        if(fitStrategy==2 || fitStrategy==4 || fitStrategy==5){
            modelBkgSS->freeze_all_parameters();
        }

        // --- fit on SS + OS together
        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            //RooFit::Extended(),
            RooFit::Strategy(strategy), // MINOS STRATEGY
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            //RooFit::Extended(),
            RooFit::Strategy(strategy), // MINOS STRATEGY
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());
    }


    //params = modelBkgSS->model->getVariables();
    //params->Print("v");

    Double_t resfr_lowPU = fr_lowPU.getVal();
    Double_t errfr_lowPU = fr_lowPU.getPropagatedError(*fitResult);
    Double_t restf_lowPU = tf_lowPU.getVal();
    Double_t errtf_lowPU = tf_lowPU.getPropagatedError(*fitResult);

    Double_t resfr_highPU = fr_highPU.getVal();
    Double_t errfr_highPU = fr_highPU.getPropagatedError(*fitResult);
    Double_t restf_highPU = tf_highPU.getVal();
    Double_t errtf_highPU = tf_highPU.getPropagatedError(*fitResult);


    // Plotting

    char ylabel[50];
    char yield[50];
    char nsigstr[50];
    char nbkgstr[50];
    char chi2str[50];
    char frstr[50];
    char tfstr[50];
    char sigstr[50];


    // --- plot opposite sign low PU region
    TCanvas *cOSL = MakeCanvas("cOSL","cOSL",720,540);
    cOSL->SetWindowPosition(cOSL->GetWindowTopX()+cOSL->GetBorderSize()+800,0);


    RooPlot *mframeOSL = m.frame(Bins(massBin));
    dataOSLowPU->plotOn(mframeOSL,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelOSLowPU->plotOn(mframeOSL, Components(*(modelSigOS->model)), LineColor(9));
    modelOSLowPU->plotOn(mframeOSL, Components(*(modelBkgSS->model)), LineColor(8));

    modelOSLowPU->plotOn(mframeOSL, LineColor(kRed));

    sprintf(yield,"%u Events",(Int_t)hYieldOSLowPU->GetEntries());
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",nsigOSLowPU.getVal(),nsigOSLowPU.getPropagatedError(*fitResult));
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgOSLowPU.getVal(),nbkgOSLowPU.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeOSL->chiSquare(nflOS+1));
    sprintf(frstr,"fr = %.3f #pm %.3f",resfr_lowPU, errfr_lowPU);
    sprintf(tfstr,"tf = %.3f #pm %.3f",restf_lowPU, errtf_lowPU);


    CPlot plotOSL("plot_OSLowPU_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement),
        mframeOSL, "Opposite Sign - Low PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotOSL.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotOSL.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotOSL.AddTextBox(std::to_string(ptCut)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotOSL.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotOSL.AddTextBox(modelSigOS->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotOSL.AddTextBox(modelBkgSS->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotOSL.AddTextBox(0.70,0.62,0.94,0.90,0,kBlack,-1,5,chi2str, nbkgstr, nsigstr, frstr, tfstr);

    plotOSL.Draw(cOSL,kTRUE,"png");

    plotOSL.SetName("plot_OSLowPU_logscale_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS));

    double ymin = hYieldOSLowPU->GetMinimum();
    if(ymin <= 0.) ymin = 1.;

    plotOSL.SetYRange(0.5 * ymin, hYieldOSLowPU->GetMaximum()*2);
    plotOSL.SetLogy();
    plotOSL.Draw(cOSL,kTRUE,"png");


    // --- plot opposite sign region high PU
    TCanvas *cOSH = MakeCanvas("cOSH","cOSH",720,540);
    cOSH->SetWindowPosition(cOSH->GetWindowTopX()+cOSH->GetBorderSize()+800,0);


    RooPlot *mframeOSH = m.frame(Bins(massBin));
    dataOSHighPU->plotOn(mframeOSH,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelOSHighPU->plotOn(mframeOSH, Components(*(modelSigOS->model)), LineColor(9));
    modelOSHighPU->plotOn(mframeOSH, Components(*(modelBkgSS->model)), LineColor(8));

    modelOSHighPU->plotOn(mframeOSH, LineColor(kRed));

    sprintf(yield,"%u Events",(Int_t)hYieldOSHighPU->GetEntries());
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",nsigOSHighPU.getVal(),nsigOSHighPU.getPropagatedError(*fitResult));
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgOSHighPU.getVal(),nbkgOSHighPU.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeOSH->chiSquare(nflOS+1));
    sprintf(frstr,"fr = %.3f #pm %.3f",resfr_highPU, errfr_highPU);
    sprintf(tfstr,"tf = %.3f #pm %.3f",restf_highPU, errtf_highPU);


    CPlot plotOSH("plot_OSHighPU_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement),
        mframeOSH,"Opposite Sign - High PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotOSH.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotOSH.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotOSH.AddTextBox(std::to_string(ptCut)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotOSH.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotOSH.AddTextBox(modelSigOS->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotOSH.AddTextBox(modelBkgSS->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotOSH.AddTextBox(0.70,0.62,0.94,0.90,0,kBlack,-1,5,chi2str, nbkgstr, nsigstr, frstr, tfstr);

    plotOSH.Draw(cOSH,kTRUE,"png");

    plotOSH.SetName("plot_OSHighPU_logscale_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement));

    ymin = hYieldOSHighPU->GetMinimum();
    if(ymin <= 0.) ymin = 1.;

    plotOSH.SetYRange(0.5 * ymin, hYieldOSHighPU->GetMaximum()*2);
    plotOSH.SetLogy();
    plotOSH.Draw(cOSH,kTRUE,"png");


    // --- plot pulls in OS low PU region
    TCanvas *cPullsOS = MakeCanvas("cResiOS","cResiOS",720,540);
    cPullsOS->SetWindowPosition(cPullsOS->GetWindowTopX()+cPullsOS->GetBorderSize()+800,0);

    RooHist* pullsHist = mframeOSL->pullHist();
    pullsHist->SetPointError(0.,0.,0.,0.);

    RooPlot *mframePullsOS = m.frame(Bins(massBin));
    mframePullsOS->addPlotable(pullsHist, "P");
    CPlot plotPullsOS("plot_ResidOS_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement),
        mframePullsOS,"Opposite Sign - Low PU","tag-probe mass [GeV/c^{2}]","pulls");

    plotPullsOS.Draw(cOSL,kTRUE,"png");


    // --- plot same sign region
    TCanvas *cSS = MakeCanvas("cSS","cSS",720,540);
    cSS->SetWindowPosition(cSS->GetWindowTopX()+cSS->GetBorderSize()+800,0);

    RooPlot *mframeSS = m.frame(Bins(massBin));
    dataSSLowPU->plotOn(mframeSS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelSSLowPU->plotOn(mframeSS, Components(*(modelBkgSS->model)), LineColor(8));

    sprintf(yield,"%u Events",(Int_t)hYieldSSLowPU->GetEntries());
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgSSLowPU.getVal(),nbkgSSLowPU.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeSS->chiSquare(nflSS));

    CPlot plotSS("plot_SSLowPU_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement),
        mframeSS,"Same Sign - Low PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotSS.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotSS.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotSS.AddTextBox(std::to_string(ptCut)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotSS.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotSS.AddTextBox(modelBkgSS->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotSS.AddTextBox(0.70,0.79,0.94,0.90,0,kBlack,-1,2,chi2str, nbkgstr);

    plotSS.Draw(cSS,kTRUE,"png");

    // --- plot same sign high pileup region
    TCanvas *cSSHighPU = MakeCanvas("cSSHighPU","cSSHighPU",720,540);
    cSSHighPU->SetWindowPosition(cSSHighPU->GetWindowTopX()+cSSHighPU->GetBorderSize()+800,0);

    RooPlot *mframeSSHighPU = m.frame(Bins(massBin));
    dataSSHighPU->plotOn(mframeSSHighPU,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelSSHighPU->plotOn(mframeSSHighPU, Components(*(modelBkgSS->model)), LineColor(8));

    sprintf(yield,"%u Events",(Int_t)hYieldSSHighPU->GetEntries());
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgSSHighPU.getVal(),nbkgSSHighPU.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeSSHighPU->chiSquare(nflSS));

    CPlot plotSSHighPU("plot_SSHighPU_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS)+"_"+std::to_string(iMeasurement),
        mframeSSHighPU,"Same Sign - High PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotSSHighPU.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotSSHighPU.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotSSHighPU.AddTextBox(std::to_string(ptCut)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotSSHighPU.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotSSHighPU.AddTextBox(modelBkgSS->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotSSHighPU.AddTextBox(0.70,0.79,0.94,0.90,0,kBlack,-1,2,chi2str, nbkgstr);

    const Double_t yMax = hYieldSSHighPU->GetMaximum() + hYieldSSHighPU->GetBinError(hYieldSSHighPU->GetMaximumBin());
    const Double_t yMin = std::max(0., hYieldSSHighPU->GetMinimum() - hYieldSSHighPU->GetBinError(hYieldSSHighPU->GetMinimumBin()));
    plotSSHighPU.SetYRange(yMin, yMax + 0.4*(yMax-yMin));
    plotSSHighPU.Draw(cSSHighPU,kTRUE,"png");

    // prepare return stuff
    std::vector<float> result = {};

    result.push_back(restf_lowPU);
    result.push_back(errtf_lowPU);
    result.push_back(errtf_lowPU);
    result.push_back(mframeOSL->chiSquare(4));
    result.push_back(mframeSS->chiSquare(3));

    return result;

}
