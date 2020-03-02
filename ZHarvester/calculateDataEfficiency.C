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


using namespace std;

// Set up fiducial region
const Float_t massLo  = 66.;
const Float_t massHi  = 116.;
const UInt_t  massBin = 50;

const Float_t etaCutTag   = 2.4;
const Float_t etaCutProbe = 2.4;
const Float_t etaBound    = 0.9;

// Set up pile-up bounds available get MC template for fitting
const Int_t minPU = 1;
const Int_t maxPU = 60;

void generateTemplate(
        const TString mcfilename,
        const Float_t ptCutTag, const Float_t ptCutProbe, TH1D *hPV=0, const TString outputDir="");

void performCount(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
        const Float_t ptCutTag, const Float_t ptCutProbe,
        TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

void performFit(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2Pass, Double_t &resChi2Fail, TH1D *passHist, TH1D *failHist,
        const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
        const Float_t ptCutTag, const Float_t ptCutProbe, const TString outputDir, const TString bkgQCDTemplate, const TString bkgTTTemplate,
        TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

std::vector<double> preFit(TH1D *failHist);

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
                const TString format="png"  // plot format
){

  CPlot::sOutDir = outputDir;

  TCanvas *cyield = MakeCanvas("cyield","cyield",720,540);
  cyield->SetWindowPosition(cyield->GetWindowTopX()+cyield->GetBorderSize()+800,0);


  Double_t NsigMax = h_yield->GetEntries();
  if(sigMod == 0){      // perform count - expect 1% faces
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

  if(sigMod == 1){     // perform fit
    sigModel = new CBreitWignerConvCrystalBall(m, kTRUE, 0);
    nfl += 4;
    if(bkgMod == 1){
      bkgModel = new CExponential(m, kTRUE, 0);
      nfl += 6;
    }
    if(bkgMod == 5){
      bkgModel = new CDasPlusExp(m, kTRUE, 0);
      nfl += 6;
    }
    else{
      std::cout <<"ERROR: unvalid background mode for computing Z yield"<< std::endl;
    }
  }
  else{
    std::cout <<"ERROR: unvalid signal mode for computing Z yield"<< std::endl;
  }

  RooAbsData *data = 0;
  data = new RooDataHist("ZReco","ZReco",RooArgList(m),h_yield);

  RooRealVar Nsig("Nsig","sigYield",NsigMax,0.,1.5*NsigMax);
  RooRealVar Nbkg("Nbkg","bkgYield",0.01*NsigMax,0.,NsigMax);
  RooAddPdf modelPdf("model","Z sig+bkg",RooArgList(*(sigModel->model),*(bkgModel->model)),RooArgList(Nsig,Nbkg));

  RooFitResult *fitResult=0;
  fitResult = modelPdf.fitTo(*data,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(1), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet()),
                             RooFit::Save());

  double resNsig  = Nsig.getVal();
  double resErrl = fabs(Nsig.getErrorLo());
  double resErrh = Nsig.getErrorHi();
  double resChi2 = 0.;

  char binlabelx[100];
  char binlabely[100];
  char effstr[100];
  char lumitext[100];
  char ylabel[50];
  char pname[50];
  char yield[50];
  char nsigstr[100];
  char nbkgstr[100];
  char chi2str[100];

  sprintf(binlabelx, "0.0 < |#eta| < 2.4");
  sprintf(binlabely, "%i GeV/c < p_{T} < 13000 GeV/c",(Int_t)ptCutTag);
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
  sprintf(ylabel,"Events / 1 GeV/c^{2}");

  RooPlot *mframe = m.frame(Bins(massBin));
  data->plotOn(mframe,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  if(bkgMod>0)
    modelPdf.plotOn(mframe,Components("backgroundPass_0"),LineStyle(kDashed),LineColor(kRed));
  modelPdf.plotOn(mframe);

  double a = Nsig.getVal(), aErr = Nsig.getPropagatedError(*fitResult);
  double b = Nbkg.getVal(), bErr = Nbkg.getPropagatedError(*fitResult);
  double fpr = b/(a+b);
  sprintf(effstr,"#frac{bkg}{sig+bkg} = %.4f #pm %.4f",fpr,1./((a+b)*(a+b))*(std::abs(b*aErr) + std::abs(a*bErr)));
  sprintf(pname,"ZYield_inclusive_%d", iBin);
  sprintf(yield,"%u Events",(Int_t)NsigMax);
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",a,aErr);
  resChi2 = mframe->chiSquare(nfl);
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",resChi2);

  if(bkgMod>0)
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",b,bErr);
  CPlot plotPass(pname,mframe,"Z Yield","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
  plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
  plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
  if(bkgMod>0) {
    plotPass.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  } else {
    plotPass.AddTextBox(0.70,0.73,0.94,0.83,0,kBlack,-1,1,nsigstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  }
  plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);

  plotPass.Draw(cyield,kTRUE,format);

  delete sigModel;
  delete bkgModel;
  delete data;
  delete cyield;

  std::vector<float> resultEff = {};

  resultEff.push_back(resNsig);
  resultEff.push_back(resErrl);
  resultEff.push_back(resErrh);
  resultEff.push_back(resChi2);
  resultEff.push_back(fpr);

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
  // Calculate efficiency and save results to plots and histograms
  TCanvas *cpass = MakeCanvas("cpass","cpass",720,540);
  cpass->SetWindowPosition(cpass->GetWindowTopX()+cpass->GetBorderSize()+800,0);
  TCanvas *cfail = MakeCanvas("cfail","cfail",720,540);
  cfail->SetWindowPosition(cfail->GetWindowTopX()+cfail->GetBorderSize()+800,cpass->GetWindowTopX()+cfail->GetBorderSize()+540);

  Double_t eff  = 0.;
  Double_t errl = 0.;
  Double_t errh = 0.;
  Double_t chi2pass = 999.;
  Double_t chi2fail = 999.;

  if(sigModPass == 0){
    performCount(eff, errl, errh, h_mass_pass, h_mass_fail, ptCutTag, ptCutProbe,
      cpass, cfail, effType, etaRegion, iBin, lumi, format);
    //cout << effType.Data() << ": " << eff << " + " << errh << " - " << errl << endl;
  }else{
    performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail,
      sigModPass, bkgModPass, sigModFail, bkgModFail,
      ptCutTag, ptCutProbe, outputDir, bkgQCDFilename, bkgTTFilename,
      cpass, cfail, effType, etaRegion, iBin, lumi, format);
    //cout << effType.Data() << ": " << eff << " + " << errh << " - " << errl << endl;
  }
  delete cpass;
  delete cfail;

  std::vector<float> resultEff = {};

  resultEff.push_back(eff);
  resultEff.push_back(errl);
  resultEff.push_back(errh);
  resultEff.push_back(chi2pass);
  resultEff.push_back(chi2fail);

  return resultEff;
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
void performCount(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
	const Float_t ptCutTag, const Float_t ptCutProbe,
	TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format
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
	TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format
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

// Inital Pdfs for models
// Signal Model:
// 	    1: Breit-Wigner convolved with Crystal Ball function
// 	    2: MC template convolved with Gaussian
// Background Model:
//      1: exponential model
//      2: quadratic model
//      3: exponential + quadratic model
//      4: Das (wide gaussian with exponential tails) model
//      5: Das + decaying exponential model
//      6: bkg QCD template
//      7: bkg QCD + ttbar template


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

  // for plotting
  char fsigpass[50];
  char fbkgpass[50];
  char fsigfail[50];
  char fbkgfail[50];

  switch(sigpass) {
    case 1:
      sigPass = new CBreitWignerConvCrystalBall(m, kTRUE, etaRegion);
      nflpass += 4; break;
    case 2:
      TH1D *h = (TH1D*)histfile->Get(Form("h_mass_pass_%s", etaRegion ? "forward" : "central"));
      assert(h);
      sigPass = new CMCTemplateConvGaussian(m,h,kTRUE,etaRegion);
      nflpass += 2; break;
  }

  switch(bkgpass) {
    case 1:
      bkgPass = new CExponential(m, kTRUE, etaRegion);
      nflpass += 1; break;
    case 2:
      bkgPass = new CQuadratic(m, kTRUE, etaRegion, 0.,0.,0.,0.,0.,0.);
      nflpass += 3; break;
    case 3:
      bkgPass = new CQuadPlusExp(m, kTRUE, etaRegion, 0.,0.,0.,0.,0.,0.);
      nflpass += 4; break;
    case 4:
      bkgPass = new CDas(m, kTRUE, etaRegion);
      nflpass += 4; break;
    case 5:
      bkgPass = new CDasPlusExp(m, kTRUE, etaRegion);
      nflpass += 6; break;
    case 6:
      bkgPass = new CQCD(m, hbkgQCDPass, kTRUE, etaRegion);
      nflpass += 1; break;
    case 7:
      bkgPass = new CQCDPlusTT(m, hbkgQCDPass, hbkgTTPass, kTRUE, etaRegion);
      nflpass += 2; break;
  }

  switch(sigfail) {
    case 1:
      sigFail = new CBreitWignerConvCrystalBall(m, kFALSE, etaRegion);
      nflfail += 4; break;
    case 2:
      TH1D *h = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
      assert(h);
      sigFail = new CMCTemplateConvGaussian(m,h,kFALSE,etaRegion);
      nflfail += 2; break;
  }

  switch(bkgfail) {
    case 1:
      bkgFail = new CExponential(m, kFALSE, etaRegion);
      nflfail += 1; break;
    case 2:
      bkgFail = new CQuadratic(m, kFALSE, etaRegion, vBkgPars[0], vBkgPars[1], vBkgPars[2], vBkgPars[3], vBkgPars[4], vBkgPars[5]);
      nflfail += 3; break;
    case 3:
      bkgFail = new CQuadPlusExp(m, kFALSE, etaRegion, vBkgPars[0], vBkgPars[1], vBkgPars[2], vBkgPars[3], vBkgPars[4], vBkgPars[5]);
      nflfail += 4; break;
    case 4:
      bkgFail = new CDas(m, kFALSE, etaRegion);
      nflfail += 4; break;
    case 5:
      bkgFail = new CDasPlusExp(m, kFALSE, etaRegion);
      nflfail += 6; break;
    case 6:
      bkgFail = new CQCD(m, hbkgQCDFail, kFALSE, etaRegion);
      nflfail += 1; break;
    case 7:
      bkgFail = new CQCDPlusTT(m, hbkgQCDFail, hbkgTTFail, kFALSE, etaRegion);
      nflfail += 2; break;
  }


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
  RooAddPdf *modelPass=0, *modelFail=0;

  modelPass = new RooAddPdf("modelPass","Model for PASS sample",
                              (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
                              (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
  modelFail = new RooAddPdf("modelFail","Model for FAIL sample",RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));

  RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
  totalPdf.addPdf(*modelPass,"Pass");
  totalPdf.addPdf(*modelFail,"Fail");

  RooFitResult *fitResult=0;
  Int_t strategy = 2;
  if(effType == "Sel" or effType == "HLT") {
    Nsig.setRange(0,2.0*NsigMax);
    strategy=1;
  }
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             //RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();

  // Plot tag and passing- and failing- probe mass distribution
  char bkgpassname[50];
  char bkgfailname[50];
  char sigpassname[50];
  char sigfailname[50];
  char binlabelx[100];
  char binlabely[100];
  char effstr[100];
  char lumitext[100];
  char ylabel[50];
  char pname[50];
  char yield[50];
  char nsigstr[100];
  char nbkgstr[100];
  char chi2str[100];
  char expstr[50];

  if(!etaRegion) sprintf(binlabelx, "0.0 < |#eta| < 0.9");
  else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");

  sprintf(binlabely, "%i GeV/c < p_{T} < 13000 GeV/c",(Int_t)ptCutProbe);
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
  sprintf(ylabel,"Events / 1 GeV/c^{2}");

  RooPlot *mframePass = m.frame(Bins(massBin));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  if(bkgpass>0){
    sprintf(bkgpassname, "backgroundPass_%d",etaRegion);
    modelPass->plotOn(mframePass,Components(bkgpassname),LineColor(8));
    if(bkgpass == 3 or bkgpass == 5 or bkgpass == 7){
      sprintf(bkgpassname, "bkg_pdf1Pass_%d",etaRegion);
      modelPass->plotOn(mframePass,Components(bkgpassname),LineStyle(kDashed),LineColor(8));
      sprintf(bkgpassname, "bkg_pdf2Pass_%d",etaRegion);
      modelPass->plotOn(mframePass,Components(bkgpassname),LineStyle(kDashed),LineColor(8));
    }
  }
  if(sigpass>0){
    sprintf(sigpassname, "signalPass_%d",etaRegion);
    modelPass->plotOn(mframePass,Components(sigpassname),LineColor(9));
  }
  modelPass->plotOn(mframePass, LineColor(2));

  double a = NsigPass.getVal(), b = NbkgPass.getVal();
  sprintf(yield,"%u Events",(Int_t)passHist->GetEntries());
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflpass));resChi2Pass = mframePass->chiSquare(nflpass);
  if(bkgpass>0)
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgPass.getVal(),NbkgPass.getPropagatedError(*fitResult));

  for(int i=0; i <=1; i++){
    // do two plots, one with linear scale and one with log scale
    if(i==1) sprintf(expstr, "_exp");
    else     sprintf(expstr, "");

    sprintf(pname,"%s_%s_pass_%d%s", effType.Data(), etaRegion ? "forward" : "central", iBin, expstr);
    CPlot plotPass(pname,mframePass,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);

    plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotPass.AddTextBox(sigPass->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotPass.AddTextBox(bkgPass->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90, 0,kBlack,-1);
    if(bkgpass>0) {
      plotPass.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
      plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
    } else {
      plotPass.AddTextBox(0.70,0.73,0.94,0.83,0,kBlack,-1,1,nsigstr);
      plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
    }
    plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);

    if(i==1){
      double ymin = passHist->GetMinimum();
      if(ymin <= 0.) ymin = 1.;

      plotPass.SetYRange(0.1 * ymin, passHist->GetMaximum()*10);
      plotPass.SetLogy();
    }

    plotPass.Draw(cpass,kTRUE,format);
  }

  RooPlot *mframeFail = m.frame(Bins(massBin));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  if(bkgfail>0){
    sprintf(bkgfailname, "backgroundFail_%d",etaRegion);
    modelFail->plotOn(mframeFail,Components(bkgfailname),LineColor(8));
    if(bkgfail == 3 or bkgfail == 5 or bkgfail == 7){
      sprintf(bkgfailname, "bkg_pdf1Fail_%d",etaRegion);
      modelFail->plotOn(mframeFail,Components(bkgfailname),LineStyle(kDashed),LineColor(8));
      sprintf(bkgfailname, "bkg_pdf2Fail_%d",etaRegion);
      modelFail->plotOn(mframeFail,Components(bkgfailname),LineStyle(kDashed),LineColor(8));
    }
  }
  if(sigfail>0){
    sprintf(sigfailname, "signalFail_%d",etaRegion);
    modelFail->plotOn(mframeFail,Components(sigfailname),LineColor(9));
  }
  modelFail->plotOn(mframeFail,LineColor(2));

  double f = NsigFail.getVal(), d = NbkgFail.getVal();
  sprintf(yield,"%u Events",(Int_t)failHist->GetEntries());
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
  sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeFail->chiSquare(nflfail));resChi2Fail = mframeFail->chiSquare(nflfail);

  for(int i=0; i <=1; i++){
    // do two plots, one with linear scale and one with log scale
    if(i==1) sprintf(expstr, "_exp");
    else     sprintf(expstr, "");

    sprintf(pname,"%s_%s_fail_%d%s", effType.Data(), etaRegion ? "forward" : "central", iBin, expstr);
    CPlot plotFail(pname,mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);

    plotFail.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotFail.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotFail.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotFail.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
    plotFail.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
    plotFail.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,-1);
    plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
    plotFail.AddTextBox(sigFail->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotFail.AddTextBox(bkgFail->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);

    if(i==1){
      double ymin = failHist->GetMinimum();
      if(ymin <= 0.) ymin = 1.;

      plotFail.SetYRange(0.5 * ymin,failHist->GetMaximum()*2);
      plotFail.SetLogy();
    }

    plotFail.Draw(cfail,kTRUE,format);

  }

  delete modelPass;
  delete modelFail;
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

//  std::vector<float> v = {1.,0.,1.,0.,1.,0.};return v;
  TH1D *h = new TH1D("h", "", massBin, massLo, massHi);
  TF1 *fq = new TF1("fq", "[0]+[1]*x+[2]*x*x", massLo, massHi);
  TF1 *fl = new TF1("fl", "[0]+[1]*x", massLo, massHi);

  for(int i = 0; i < 15; i++){
    h->SetBinContent(i+1, failHist->GetBinContent(i+1));
    h->SetBinError(i+1, failHist->GetBinError(i+1));
  }
  for(int i = 35; i < 50; i++){
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
        TH1D           *hYieldOS,             // histogram with reconstructed z candidates
        TH1D           *hYieldSS,             // histogram with reconstructed z candidates in same charge control region
        TH1D           *hYieldSSHighPU,       // histogram with reconstructed z candidates in same charge control region
        Int_t          sigOS=1,
        Int_t          bkgSS=2,
        Int_t          fitStrategy=1,
        Double_t       ptCutTag=27.,
        Double_t       ptCutProbe=27.,
        TH1D           *hPV=0,
        const TString  mcfilename="",
        const TString  outputDir="./"
){
    CPlot::sOutDir = outputDir;

    Int_t ptCut = 30;

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("OS",1);
    sample.defineType("SS",2);
    sample.defineType("SSHighPU",3);

    RooAbsData *dataOS       = new RooDataHist("dataOS","dataOS",RooArgSet(m),hYieldOS);
    RooAbsData *dataSS       = new RooDataHist("dataSS","dataSS",RooArgSet(m),hYieldSS);
    RooAbsData *dataSSHighPU = new RooDataHist("dataSSHighPU","dataSSHighPU",RooArgSet(m),hYieldSSHighPU);
    RooAbsData *dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("OS",*((RooDataHist*)dataOS)),
        RooFit::Import("SS",*((RooDataHist*)dataSS)),
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

    Double_t NsigMax = hYieldOS->Integral();
    Double_t NbkgMax = hYieldSS->Integral();
    Double_t NbkgHighPUMax = hYieldSSHighPU->Integral();

    RooRealVar nsigOS("nsigOS","signalfraction in OS", 0.99*NsigMax, 0., 1.5*NsigMax);
    RooRealVar nbkgOS("nbkgOS","backgroundfraction in OS", NbkgMax, 0., NsigMax);
    RooRealVar nbkgSS("nbkgSS","backgroundfraction in SS", NbkgMax, 0., 1.5*NbkgMax);
    RooRealVar nbkgSSHighPU("nbkgSSHighPU","backgroundfraction in SS high PU", NbkgHighPUMax, 0., 1.5*NbkgHighPUMax);

    RooFormulaVar fr("fr","@0/(@0 + @1)",RooArgList(nbkgOS,nsigOS));
    RooFormulaVar tf("tf","@0/@1",RooArgList(nbkgOS,nbkgSS));

    RooAddPdf *modelOS = new RooAddPdf("modelOS","Model for OS sample", RooArgList(*(modelSigOS->model),*(modelBkgSS->model)), RooArgList(nsigOS, nbkgOS));
    RooAddPdf *modelSS = new RooAddPdf("modelSS","Model for SS sample", RooArgList(*(modelBkgSS->model)), RooArgList(nbkgSS));
    RooAddPdf *modelSSHighPU = new RooAddPdf("modelSSHighPU","Model for SS high PU sample", RooArgList(*(modelBkgSS->model)), RooArgList(nbkgSSHighPU));

    Int_t strategy = 1;
    RooFitResult *fitResult=0;

    if(fitStrategy == 2){
        // --- fit on SS region
        fitResult = modelSS->fitTo(*dataSS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(1),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = modelSS->fitTo(*dataSS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(2),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());
    }
    else if(fitStrategy == 4){
        // --- fit on SS region
        fitResult = modelSSHighPU->fitTo(*dataSSHighPU,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(1),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = modelSSHighPU->fitTo(*dataSSHighPU,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(2),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());
    }


    RooArgSet* params = modelBkgSS->model->getVariables();
    params->Print("v");

    if(fitStrategy == 1){
        // --- fit on OS region
        fitResult = modelOS->fitTo(*dataOS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(1),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());

        fitResult = modelOS->fitTo(*dataOS,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(2),
            //RooFit::Minos(RooArgSet(eff)),
            RooFit::Save());
    }


    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(*modelOS,"OS");
    totalPdf.addPdf(*modelSS,"SS");

    if(fitStrategy==2 || fitStrategy==3 || fitStrategy==4){
        if(fitStrategy==2 || fitStrategy==4){
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


    params = modelBkgSS->model->getVariables();
    params->Print("v");

    Double_t resfr = fr.getVal();
    Double_t errfr = fr.getPropagatedError(*fitResult);
    Double_t restf = tf.getVal();
    Double_t errtf = tf.getPropagatedError(*fitResult);



    // Plotting

    char ylabel[50];
    char yield[50];
    char nsigstr[50];
    char nbkgstr[50];
    char chi2str[50];
    char frstr[50];
    char tfstr[50];
    char sigstr[50];


    // --- plot opposite sign region
    TCanvas *cOS = MakeCanvas("cOS","cOS",720,540);
    cOS->SetWindowPosition(cOS->GetWindowTopX()+cOS->GetBorderSize()+800,0);


    RooPlot *mframeOS = m.frame(Bins(massBin));
    dataOS->plotOn(mframeOS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelOS->plotOn(mframeOS, Components(*(modelSigOS->model)), LineColor(9));
    modelOS->plotOn(mframeOS, Components(*(modelBkgSS->model)), LineColor(8));

    modelOS->plotOn(mframeOS, LineColor(kRed));

    sprintf(yield,"%u Events",(Int_t)hYieldOS->GetEntries());
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",nsigOS.getVal(),nsigOS.getPropagatedError(*fitResult));
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgOS.getVal(),nbkgOS.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeOS->chiSquare(nflOS+1));
    sprintf(frstr,"fr = %.3f #pm %.3f",resfr, errfr);
    sprintf(tfstr,"tf = %.3f #pm %.3f",restf, errtf);


    CPlot plotOS("plot_OS_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS),mframeOS,"Opposite Sign - Low PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotOS.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotOS.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotOS.AddTextBox(std::to_string(ptCut)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotOS.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotOS.AddTextBox(modelSigOS->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotOS.AddTextBox(modelBkgSS->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotOS.AddTextBox(0.70,0.62,0.94,0.90,0,kBlack,-1,5,chi2str, nbkgstr, nsigstr, frstr, tfstr);

    plotOS.Draw(cOS,kTRUE,"png");

    plotOS.SetName("plot_OS_logscale_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS));

    double ymin = hYieldOS->GetMinimum();
    if(ymin <= 0.) ymin = 1.;

    plotOS.SetYRange(0.5 * ymin, hYieldOS->GetMaximum()*2);
    plotOS.SetLogy();
    plotOS.Draw(cOS,kTRUE,"png");


    // --- plot pulls in OS region
    TCanvas *cPullsOS = MakeCanvas("cResiOS","cResiOS",720,540);
    cPullsOS->SetWindowPosition(cPullsOS->GetWindowTopX()+cPullsOS->GetBorderSize()+800,0);

    RooHist* pullsHist = mframeOS->pullHist();
    pullsHist->SetPointError(0.,0.,0.,0.);

    RooPlot *mframePullsOS = m.frame(Bins(massBin));
    mframePullsOS->addPlotable(pullsHist, "P");
    CPlot plotPullsOS("plot_ResidOS_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS),mframePullsOS,"Opposite Sign - Low PU","tag-probe mass [GeV/c^{2}]","pulls");

    plotPullsOS.Draw(cOS,kTRUE,"png");


    // --- plot same sign region
    TCanvas *cSS = MakeCanvas("cSS","cSS",720,540);
    cSS->SetWindowPosition(cSS->GetWindowTopX()+cSS->GetBorderSize()+800,0);

    RooPlot *mframeSS = m.frame(Bins(massBin));
    dataSS->plotOn(mframeSS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelSS->plotOn(mframeSS, Components(*(modelBkgSS->model)), LineColor(8));

    sprintf(yield,"%u Events",(Int_t)hYieldSS->GetEntries());
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgSS.getVal(),nbkgSS.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeSS->chiSquare(nflSS));

    CPlot plotSS("plot_SS_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS),mframeSS,"Same Sign - Low PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

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

    CPlot plotSSHighPU("plot_SSHighPU_"+std::to_string(100*fitStrategy+10*bkgSS+sigOS),mframeSSHighPU,"Same Sign - High PU","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

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

    result.push_back(restf);
    result.push_back(errtf);
    result.push_back(errtf);
    result.push_back(mframeOS->chiSquare(4));
    result.push_back(mframeSS->chiSquare(3));

    return result;

}
