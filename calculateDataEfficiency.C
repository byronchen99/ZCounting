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

const Float_t ptCut    = 27.;
const Float_t etaCut   = 2.4;
const Float_t etaBound = 0.9;

// Set up pile-up bounds available get MC template for fitting
const Int_t minPU = 1;
const Int_t maxPU = 60;

void generateTemplate(
        const TString mcfilename, const TString purwDir, const Int_t meanPU);

void performCount(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
        TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

void performFit(
        Double_t &resEff, Double_t &resErrl, Double_t &resErrh, Double_t &resChi2Pass, Double_t &resChi2Fail, TH1D *passHist, TH1D *failHist,
        const Int_t sigpass, const Int_t bkgpass, const Int_t sigfail, const Int_t bkgfail,
        TCanvas *cpass, TCanvas *cfail, const TString effType, const Bool_t etaRegion, const Int_t iBin, const Float_t lumi, const TString format);

std::vector<double> preFit(TH1D *failHist);


float getZyield(
                const TString inputFile,
                const TString hName,        //Name of the histogram to take the Z yield from
                const TString runNum,
                const Int_t   startLS,
                const Int_t   endLS,
                const Int_t   sigMod=0,
                const Int_t   bkgMod=0
){
  TFile *infile = new TFile(inputFile);


  TH1F *h_yield = (TH1F*)infile->Get("DQMData/Run "+runNum+"/ZCounting/Run summary/Histograms/"+hName);

  //if(sigMod == 0){          // perform count
  return h_yield->Integral(startLS,endLS);
  /*}                         
  /else if(sigMod == 1){     // perform fit -- under construction
    CSignalModel     *sigModel = 0;
    CBackgroundModel *bkgModel = 0;

    sigModel = new CBreitWignerConvCrystalBall(m,kTRUE);
    if(bkgMod == 1){
      bkgModel = new CExponential(m, kTRUE, etaRegion);
    }
    else{
      std::cout <<"ERROR: unvalid background mode for computing Z yield"<< std::endl;
    }
  }
  else{
    std::cout <<"ERROR: unvalid signal mode for computing Z yield"<< std::endl;
    return 0;
  }*/
}


std::vector<float> calculateDataEfficiency(
		const TString inputFile,    // DQMIO file 
		const TString outputDir,    // output directory
		const TString runNum,       // Runnum
		const Int_t   iBin,         // Label of measurement in currect run
		const Int_t   startLS,      // starting LS for currect measurement
		const Int_t   endLS,        // ending LS for currect measurement
		const Float_t meanPU,       // mean #PU, only needed if use MC template as fitting model
		const TString effType,      // "HLT" or "SIT" or "Glo" or "Sta" or "Trk"
		const Bool_t  etaRegion,    // 0 for barrel, 1 for forward 
		const Int_t   sigModPass,   // signal extraction method for PASS sample
		const Int_t   bkgModPass,   // background model for PASS sample
		const Int_t   sigModFail,   // signal extraction method for FAIL sample      
		const Int_t   bkgModFail,   // background model for FAIL sample
		const Float_t lumi=10.,     // luminosity for plot label
                const TString mcfilename="",// ROOT file containing MC events to generate templates from
                const TString purwDir="",
		const TString format="png" // plot format
){

  gSystem->mkdir(outputDir,kTRUE);
  CPlot::sOutDir = outputDir + TString("/Run") + runNum + TString("/plots");

  // Generate histogram templates from MC if necessary
  if(sigModPass==2 || sigModFail==2) generateTemplate(mcfilename, purwDir, Int_t(meanPU));

  // Load the DQMIO file and corresponding histograms 
  TFile *infile = new TFile(inputFile);

  TString region = etaRegion ? "forward" : "central";
  TH2F *h_mass_lumi_pass = (TH2F*)infile->Get("DQMData/Run "+runNum+"/ZCounting/Run summary/Histograms/h_mass_"+effType+"_pass_"+region);
  TH2F *h_mass_lumi_fail = (TH2F*)infile->Get("DQMData/Run "+runNum+"/ZCounting/Run summary/Histograms/h_mass_"+effType+"_fail_"+region);

  TH1D *h_mass_pass = h_mass_lumi_pass->ProjectionY("h_mass_pass", startLS, endLS, "e"); 
  TH1D *h_mass_fail = h_mass_lumi_fail->ProjectionY("h_mass_fail", startLS, endLS, "e");

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
    performCount(eff, errl, errh, h_mass_pass, h_mass_fail, cpass, cfail, effType, etaRegion, iBin, lumi, format);
    //cout << effType.Data() << ": " << eff << " + " << errh << " - " << errl << endl;
  }else{
    performFit(eff, errl, errh, chi2pass, chi2fail, h_mass_pass, h_mass_fail, sigModPass, bkgModPass, sigModFail, bkgModFail, cpass, cfail, effType, etaRegion, iBin, lumi, format);
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
	const TString purwDir,
	const Int_t   meanPU 
){
  cout << "Creating histogram templates... "; cout.flush();

  TFile *infile    = new TFile(mcfilename);
  TTree *eventTree = (TTree*)infile->Get("Events");

  Int_t modPU = meanPU;
  if(modPU < minPU) modPU = minPU;
  else if(modPU > maxPU) modPU = maxPU;

  TFile *f_rw = TFile::Open(Form("%s/ReweightHistogram/pileup_rw_%d.root", purwDir.Data(), modPU), "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("h_rw");

  Float_t mass, pt, eta, phi;
  Double_t weight, wgt;
  Int_t q;
  UInt_t npv, npu, pass, runNum, lumiSec, evtNum;

  eventTree->SetBranchAddress("mass",   &mass);
  eventTree->SetBranchAddress("pt",     &pt);
  eventTree->SetBranchAddress("eta",    &eta);
  eventTree->SetBranchAddress("phi",    &phi);
  eventTree->SetBranchAddress("weight", &weight);
  eventTree->SetBranchAddress("q",      &q);
  eventTree->SetBranchAddress("npv",    &npv);
  eventTree->SetBranchAddress("npu",    &npu);
  eventTree->SetBranchAddress("pass",   &pass);
  eventTree->SetBranchAddress("runNum", &runNum);
  eventTree->SetBranchAddress("lumiSec",&lumiSec);
  eventTree->SetBranchAddress("evtNum", &evtNum);

  TH1D *h_mass_pass_central = new TH1D("h_mass_pass_central", "", massBin, massLo, massHi);
  TH1D *h_mass_fail_central = new TH1D("h_mass_fail_central", "", massBin, massLo, massHi);
  TH1D *h_mass_pass_forward = new TH1D("h_mass_pass_forward", "", massBin, massLo, massHi);
  TH1D *h_mass_fail_forward = new TH1D("h_mass_fail_forward", "", massBin, massLo, massHi);

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
    eventTree->GetEntry(ientry);

    if(mass < massLo)  continue;
    if(mass > massHi)  continue;
    if(pt   < ptCut)   continue;
    if(fabs(eta) > etaCut) continue;

    wgt = weight;
    //no pileup reweighting:  wgt *= h_rw->GetBinContent(h_rw->FindBin(npu)); 
    if(fabs(eta) < etaBound){
      if(pass) h_mass_pass_central->Fill(mass, wgt);
      else     h_mass_fail_central->Fill(mass, wgt);
    }else{
      if(pass) h_mass_pass_forward->Fill(mass, wgt);
      else     h_mass_fail_forward->Fill(mass, wgt);
    }
  }
  delete infile;
  infile=0, eventTree=0;

  TFile outfile("histTemplates.root", "RECREATE");
  h_mass_pass_central->Write();
  h_mass_fail_central->Write();
  h_mass_pass_forward->Write();
  h_mass_fail_forward->Write();
  outfile.Write();
  outfile.Close();

  f_rw->Close();
  delete f_rw;

  cout << "Done!" << endl;
}
//--------------------------------------------------------------------------------------------------
void performCount(
	Double_t &resEff, Double_t &resErrl, Double_t &resErrh, TH1D *passHist, TH1D *failHist,
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

  sprintf(binlabely, "27 GeV/c < p_{T} < 13000 GeV/c");
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
// 	1: Breit-Wigner convolved with Crystal Ball function
// 	2: MC template convolved with Gaussian
// Background Model:
//  	1: exponential model
//  	2: quadratic model
//  	3: exponential + quadratic model
//  	4: Das (wide gaussian with exponential tails) model
//  	5: Das + decaying exponential model

  TFile *histfile = 0;
  if(sigpass==2 || sigfail==2) {
    histfile = new TFile("histTemplates.root");
    assert(histfile);
  }
  std::vector<double> vBkgPars;
  if(bkgfail== 2 or bkgfail==3){
    vBkgPars = preFit(failHist);
  }

  CSignalModel     *sigPass = 0;
  CBackgroundModel *bkgPass = 0;
  CSignalModel     *sigFail = 0;
  CBackgroundModel *bkgFail = 0;

  Int_t nflpass=0, nflfail=0;

  switch(sigpass) {
    case 1:
      sigPass = new CBreitWignerConvCrystalBall(m,kTRUE);
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
  }

  switch(sigfail) {
    case 1:
      sigFail = new CBreitWignerConvCrystalBall(m,kFALSE);
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
  if(effType == "SIT" or effType == "HLT") {Nsig.setRange(0,2.0*NsigMax);strategy=1;} 
  RooMsgService::instance().setSilentMode(kTRUE);
  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  fitResult = totalPdf.fitTo(*dataCombined,
                             RooFit::PrintEvalErrors(-1),
                             RooFit::PrintLevel(-1),
                             RooFit::Warnings(0),
                             RooFit::Extended(),
                             RooFit::Strategy(strategy), // MINOS STRATEGY
                             RooFit::Minos(RooArgSet(eff)),
                             RooFit::Save());

  if((fabs(eff.getErrorLo())<5e-4) || (eff.getErrorHi()<5e-4))
    fitResult = totalPdf.fitTo(*dataCombined, RooFit::PrintEvalErrors(-1), RooFit::Extended(), RooFit::Strategy(1), RooFit::Save());

  resEff  = eff.getVal();
  resErrl = fabs(eff.getErrorLo());
  resErrh = eff.getErrorHi();

  // Plot tag and passing- and failing- probe mass distribution 
  char bkgpassname[50];
  char bkgfailname[50];
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

  sprintf(bkgpassname, "backgroundPass_%d",etaRegion);
  sprintf(bkgfailname, "backgroundFail_%d",etaRegion);

  if(!etaRegion) sprintf(binlabelx, "0.0 < |#eta| < 0.9");
  else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");

  sprintf(binlabely, "27 GeV/c < p_{T} < 13000 GeV/c");
  sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",resEff,resErrl,resErrh);
  sprintf(lumitext,"%.1f pb^{-1}  at  #sqrt{s} = 13 TeV",lumi);
//  sprintf(lumitext,"  %.0f lumi-sections at  #sqrt{s} = 13 TeV",lumi);
  sprintf(ylabel,"Events / 1 GeV/c^{2}");

  RooPlot *mframePass = m.frame(Bins(massBin));
  dataPass->plotOn(mframePass,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  if(bkgpass>0)
    modelPass->plotOn(mframePass,Components(bkgpassname),LineStyle(kDashed),LineColor(kRed));
  modelPass->plotOn(mframePass);

  RooPlot *mframeFail = m.frame(Bins(massBin));
  dataFail->plotOn(mframeFail,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
  modelFail->plotOn(mframeFail,Components(bkgfailname),LineStyle(kDashed),LineColor(kRed));
  modelFail->plotOn(mframeFail);

  double a = NsigPass.getVal(), b = NbkgPass.getVal();
  sprintf(pname,"%s_%s_pass_%d", effType.Data(), etaRegion ? "forward" : "central", iBin);
  sprintf(yield,"%u Events",(Int_t)passHist->GetEntries());
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigPass.getVal(),NsigPass.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframePass->chiSquare(nflpass));resChi2Pass = mframePass->chiSquare(nflpass);
  if(bkgpass>0)
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgPass.getVal(),NbkgPass.getPropagatedError(*fitResult));
  CPlot plotPass(pname,mframePass,"Passing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotPass.AddTextBox(binlabelx,0.21,0.78,0.51,0.83,0,kBlack,-1);
  plotPass.AddTextBox(binlabely,0.21,0.73,0.51,0.78,0,kBlack,-1);
  plotPass.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
  plotPass.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
  if(bkgpass>0) {
    plotPass.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  } else {
    plotPass.AddTextBox(0.70,0.73,0.94,0.83,0,kBlack,-1,1,nsigstr);
    plotPass.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,0);
  }
  plotPass.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotPass.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);

  plotPass.Draw(cpass,kTRUE,format);

  double f = NsigFail.getVal(), d = NbkgFail.getVal();
  sprintf(pname,"%s_%s_fail_%d", effType.Data(), etaRegion ? "forward" : "central", iBin);
  sprintf(yield,"%u Events",(Int_t)failHist->GetEntries());
  sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",NsigFail.getVal(),NsigFail.getPropagatedError(*fitResult));
  sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",NbkgFail.getVal(),NbkgFail.getPropagatedError(*fitResult));
  sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeFail->chiSquare(nflfail));resChi2Fail = mframeFail->chiSquare(nflfail);
  CPlot plotFail(pname,mframeFail,"Failing probes","tag-probe mass [GeV/c^{2}]",ylabel);
  plotFail.AddTextBox(binlabelx,0.21,0.75,0.51,0.80,0,kBlack,-1);
  plotFail.AddTextBox(binlabely,0.21,0.70,0.51,0.75,0,kBlack,-1);
  plotFail.AddTextBox(yield,0.21,0.66,0.51,0.70,0,kBlack,-1);
  plotFail.AddTextBox(effstr,0.70,0.85,0.94,0.90,0,kBlack,-1);
  plotFail.AddTextBox(0.70,0.68,0.94,0.83,0,kBlack,-1,2,nsigstr,nbkgstr);
  plotFail.AddTextBox(chi2str,0.70,0.62,0.94,0.67,0,kBlack,-1);
  plotFail.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
  plotFail.AddTextBox(lumitext,0.62,0.92,0.94,0.99,0,kBlack,-1);
  plotFail.Draw(cfail,kTRUE,format);

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

