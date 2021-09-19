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
    const TString effType,
    TH1D *hPV=0
);

//--------------------------------------------------------------------------------------------------
// generate template for extraction of Z yield in {BB, BE, EE} region
void generateTemplate_ZYield(
	const TString mcfilename,
	TH1D          *hPV
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
    const TString effType,
	TH1D          *hPV
){

    TFile *outfile = TFile::Open(outputDir+"/histTemplates_"+effType+".root","CREATE");
    if(!outfile){
        cout << "Use existing template "<< endl;
        return;
    }
    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get(effType);
    TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

    if(hPV)
        hPV->Divide(hPVtemplate);

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    Int_t npv;
    Bool_t pass;

    eventTree->SetBranchAddress("mass",           &mass);
    eventTree->SetBranchAddress("ptTag",          &ptTag);
    eventTree->SetBranchAddress("ptProbe",        &ptProbe);
    eventTree->SetBranchAddress("etaTag",         &etaTag);
    eventTree->SetBranchAddress("etaProbe",       &etaProbe);
    eventTree->SetBranchAddress("nPV",            &npv);
    eventTree->SetBranchAddress("pass",           &pass);
    eventTree->SetBranchAddress("eventWeight",    &wgt);

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

    // set negative bin entries to 0
    for(int i=1; i <= massBin; i++){
        if(h_mass_pass_central->GetBinContent(i) < 0.)
            h_mass_pass_central->SetBinContent(i, 0.);

        if(h_mass_fail_central->GetBinContent(i) < 0.)
            h_mass_fail_central->SetBinContent(i, 0.);

        if(h_mass_pass_forward->GetBinContent(i) < 0.)
            h_mass_pass_forward->SetBinContent(i, 0.);

        if(h_mass_fail_forward->GetBinContent(i) < 0.)
            h_mass_fail_forward->SetBinContent(i, 0.);

    }
    outfile->cd();
    h_mass_pass_central->Write();
    h_mass_fail_central->Write();
    h_mass_pass_forward->Write();
    h_mass_fail_forward->Write();
    outfile->Write();
    outfile->Close();

    infile->Close();
    delete infile;

    cout << "Done!" << endl;
}

//--------------------------------------------------------------------------------------------------
void generateTemplate_ZYield(
	const TString mcfilename,
	TH1D          *hPV
){
    TFile *outfile = TFile::Open(outputDir+"/histTemplates_HLT.root","CREATE");
    if(!outfile){
        cout << "Use existing template "<< endl;
        return;
    }
    cout << "Creating histogram templates... "; cout.flush();

    TFile *infile    = new TFile(mcfilename);
    TTree *eventTree = (TTree*)infile->Get("HLT");
    TH1D *hPVtemplate = (TH1D*)infile->Get("hPV");

    if(hPV)
        hPV->Divide(hPVtemplate);

    Double_t mass, ptTag, etaTag, ptProbe, etaProbe;
    Double_t wgt;
    Int_t npv;
    Bool_t pass;

    eventTree->SetBranchAddress("mass",           &mass);
    eventTree->SetBranchAddress("ptTag",          &ptTag);
    eventTree->SetBranchAddress("ptProbe",        &ptProbe);
    eventTree->SetBranchAddress("etaTag",         &etaTag);
    eventTree->SetBranchAddress("etaProbe",       &etaProbe);
    eventTree->SetBranchAddress("nPV",            &npv);
    eventTree->SetBranchAddress("pass",           &pass);
    eventTree->SetBranchAddress("eventWeight",    &wgt);

    TH1D *h_mass_pass_BB = new TH1D("h_mass_pass_BB", "", massBin, massLo, massHi);
    TH1D *h_mass_pass_BE = new TH1D("h_mass_pass_BE", "", massBin, massLo, massHi);
    TH1D *h_mass_pass_EE = new TH1D("h_mass_pass_EE", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_BB = new TH1D("h_mass_fail_BB", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_BE = new TH1D("h_mass_fail_BE", "", massBin, massLo, massHi);
    TH1D *h_mass_fail_EE = new TH1D("h_mass_fail_EE", "", massBin, massLo, massHi);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        eventTree->GetEntry(ientry);

        if(mass < massLo)  continue;
        if(mass > massHi)  continue;
        if(ptTag   < ptCutTag)   continue;
        if(ptProbe   < ptCutProbe)   continue;
        if(fabs(etaTag) > etaCutTag) continue;
        if(fabs(etaProbe) > etaCutProbe) continue;

        if(hPV)
            wgt *= hPV->GetBinContent(hPV->FindBin(npv));

        if(fabs(etaProbe) < etaBound && fabs(etaTag) < etaBound){
            if(pass) h_mass_pass_BB->Fill(mass, wgt);
            else     h_mass_fail_BB->Fill(mass, wgt);
        }
        else if(fabs(etaProbe) >= etaBound && fabs(etaTag) >= etaBound){
            if(pass) h_mass_pass_EE->Fill(mass, wgt);
            else     h_mass_fail_EE->Fill(mass, wgt);
        }
        else{
            if(pass) h_mass_pass_BE->Fill(mass, wgt);
            else     h_mass_fail_BE->Fill(mass, wgt);
        }
    }

    // set negative bin entries to 0
    for(int i=1; i <= massBin+1; i++){
        if(h_mass_pass_BB->GetBinContent(i) < 0.)
            h_mass_pass_BB->SetBinContent(i, 0.);

        if(h_mass_pass_BE->GetBinContent(i) < 0.)
            h_mass_pass_BE->SetBinContent(i, 0.);

        if(h_mass_pass_EE->GetBinContent(i) < 0.)
            h_mass_pass_EE->SetBinContent(i, 0.);

        if(h_mass_fail_BB->GetBinContent(i) < 0.)
            h_mass_fail_BB->SetBinContent(i, 0.);

        if(h_mass_fail_BE->GetBinContent(i) < 0.)
            h_mass_fail_BE->SetBinContent(i, 0.);

        if(h_mass_fail_EE->GetBinContent(i) < 0.)
            h_mass_fail_EE->SetBinContent(i, 0.);
    }

    outfile->cd();
    h_mass_pass_BB->Write();
    h_mass_pass_BE->Write();
    h_mass_pass_EE->Write();
    h_mass_fail_BB->Write();
    h_mass_fail_BE->Write();
    h_mass_fail_EE->Write();
    outfile->Write();
    outfile->Close();

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
    RooWorkspace *w,
    T &Nsig,
    RooRealVar &Nbkg,
    const Int_t nEntries,
    const Int_t iBin,
    const TString effType="yield",
    const Bool_t etaRegion=kTRUE,
    const Bool_t passRegion=kTRUE,
    const Double_t yMax=0.
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
    char effstr[100] = "";

    RooRealVar *eff = (RooRealVar*)w->var("eff");

    if(effType == "BB" || effType == "BE" || effType == "EE"){
        sprintf(pname,"yield_%s_%s_%d", effType.Data(), passRegion ? "2HLT" : "1HLT", iBin);
        sprintf(ctitle,"%s %s %d ", passRegion ? "2HLT" : "1HLT", effType.Data(), iBin);
        if(effType == "BB"){
            sprintf(binlabelx, "|#eta| < %.1f",etaBound);
            if(eff==0)
                eff = (RooRealVar*)w->var("effB");
        }
        if(effType == "BE"){
            sprintf(binlabelx, "|#eta| < 2.4");
        }
        if(effType == "EE"){
            sprintf(binlabelx, "%.1f < |#eta| < 2.4",etaBound);
            if(eff==0)
                eff = (RooRealVar*)w->var("effE");
        }
    }
    else{
        // try to find efficiency in fitResult
        if(effType=="Sta"){
            eff = (RooRealVar*)w->var("effSta");
        }
        else if(effType=="Trk"){
            eff = (RooRealVar*)w->var("effTrk");
        }
        else if(effType=="Glo"){
            RooFormulaVar *effGlo = (RooFormulaVar*)w->arg("effGlo");
            if(effGlo!=0){
                eff = new RooRealVar("eff","eff",effGlo->getVal());
                eff->setError(effGlo->getPropagatedError(*fitResult));
            }
        }

        sprintf(pname,"%s_%s_%s_%d", effType.Data(), etaRegion ? "forward" : "central", passRegion ? "pass" : "fail", iBin);
        sprintf(ctitle,"%s %s (%d) ", effType.Data(), passRegion ? "pass" : "fail", iBin);

        if(!etaRegion) sprintf(binlabelx, "|#eta| < 0.9");
        else           sprintf(binlabelx, "0.9 < |#eta| < 2.4");
    }
    if(eff != 0){
        sprintf(effstr,"#varepsilon = %.4f_{ -%.4f}^{ +%.4f}",eff->getVal(),std::abs(eff->getErrorLo()),eff->getErrorHi());
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
// perform count for tag and probe efficiency
std::vector<float> performCount(
	TH1D *passHist,
    TH1D *failHist,
	const TString effType,
    const Bool_t etaRegion,
    const Int_t iBin
){

    const Double_t npass  = passHist->Integral();
    const Double_t ntotal = failHist->Integral() + npass;
    const Double_t resEff  = (ntotal>0) ? npass/ntotal : 0;

    // Calculate the boundaries for the frequentist Clopper-Pearson interval
    const Double_t resErrl = resEff - TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kFALSE);
    const Double_t resErrh = TEfficiency::ClopperPearson((UInt_t)ntotal, (UInt_t)npass, 0.68269, kTRUE) - resEff;

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
    plotPass.Draw(cpass,kTRUE,"png");

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

    plotFail.Draw(cfail,kTRUE,"png");

    std::vector<float> resultEff = {};

    resultEff.push_back(resEff);
    resultEff.push_back(resErrl);
    resultEff.push_back(resErrh);

    return resultEff;
}

//--------------------------------------------------------------------------------------------------
// perform fit on one histogram for yield extraxtion
void getZyield(
    TH1D *h_yield,
    const Int_t    iBin,         // Label of measurement in currect run
    const Int_t    sigMod,
    const Int_t    bkgMod,
    const TString  etaRegion="",    // {BB, BE or EE}
    const Bool_t   passRegion=true,
    const TString  mcfilename="",
    TH1D*          _hQCD=0,
    TH1D*          hPV=0
){
    std::cout<<">>> Do fit in "<< etaRegion <<" for Z yield extraction"<<std::endl;

    RooRealVar m("m", "mass", massLo, massHi);
    m.setBins(10000);

    CSignalModel     *sigModel = 0;
    CBackgroundModel *bkgModel = 0;

    Int_t nfl=0;

    TFile *histfile = 0;
    TH1D *h=0;
    if(sigMod==2) {
        generateTemplate_ZYield(mcfilename, hPV);
        histfile = new TFile(outputDir+"/histTemplates_HLT.root");
        if(passRegion)
            h = (TH1D*)histfile->Get("h_mass_pass_"+etaRegion);
        else
            h = (TH1D*)histfile->Get("h_mass_fail_"+etaRegion);

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
    RooFitResult *best_fitResult=0;
    Int_t best_status = 0;
    Double_t best_minNll = 0;
    int i = 0;
    do {
        std::cout<<">>> Fit with strategy "<<i<<std::endl;
        // different fit strategies:
        // Strategy 0: just fit fail and pass region simultaneously in full range
        // Strategy 1: first fit fail pdf in sideband range,
        //    then in complete fail range, then fit fail and pass region simultaneously in full range
        // Strategy 2: Same as 1 but with wider central region instead of sideband regions

        if(i == 1){
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
        if(i == 2){
            // fit signal shape to central region only
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

        fitResult = modelPdf.fitTo(*data,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(),
            RooFit::Strategy(1), // MINOS STRATEGY
            //RooFit::Minos(RooArgSet()),
            RooFit::Save());

        if(fitResult->minNll() < best_minNll && fitResult->status() >= 2){
            best_minNll = fitResult->minNll();
            best_fitResult = fitResult;
            best_status = fitResult->status();
        }

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ strategy = " << i << std::endl;
        std::cout<<"------ -log(L)  = " << fitResult->minNll() <<std::endl;
        std::cout<<"------ status   = " << fitResult->status() <<std::endl;
        std::cout<<"---------------------------------------" <<std::endl;

        i++;

    } while(i < 3 && best_status != 4);

    RooFormulaVar purity("purity","Nsig/(Nsig+Nbkg)",RooArgList(Nsig,Nbkg));

    TFile *fFit = new TFile(
        outputDir+"/workspace_yield_"+etaRegion+"_"+iBin+".root",
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(purity);
    w->import(*data);
    w->import(modelPdf);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    const Double_t chi2 = make_plot(nfl, m, data, modelPdf, sigModel, bkgModel,
        best_fitResult, w, Nsig, Nbkg, h_yield->Integral(), iBin, etaRegion, 0, passRegion);

    w->factory(Form("chi2pass[%f]",chi2));

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2 <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    delete sigModel;
    delete bkgModel;
    delete data;
}

//--------------------------------------------------------------------------------------------------
// perform fit on two histograms for tag and probe efficiency
void calculateDataEfficiency(
        TH1D          *passHist,            // histogram with passing probes
        TH1D          *failHist,            // histogram with failing probes
		const Int_t   iBin,                 // Label of measurement number in currect run
		const TString effType,              // "HLT" or "Sel" or "Glo" or "Sta" or "Trk"
        const Bool_t  etaRegion,            // Barrel (0) or Endcap (1)
		const Int_t   sigpass,              // signal model for PASS sample
		const Int_t   bkgpass,              // background model for PASS sample
		const Int_t   sigfail,              // signal model for FAIL sample
		const Int_t   bkgfail,              // background model for FAIL sample
        TH1D          *hPV=0,
        const TString mcfilename="",        // ROOT file containing MC events to generate templates from
        const TString bkgQCDFilename="",    // ROOT file containing bkg template
        const TString bkgTTFilename=""
){

    std::cout<<">>> Do fit in "<< (etaRegion ? "B" : "E") <<" for "<<effType<<std::endl;

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("Pass",1);
    sample.defineType("Fail",2);

    RooAbsData* dataPass=0;
    RooAbsData* dataFail=0;
    RooAbsData* dataCombined=0;

    dataPass     = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
    dataFail     = new RooDataHist("dataFail","dataFail",RooArgSet(m),failHist);
    dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("Pass",*((RooDataHist*)dataPass)),
        RooFit::Import("Fail",*((RooDataHist*)dataFail)));

    TFile *histfile = 0;
    if(sigpass==2 || sigfail==2) {
        generateTemplate(mcfilename, effType, hPV);
        histfile = new TFile(outputDir+"/histTemplates_"+effType+".root");
        assert(histfile);
    }
    std::vector<double> vBkgPars;
    if(bkgfail== 2 or bkgfail==3){
        vBkgPars = preFit(failHist);
    }
    TFile *histbkgQCDfile = 0;
    // TFile *histbkgTTfile = 0;
    if(bkgfail== 7 or bkgpass==7){
        histbkgQCDfile = new TFile(bkgQCDFilename);
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
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(
            Form("h_mass_%s_pass_%s", effType.Data(), etaRegion ? "forward" : "central"));
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
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get(
            Form("h_mass_%s_fail_%s", effType.Data(), etaRegion ? "forward" : "central"));
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

    TH1D *hPass=0;
    if(sigpass==2) {
        hPass = (TH1D*)histfile->Get(Form("h_mass_pass_%s", etaRegion ? "forward" : "central"));
        hPass->SetDirectory(0);
    }

    nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 0, hPass);
    nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 0, hbkgQCDPass);

    TH1D *hFail=0;
    if(sigfail==2) {
        hFail = (TH1D*)histfile->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
        hFail->SetDirectory(0);
    }

    nflfail += set_signal_model(sigfail, sigFail, m, kFALSE, 0, hFail);
    nflfail += set_background_model(bkgfail, bkgFail, m, kFALSE, 0, hbkgQCDFail, &vBkgPars);

    Double_t NsigMax     = passHist->Integral()+failHist->Integral();
    Double_t NbkgFailMax = failHist->Integral();
    Double_t NbkgPassMax = passHist->Integral();
    RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0,1.5*NsigMax);
    RooRealVar eff("eff","Efficiency",0.98,0.,1.0);
    RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);
    if(bkgpass==0)
        NbkgPass.setVal(0);
    RooRealVar NbkgFail("NbkgFail","Background count in FAIL sample",0.1*NbkgFailMax, 0.0, NbkgFailMax);

    RooFormulaVar NsigPass("NsigPass","eff*Nsig",RooArgList(eff,Nsig));
    RooFormulaVar NsigFail("NsigFail","(1.0-eff)*Nsig",RooArgList(eff,Nsig));

    RooAddPdf modelPass("modelPass","Model for PASS sample",
        (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
        (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
    RooAddPdf modelFail("modelFail","Model for FAIL sample",
        RooArgList(*(sigFail->model),*(bkgFail->model)),RooArgList(NsigFail,NbkgFail));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPass,"Pass");
    totalPdf.addPdf(modelFail,"Fail");

    Int_t strategy = 2;
    if(effType == "Sel" or effType == "HLT") {
        Nsig.setRange(0,2.0*NsigMax);
        strategy=1;
    }

    RooMsgService::instance().setSilentMode(kTRUE);

    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Int_t best_status = 0;
    Double_t best_minNll = 0;
    int i = 0;
    do {
        std::cout<<">>> Fit with strategy "<<i<<std::endl;
        // different fit strategies:
        // Strategy 0: just fit fail and pass region simultaneously in full range
        // Strategy 1: first fit fail pdf in sideband range,
        //    then in complete fail range, then fit fail and pass region simultaneously in full range
        // Strategy 2: Same as 1 but with wider sideband regions

        if(i > 0) {

            if(i==1) {
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

        if(fitResult->minNll() < best_minNll && fitResult->status() >= 2){
            best_minNll = fitResult->minNll();
            best_fitResult = fitResult;
            best_status = fitResult->status();
        }

        std::cout<<"---------------------------------------" <<std::endl;
        std::cout<<"------ strategy = " << i << std::endl;
        std::cout<<"------ -log(L)  = " << fitResult->minNll() <<std::endl;
        std::cout<<"------ status   = " << fitResult->status() <<std::endl;
        std::cout<<"---------------------------------------" <<std::endl;

        i++;

    } while(i < 3 && best_status != 4);

    TFile *fFit = new TFile(
        Form(outputDir+"/workspace_%s_%s_%i.root", effType.Data(), etaRegion ? "E" : "B", iBin),
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    bkgFail->Print();
    bkgPass->Print();

    const Double_t chi2pass = make_plot(nflpass, m, dataPass, modelPass,
        sigPass, bkgPass, best_fitResult, w, NsigPass, NbkgPass, passHist->Integral(), iBin,
        effType.Data(), etaRegion, kTRUE);

    Double_t yMax = 0.;
    if(effType == "Glo")
        yMax = 1.5*failHist->GetMaximum();

    const Double_t chi2fail = make_plot(nflfail, m, dataFail, modelFail,
        sigFail, bkgFail, best_fitResult, w, NsigFail, NbkgFail, failHist->Integral(), iBin,
        effType.Data(), etaRegion, kFALSE, yMax);

    w->factory(Form("chi2pass[%f]",chi2pass));
    w->factory(Form("chi2fail[%f]",chi2fail));

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
    std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

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
// perform simultaneous fit in three histograms for global muon efficiency
void calculateGloEfficiency(
        TH1D          *passHist,            // hist with global muon passing probes
        TH1D          *failHistTrk,         // hist with inner track failing probes
        TH1D          *failHistSta,         // hist with outer track failing probes
		const Int_t   iBin,                 // Label of measurement number in currect run
        const Bool_t  etaRegion,            // Barrel (1) or Endcap (0)
		const Int_t   sigpass,              // signal extraction method for PASS sample
		const Int_t   bkgpass,              // background model for PASS sample
		const Int_t   sigfail,              // signal extraction method for FAIL sample
		const Int_t   bkgfail,              // background model for FAIL sample
        TH1D          *hPV=0,
        const TString mcfilename="",        // ROOT file containing DY MC events to generate templates from
        const TString bkgQCDFilename="",    // ROOT file containing bkg template
        const TString bkgTTFilename=""
){

    std::cout<<">>> Do fit in "<< (etaRegion ? "B" : "E") <<" for "<<std::endl;

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("Pass",1);
    sample.defineType("FailTrk",2);
    sample.defineType("FailSta",3);

    RooAbsData *dataPass=0;
    RooAbsData *dataFailTrk=0;
    RooAbsData *dataFailSta=0;
    RooAbsData *dataCombined=0;

    dataPass     = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
    dataFailTrk  = new RooDataHist("dataFailTrk","dataFailTrk",RooArgSet(m),failHistTrk);
    dataFailSta  = new RooDataHist("dataFailSta","dataFailSta",RooArgSet(m),failHistSta);
    dataCombined = new RooDataHist("dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("Pass",*((RooDataHist*)dataPass)),
        RooFit::Import("FailTrk",*((RooDataHist*)dataFailTrk)),
        RooFit::Import("FailSta",*((RooDataHist*)dataFailSta))
        );

    TFile *histfileTrk = 0;
    TFile *histfileSta = 0;
    if(sigpass==2 || sigfail==2) {
        // Generate histogram templates from MC (if necessary)
        generateTemplate(mcfilename, "Sta", hPV);
        generateTemplate(mcfilename, "Trk", hPV);
        histfileTrk = new TFile(outputDir+"/histTemplates_Trk.root");
        histfileSta = new TFile(outputDir+"/histTemplates_Sta.root");
        assert(histfileTrk);
        assert(histfileSta);
    }

    CSignalModel     *sigPass = 0;
    CBackgroundModel *bkgPass = 0;
    CSignalModel     *sigFailTrk = 0;
    CBackgroundModel *bkgFailTrk = 0;
    CSignalModel     *sigFailSta = 0;
    CBackgroundModel *bkgFailSta = 0;

    Int_t nflpass=0, nflfailTrk=0, nflfailSta=0;

    TH1D *h=0;
    if(sigpass==2) {
        // Sta pass == Trk pass =: Glo
        h = (TH1D*)histfileSta->Get(Form("h_mass_pass_%s", etaRegion ? "forward" : "central"));
        h->SetDirectory(0);
    }

    nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 0, h);
    nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 0);

    TH1D *hSta = 0;
    TH1D *hTrk = 0;
    if(sigfail==2) {
        hTrk = (TH1D*)histfileTrk->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
        hSta = (TH1D*)histfileSta->Get(Form("h_mass_fail_%s", etaRegion ? "forward" : "central"));
        hTrk->SetDirectory(0);
        hSta->SetDirectory(0);
    }

    nflfailTrk += set_signal_model(sigfail, sigFailTrk, m, kFALSE, 0, hTrk);
    nflfailTrk += set_background_model(bkgfail, bkgFailTrk, m, kFALSE, 0);
    nflfailSta += set_signal_model(sigfail, sigFailSta, m, kFALSE, 1, hSta);
    nflfailSta += set_background_model(bkgfail, bkgFailSta, m, kFALSE, 1);

    Double_t NsigMax = passHist->Integral() + failHistTrk->Integral() + failHistSta->Integral();
    Double_t NbkgFailTrkMax = failHistTrk->Integral();
    Double_t NbkgFailStaMax = failHistSta->Integral();
    Double_t NbkgPassMax = passHist->Integral();

    RooRealVar Nsig("Nsig","Signal Yield",NsigMax,0,1.5*NsigMax);
    RooRealVar effTrk("effTrk","Tracking Efficiency",0.98,0.,1.0);
    RooRealVar effSta("effSta","Standalone Efficiency",0.98,0.,1.0);

    RooRealVar NbkgPass("NbkgPass","Background count in PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);

    if(bkgpass==0)
        NbkgPass.setVal(0);
    RooRealVar NbkgFailTrk("NbkgFailTrk", "Background count in Trk FAIL sample",
        0.8*NbkgFailTrkMax, 0.0, NbkgFailTrkMax);
    RooRealVar NbkgFailSta("NbkgFailSta", "Background count in Sta FAIL sample",
        0.9*NbkgFailStaMax, 0.0, NbkgFailStaMax);

    RooFormulaVar NsigPass("NsigPass","effTrk*effSta*Nsig",RooArgList(effTrk, effSta, Nsig));
    RooFormulaVar NsigFailTrk("NsigFailTrk","effSta*(1.0-effTrk)*Nsig",RooArgList(effSta,effTrk,Nsig));
    RooFormulaVar NsigFailSta("NsigFailSta","effTrk*(1.0-effSta)*Nsig",RooArgList(effTrk,effSta,Nsig));

    RooAddPdf modelPass("modelPass","Model for PASS sample",
        (bkgpass>0) ? RooArgList(*(sigPass->model),*(bkgPass->model)) :  RooArgList(*(sigPass->model)),
        (bkgpass>0) ? RooArgList(NsigPass,NbkgPass) : RooArgList(NsigPass));
    RooAddPdf modelFailTrk("modelFailTrk","Model for Trk FAIL sample",
        RooArgList(*(sigFailTrk->model),*(bkgFailTrk->model)),RooArgList(NsigFailTrk, NbkgFailTrk));
    RooAddPdf modelFailSta("modelFailSta","Model for Sta FAIL sample",
        RooArgList(*(sigFailSta->model),*(bkgFailSta->model)),RooArgList(NsigFailSta, NbkgFailSta));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPass,"Pass");
    totalPdf.addPdf(modelFailTrk,"FailTrk");
    totalPdf.addPdf(modelFailSta,"FailSta");

    RooMsgService::instance().setSilentMode(kTRUE);

    Int_t strategy = 2;
    RooFitResult *fitResult=0;
    RooFitResult *best_fitResult=0;
    Int_t best_status = 0;
    Double_t best_minNll = 0;
    int i = 0;
    do {
        std::cout<<">>> Fit with strategy "<<i<<std::endl;
        // different fit strategies:
        // Strategy 0: just fit fail and pass region simultaneously in full range
        // Strategy 1: first fit fail pdf in sideband range,
        //    then in complete fail range, then fit fail and pass region simultaneously in full range
        // Strategy 2: Same as 1 but with wider sideband regions

        if(i > 0) {

            if(i==1) {
                m.setRange("rangeLow", massLo, 81);
                m.setRange("rangeHigh", 101, massHi);
            }
            else{
                m.setRange("rangeLow", massLo, 76);
                m.setRange("rangeHigh", 106, massHi);
            }

            std::cout<<">>> Prefit Trk"<<std::endl;
            // Tracking
            effSta.setConstant(kTRUE);
            // fit background pdf in Fail sideband region
            fitResult = bkgFailTrk->model->fitTo(*dataFailTrk,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                RooFit::Minos(RooArgSet(effTrk)),
                RooFit::Save());

            // fit total pdf in Fail
            fitResult = modelFailTrk.fitTo(*dataFailTrk,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                RooFit::Minos(RooArgSet(effTrk)),
                RooFit::Save());

            effSta.setConstant(kFALSE);

            std::cout<<">>> Prefit Sta"<<std::endl;

            effTrk.setConstant(kTRUE);

            fitResult = bkgFailSta->model->fitTo(*dataFailSta,
                RooFit::Range("rangeLow,rangeHigh"),
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                RooFit::Minos(RooArgSet(effSta)),
                RooFit::Save());

            fitResult = modelFailSta.fitTo(*dataFailSta,
                RooFit::PrintEvalErrors(-1),
                RooFit::PrintLevel(-1),
                RooFit::Warnings(0),
                RooFit::Extended(1),
                RooFit::Strategy(strategy), // MINUIT STRATEGY
                RooFit::Minos(RooArgSet(effSta)),
                RooFit::Save());

            effTrk.setConstant(kFALSE);
        }

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy), // MINUIT STRATEGY
            RooFit::Minos(RooArgSet(effTrk, effSta)),
            RooFit::Save());

        fitResult = totalPdf.fitTo(*dataCombined,
            RooFit::PrintEvalErrors(-1),
            RooFit::PrintLevel(-1),
            RooFit::Warnings(0),
            RooFit::Extended(1),
            RooFit::Strategy(strategy), // MINUIT STRATEGY
            RooFit::Minos(RooArgSet(effTrk, effSta)),
            RooFit::Save());

        if(fitResult->minNll() < best_minNll && fitResult->status() >= 2){
            best_minNll = fitResult->minNll();
            best_fitResult = fitResult;
            best_status = fitResult->status();
        }

        std::cout<<"---------------------------------------"<<std::endl;
        std::cout<<"------ strategy = " << i <<std::endl;
        std::cout<<"------ -log(L)  = " << fitResult->minNll() <<std::endl;
        std::cout<<"------ status   = " << fitResult->status() <<std::endl;
        std::cout<<"---------------------------------------"<<std::endl;

        i++;

    } while(i < 3 && best_status != 4);

    TFile *fFit = new TFile(
        Form(outputDir+"/workspace_Glo_%s_%i.root", etaRegion ? "E" : "B", iBin),
        "RECREATE");

    // calculate global muon efficiency with error propagation
    RooFormulaVar effGlo("effGlo","effTrk*effSta",RooArgList(effTrk, effSta));

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(effGlo);
    w->import(*dataCombined);
    w->import(totalPdf);
    w->Write();

    best_fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    const Double_t chi2pass = make_plot(nflpass, m, dataPass, modelPass,
        sigPass, bkgPass, best_fitResult, w, NsigPass, NbkgPass,
        passHist->Integral(), iBin, "Glo", etaRegion, kTRUE);

    Double_t yMax = 1.5*failHistTrk->GetMaximum();

    const Double_t chi2failTrk = make_plot(nflfailTrk, m, dataFailTrk, modelFailTrk,
        sigFailTrk, bkgFailTrk, best_fitResult, w, NsigFailTrk, NbkgFailTrk,
        failHistTrk->Integral(), iBin, "Trk", etaRegion, kFALSE);

    yMax = 1.5*failHistSta->GetMaximum();

    const Double_t chi2failSta = make_plot(nflfailSta, m, dataFailSta, modelFailSta,
        sigFailSta, bkgFailSta, best_fitResult, w, NsigFailSta, NbkgFailSta,
        failHistSta->Integral(), iBin, "Sta", etaRegion, kFALSE);

    w->factory(Form("chi2pass[%f]",chi2pass));
    w->factory(Form("chi2failSta[%f]",chi2failSta));
    w->factory(Form("chi2failTrk[%f]",chi2failTrk));

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
    std::cout<<"------ chi2failSta = " << chi2failSta <<std::endl;
    std::cout<<"------ chi2failTrk = " << chi2failTrk <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    delete dataCombined;
    delete dataPass;
    delete dataFailTrk;
    delete dataFailSta;
    delete sigPass;
    delete bkgPass;
    delete sigFailTrk;
    delete bkgFailTrk;
    delete sigFailSta;
    delete bkgFailSta;
    delete histfileSta;
    delete histfileTrk;
    // delete histbkgTTfile;
}

//--------------------------------------------------------------------------------------------------
// perform fit in 2 region (1HLT and 2HLT) to extract the Z yield and efficiency together
void calculateHLTEfficiencyAndYield(
    TH1D          *passHist,            // histogram with events where both muon pass HLT
    TH1D          *failHist,            // histogram with events where one muon pass HLT
	const Int_t   iBin,                 // Label of measurement number in currect run
    const TString etaRegion,            // {BB, BE, EE}
	const Int_t   sigpass,              // signal model for PASS sample
	const Int_t   bkgpass,              // background model for PASS sample
	const Int_t   sigfail,              // signal model for FAIL sample
	const Int_t   bkgfail,              // background model for FAIL sample
    const Double_t corr = 1.,           // correlation factors for second muon
    TH1D          *hPV=0,
    const TString mcfilename="",        // ROOT file containing MC events to generate templates from
    const TString bkgQCDFilename="",    // ROOT file containing bkg template
    const TString bkgTTFilename=""
){

    RooRealVar m("m","mass",massLo,massHi);
    m.setBins(10000);

    RooCategory sample("sample","");
    sample.defineType("Pass",1);
    sample.defineType("Fail",2);

    RooAbsData *dataPass = new RooDataHist("dataPass","dataPass",RooArgSet(m),passHist);
    RooAbsData *dataFail = new RooDataHist("dataFail","dataFail",RooArgSet(m),failHist);

    RooAbsData *dataCombined = new RooDataHist(
        "dataCombined","dataCombined",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("Pass",*((RooDataHist*)dataPass)),
        RooFit::Import("Fail",*((RooDataHist*)dataFail)));

    TFile *histfile = 0;
    if(sigpass==2 || sigfail==2) {
        generateTemplate_ZYield(mcfilename, hPV);
        histfile = new TFile(outputDir+"/histTemplates_HLT.root");
        assert(histfile);
    }
    std::vector<double> vBkgPars;
    if(bkgfail== 2 or bkgfail==3){
        vBkgPars = preFit(failHist);
    }
    TFile *histbkgQCDfile = 0;
    if(bkgfail== 7 or bkgpass==7){
        histbkgQCDfile = new TFile(bkgQCDFilename);
        assert(histbkgQCDfile);
    }
    TH1D *hbkgQCDPass = 0;
    TH1D *hbkgQCDFail = 0;
    if(bkgpass == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get("h_mass_HLT_pass_"+ etaRegion);
        hbkgQCDPass = hbkg->ProjectionY();
        assert(hbkgQCDPass);
    }
    if(bkgfail == 7){
        TH2D* hbkg = (TH2D*)histbkgQCDfile->Get("h_mass_HLT_fail_"+ etaRegion);
        hbkgQCDFail = hbkg->ProjectionY();
        assert(hbkgQCDFail);
    }

    CSignalModel     *sigPass = 0;
    CBackgroundModel *bkgPass = 0;
    CSignalModel     *sigFail = 0;
    CBackgroundModel *bkgFail = 0;

    Int_t nflpass=0, nflfail=0;

    TH1D *hPass=0;
    if(sigpass==2) {
        // set signal templates of fail histograms
        hPass = (TH1D*)histfile->Get("h_mass_pass_"+ etaRegion);
        hPass->SetDirectory(0);
    }

    nflpass += set_signal_model(sigpass, sigPass, m, kTRUE, 0, hPass);
    nflpass += set_background_model(bkgpass, bkgPass, m, kTRUE, 0, hbkgQCDPass);

    TH1D *hFail=0;
    if(sigfail==2) {
        // set signal templates of pass histograms
        hFail = (TH1D*)histfile->Get("h_mass_fail_"+ etaRegion);
        hFail->SetDirectory(0);
    }

    nflfail += set_signal_model(sigfail, sigFail, m, kFALSE, 0, hFail);
    nflfail += set_background_model(bkgfail, bkgFail, m, kFALSE, 0, hbkgQCDFail);

    Double_t NsigMax     = passHist->Integral() + failHist->Integral();
    Double_t NbkgFailMax = failHist->Integral();
    Double_t NbkgPassMax = passHist->Integral();

    RooRealVar eff("eff","Efficiency barrel",0.95,0.,1.0);

    RooRealVar Nsig("Nsig","Signal Yield ",NsigMax,0,1.5*NsigMax);
    RooRealVar c("c", "Correlation factor ", corr);

    RooFormulaVar NsigPass("NsigPass","eff*eff*Nsig*c",RooArgList(eff,Nsig,c));
    RooFormulaVar NsigFail("NsigFail","2*eff*(1.0-c*eff)*Nsig",RooArgList(eff,Nsig,c));

    RooRealVar NbkgPass("NbkgPass","Background count in  PASS sample",0.01*NbkgPassMax, 0.0, NbkgPassMax);
    RooRealVar NbkgFail("NbkgFail","Background count in  FAIL sample",0.05*NbkgFailMax, 0.0, NbkgFailMax);

    RooAddPdf modelPass("modelPass ","Model for  PASS sample",
        RooArgList(*(sigPass->model),*(bkgPass->model)),
        RooArgList(NsigPass, NbkgPass));
    RooAddPdf modelFail("modelFail ","Model for  FAIL sample",
        RooArgList(*(sigFail->model),*(bkgFail->model)),
        RooArgList(NsigFail, NbkgFail));

    RooSimultaneous totalPdf("totalPdf","totalPdf",sample);
    totalPdf.addPdf(modelPass,"Pass");
    totalPdf.addPdf(modelFail,"Fail");

    RooFitResult *fitResult=0;
    Int_t strategy = 2;

    RooMsgService::instance().setSilentMode(kTRUE);

    // fit all regions together
    std::cout<<"--- fit all regions together -- "<<std::endl;
    fitResult = totalPdf.fitTo(*dataCombined,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        // RooFit::Minos(RooArgSet(eff)),
        RooFit::Save());

    fitResult = totalPdf.fitTo(*dataCombined,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(1),
        RooFit::Strategy(strategy), // MINUIT STRATEGY
        RooFit::Minos(RooArgSet(Nsig)),
        RooFit::Save());

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ -log(L)  = " << fitResult->minNll() <<std::endl;
    std::cout<<"------ status   = " << fitResult->status() <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

    bkgFail->Print();
    bkgPass->Print();

    TFile *fFit = new TFile(
        outputDir+"/workspace_"+etaRegion+"_"+iBin+".root",
        "RECREATE");

    // save all information in RooWorkspace
    RooWorkspace* w = new RooWorkspace("workspace","Workspace");
    w->import(*dataCombined);
    w->import(totalPdf);
    w->Write();

    fitResult->Write("fitResult");

    fFit->Write();
    fFit->Close();

    const Double_t chi2pass = make_plot(nflpass, m, dataPass, modelPass,
        sigPass, bkgPass, fitResult, w, NsigPass, NbkgPass, passHist->Integral(),
        iBin, etaRegion, 0, kTRUE);

    const Double_t chi2fail = make_plot(nflfail, m, dataFail, modelFail,
        sigFail, bkgFail, fitResult, w, NsigFail, NbkgFail, failHist->Integral(),
        iBin, etaRegion, 0, kFALSE);

    w->factory(Form("chi2pass[%f]",chi2pass));
    w->factory(Form("chi2fail[%f]",chi2fail));

    std::cout<<"---------------------------------------"<<std::endl;
    std::cout<<"------ chi2pass = " << chi2pass <<std::endl;
    std::cout<<"------ chi2fail = " << chi2fail <<std::endl;
    std::cout<<"---------------------------------------"<<std::endl;

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
    }
    else{
        v = {fq->GetParameter(0), fq->GetParError(0), fq->GetParameter(1), fq->GetParError(1), fq->GetParameter(2), fq->GetParError(2)};
    }

    return v;
}
