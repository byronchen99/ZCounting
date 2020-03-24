#include <iostream>                 // standard I/O

#include "RooFitter.hh"


RooFitter::RooFitter(
    const Double_t  _ptCutTag,
    const Double_t  _ptCutProbe,
    const Int_t     _sigModelType,
    const Int_t     _bkgModelType,
    const TString   _outputDir,
    const TString   _mcfilename_dy,
    const TString   _mcfilename_tt,
    TH1D*           _hPV,
    const Float_t   _massLo,
    const Float_t   _massHi
):
    ptCutTag(_ptCutTag),
    ptCutProbe(_ptCutProbe),
    m("m","mass",massLo,massHi),
    sigModelType(_sigModelType),
    bkgModelType(_bkgModelType),
    massLo(_massLo),
    massHi(_massHi)
{
    std::cout<<"init RooFitter"<<std::endl;

    CPlot::sOutDir = _outputDir;
    m.setBins(10000);

    switch(_sigModelType) {
        case 1:
            modelSig = new CBreitWignerConvCrystalBall(m, kFALSE, 0);
            ndfSig += 4; break;
        case 2:{
            TH1D *h = generateTemplate_ZYield(_mcfilename_dy, _hPV);
            assert(h);
            modelSig = new CMCTemplateConvGaussian(m,h,kFALSE,0);
            ndfSig += 2; break;
        }
        case 3:{
            TH1D *h_dy = generateTemplate_ZYield(_mcfilename_dy, _hPV);
            TH1D *h_tt = generateTemplate_ZYield(_mcfilename_tt, _hPV);
            assert(h_dy);
            assert(h_tt);
            modelSig = new CMCStackConvGaussian(m,h_dy,h_tt,kFALSE,0);
            ndfSig += 3; break;
        }
    }
    switch(_bkgModelType) {
        case 1:
            modelBkg = new CExponential(m, kFALSE, 0);
            ndfBkg += 1; break;
        case 2:
            modelBkg = new CQuadratic(m, kFALSE, 0);
            ndfBkg += 3; break;
        case 3:
            modelBkg = new CQuadPlusExp(m, kFALSE, 0);
            ndfBkg += 4; break;
        case 4:
            modelBkg = new CDas(m, kFALSE, 0);
            ndfBkg += 4; break;
        case 5:
            modelBkg = new CDasPlusExp(m, kFALSE, 0);
            ndfBkg += 6; break;
    }
}

//--------------------------------------------------------------------------------------------------
void RooFitter::update_sigModel(
    const TString   _mcfilename_dy,
    const TString   _mcfilename_tt,
    TH1D*           _hPV
){
    switch(sigModelType) {
        case 1:
            modelSig = new CBreitWignerConvCrystalBall(m, kFALSE, 0);
            ndfSig = 4; break;
        case 2:{
            TH1D *h = generateTemplate_ZYield(_mcfilename_dy, _hPV);
            assert(h);
            modelSig = new CMCTemplateConvGaussian(m,h,kFALSE,0);
            ndfSig += 2; break;
        }
        case 3:{
            TH1D *h_dy = generateTemplate_ZYield(_mcfilename_dy, _hPV);
            TH1D *h_tt = generateTemplate_ZYield(_mcfilename_tt, _hPV);
            assert(h_dy);
            assert(h_tt);
            modelSig = new CMCStackConvGaussian(m,h_dy,h_tt,kFALSE,0);
            ndfSig += 3; break;
        }
    }
}

//--------------------------------------------------------------------------------------------------
void RooFitter::fit_backgroundmodel(
    TH1D           *_hYield
){

    RooAbsData *data = new RooDataHist("data","data",RooArgSet(m),_hYield);

    const Double_t yield_forBkg = _hYield->Integral();

    RooRealVar nbkg("nbkg","yield in background fit", yield_forBkg, 0., 1.5*yield_forBkg);
    RooAddPdf *modelBkg_forBkg = new RooAddPdf("modelSSLowPU","Model for SS sample Low PU", RooArgList(*(modelBkg->model)), RooArgList(nbkg));


    RooFitResult *fitResult=0;

    fitResult = modelBkg_forBkg->fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(1),
        //RooFit::Minos(RooArgSet(eff)),
        RooFit::Save());

    fitResult = modelBkg_forBkg->fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(2),
        //RooFit::Minos(RooArgSet(eff)),
        RooFit::Save());

    // --- plot same sign region
    char ylabel[50];
    char yield[50];
    char nbkgstr[50];
    char chi2str[50];

    TCanvas *cSS = MakeCanvas("cSS","cSS",720,540);
    cSS->SetWindowPosition(cSS->GetWindowTopX()+cSS->GetBorderSize()+800,0);

    RooPlot *mframeSS = m.frame(Bins(massBin));
    data->plotOn(mframeSS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelBkg_forBkg->plotOn(mframeSS, Components(*(modelBkg->model)), LineColor(8));

    sprintf(yield,"%u Events",(Int_t)_hYield->GetEntries());
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkg.getVal(),nbkg.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeSS->chiSquare(ndfBkg));

    CPlot plotSS("plot_SS_"+std::to_string(10*bkgModelType+sigModelType)+"_Inclusive",
        mframeSS,"Same Sign - Inclusive Measurement","tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotSS.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotSS.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotSS.AddTextBox(std::to_string((int)ptCutProbe)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotSS.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotSS.AddTextBox(modelBkg->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotSS.AddTextBox(0.70,0.79,0.94,0.90,0,kBlack,-1,2,chi2str, nbkgstr);

    const Double_t yMax = _hYield->GetMaximum() + _hYield->GetBinError(_hYield->GetMaximumBin());
    const Double_t yMin = std::max(0., _hYield->GetMinimum() - _hYield->GetBinError(_hYield->GetMinimumBin()));
    plotSS.SetYRange(yMin, yMax + 0.4*(yMax-yMin));
    plotSS.Draw(cSS,kTRUE,"png");

    delete data;
    delete modelBkg_forBkg;
}

//--------------------------------------------------------------------------------------------------
std::vector<float> RooFitter::fit_simultanious(
    TH1D*           _hYieldOS,
    TH1D*           _hYieldSS,
    const TString   _name
){
    /*
        performs a fit of signal and background region in the histogram "_hYieldOS"
        use histogram "_hYieldSS" to compute transfer factor only
    */

    RooCategory sample("sample","");
    sample.defineType("OS",1);
    sample.defineType("SS",2);

    RooAbsData *dataOS = new RooDataHist("dataOS","opposite sign data",RooArgSet(m),_hYieldOS);
    RooAbsData *dataSS = new RooDataHist("dataSS","same sign data",RooArgSet(m),_hYieldSS);

    RooAbsData *data = new RooDataHist("data","data",RooArgList(m),
        RooFit::Index(sample),
        RooFit::Import("OS",*((RooDataHist*)dataOS)),
        RooFit::Import("SS",*((RooDataHist*)dataSS)));

    const Double_t yield_SS = _hYieldSS->Integral();
    const Double_t yield_OS = _hYieldOS->Integral();

    RooRealVar nbkgSS("nbkgSS","yield in SS background fit", yield_SS, 0., 1.5*yield_SS);
    RooRealVar nbkgOS("nbkgOS","yield in OS background fit", yield_SS, 0., yield_OS);
    RooRealVar nsigOS("nsigOS","yield in OS backsignalground fit", yield_OS, 0., 1.5*yield_OS);

    RooFormulaVar fr("fr","@0/(@0 + @1)",RooArgList(nbkgOS,nsigOS));
    RooFormulaVar tf("tf","@0/@1",RooArgList(nbkgOS,nbkgSS));

    RooAddPdf *modelSS = new RooAddPdf("modelSSLowPU","Model for SS sample Low PU", RooArgList(*(modelBkg->model)), RooArgList(nbkgSS));
    RooAddPdf *modelOS = new RooAddPdf("modelOSLowPU","Model for OS sample Low PU", RooArgList(*(modelSig->model),*(modelBkg->model)), RooArgList(nsigOS, nbkgOS));

    modelBkg->freeze_all_parameters();

    RooFitResult *fitResult=0;

    fitResult = modelOS->fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(1),
        //RooFit::Minos(RooArgSet(eff)),
        RooFit::Save());

    fitResult = modelOS->fitTo(*data,
        RooFit::PrintEvalErrors(-1),
        RooFit::PrintLevel(-1),
        RooFit::Warnings(0),
        RooFit::Extended(),
        RooFit::Strategy(2),
        //RooFit::Minos(RooArgSet(eff)),
        RooFit::Save());


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

    // --- plot opposite signregion
    TCanvas *cOS = MakeCanvas("cOS","cOS",720,540);
    cOS->SetWindowPosition(cOS->GetWindowTopX()+cOS->GetBorderSize()+800,0);

    RooPlot *mframeOS = m.frame(Bins(massBin));
    dataOS->plotOn(mframeOS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelOS->plotOn(mframeOS, Components(*(modelSig->model)), LineColor(9));
    modelOS->plotOn(mframeOS, Components(*(modelBkg->model)), LineColor(8));

    modelOS->plotOn(mframeOS, LineColor(kRed));

    sprintf(yield,"%u Events",(Int_t)_hYieldOS->GetEntries());
    sprintf(nsigstr,"N_{sig} = %.1f #pm %.1f",nsigOS.getVal(),nsigOS.getPropagatedError(*fitResult));
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgOS.getVal(),nbkgOS.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeOS->chiSquare(ndfSig+1));
    sprintf(frstr,"fr = %.3f #pm %.3f",resfr, errfr);
    sprintf(tfstr,"tf = %.3f #pm %.3f",restf, errtf);


    CPlot plotOS("plot_OS_"+std::to_string(10*bkgModelType+sigModelType)+"_"+_name,
        mframeOS, "Opposite Sign - Measurement "+_name,"tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotOS.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotOS.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotOS.AddTextBox(std::to_string((int)ptCutProbe)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotOS.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotOS.AddTextBox(modelSig->model->GetTitle(), 0.21, 0.63, 0.41, 0.67, 0, 9, -1, 12);
    plotOS.AddTextBox(modelBkg->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotOS.AddTextBox(0.70,0.62,0.94,0.90,0,kBlack,-1,5,chi2str, nbkgstr, nsigstr, frstr, tfstr);

    plotOS.Draw(cOS,kTRUE,"png");

    plotOS.SetName("plot_OS_logscale_"+std::to_string(10*bkgModelType+sigModelType)+"_"+_name);

    double ymin = _hYieldOS->GetMinimum();
    if(ymin <= 0.) ymin = 1.;

    plotOS.SetYRange(0.5 * ymin, _hYieldOS->GetMaximum()*2);
    plotOS.SetLogy();
    plotOS.Draw(cOS,kTRUE,"png");


    // --- plot same sign region
    TCanvas *cSS = MakeCanvas("cSS","cSS",720,540);
    cSS->SetWindowPosition(cSS->GetWindowTopX()+cSS->GetBorderSize()+800,0);

    RooPlot *mframeSS = m.frame(Bins(massBin));
    dataSS->plotOn(mframeSS,MarkerStyle(kFullCircle),MarkerSize(0.8),DrawOption("ZP"));
    modelSS->plotOn(mframeSS, Components(*(modelBkg->model)), LineColor(8));

    sprintf(yield,"%u Events",(Int_t)_hYieldSS->GetEntries());
    sprintf(nbkgstr,"N_{bkg} = %.1f #pm %.1f",nbkgSS.getVal(),nbkgSS.getPropagatedError(*fitResult));
    sprintf(chi2str,"#chi^{2}/DOF = %.3f",mframeSS->chiSquare(1));

    CPlot plotSS("plot_SS_"+std::to_string(10*bkgModelType+sigModelType)+"_"+_name,
        mframeSS,"Same Sign - Measurement "+_name,"tag-probe mass [GeV/c^{2}]","Events / 1 GeV/c^{2}");

    plotSS.AddTextBox("CMS Preliminary",0.19,0.83,0.54,0.89,0);
    plotSS.AddTextBox("|#eta| < 2.4",0.21,0.78,0.51,0.83,0,kBlack,-1);
    plotSS.AddTextBox(std::to_string((int)ptCutProbe)+" GeV/c < p_{T} < 13000 GeV/c",0.21,0.73,0.51,0.78,0,kBlack,-1);
    plotSS.AddTextBox(yield,0.21,0.69,0.51,0.73,0,kBlack,-1);
    plotSS.AddTextBox(modelBkg->model->GetTitle(), 0.21, 0.59, 0.41, 0.63, 0, 8, -1, 12);
    plotSS.AddTextBox(0.70,0.79,0.94,0.90,0,kBlack,-1,2,chi2str, nbkgstr);

    plotSS.Draw(cSS,kTRUE,"png");


    delete dataOS;
    delete dataSS;
    delete data;
    delete modelOS;
    delete modelSS;

    std::vector<float> result = {};

    result.push_back(resfr);
    result.push_back(errfr);
    result.push_back(restf);
    result.push_back(errtf);
    result.push_back(mframeOS->chiSquare(ndfSig+1));
    result.push_back(mframeSS->chiSquare(1));

    return result;
}


//--------------------------------------------------------------------------------------------------
TH1D* RooFitter::generateTemplate_ZYield(
	const TString mcfilename,
	TH1D*         hPV
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
    h_mass_zyield->SetDirectory(0);

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

    infile->Close();
    delete infile;

    return h_mass_zyield;
}
