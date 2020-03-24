#ifndef ROOFITTER_HH
#define ROOFITTER_HH

#include "RooRealVar.h"
#include "RooCategory.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TTree.h>
#include <TStyle.h>

#include "ZMMSignals.hh"
#include "ZMMBackgrounds.hh"
#include "RooGaussDoubleSidedExp.h"
#include "CPlot.hh"
#include "MitStyleRemix.hh"



class RooFitter
{
public:
    RooFitter(){};
    RooFitter(
        const Double_t  _ptCutTag,
        const Double_t  _ptCutProbe,
        const Int_t     _sigModelType,
        const Int_t     _bkgModelType,
        const TString   _outputDir="",
        const TString   _mcfilename_dy="",
        const TString   _mcfilename_tt="",
        TH1D*           _hPV=0,
        const Float_t   _massLo  = 66.,
        const Float_t   _massHi  = 116.
    );
    ~RooFitter(){}

    void update_sigModel(const TString _mcfilename_dy="", const TString _mcfilename_tt="", TH1D* _hPV=0);

    void fit_backgroundmodel(TH1D* _hYield);
    std::vector<float> fit_simultanious(TH1D* _hYieldOS, TH1D* _hYieldSS, const TString _name="");

private:

    const Float_t ptCutTag = 30.;
    const Float_t ptCutProbe = 30.;

    const Float_t massLo  = 66.;
    const Float_t massHi  = 116.;
    const UInt_t  massBin = 50;

    const Float_t etaCutTag   = 2.4;
    const Float_t etaCutProbe = 2.4;
    const Float_t etaBound    = 0.9;


    const UInt_t sigModelType = 0;
    const UInt_t bkgModelType = 0;

    CSignalModel     *modelSig = 0;
    CBackgroundModel *modelBkg = 0;

    Int_t ndfSig=0;
    Int_t ndfBkg=0;

    RooRealVar m;

    RooCategory sample;

    TH1D* generateTemplate_ZYield(
        const TString mcfilename,
    	TH1D*         hPV
    );

};

#endif
