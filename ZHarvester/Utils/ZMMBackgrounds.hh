#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooUniform.h"
#include "RooGenericPdf.h"
#include "RooGaussDoubleSidedExp.h"
#include "RooCMSShape.h"

class CBackgroundModel
{
public:
    CBackgroundModel():model(0){}
    virtual ~CBackgroundModel() { delete model; }
    RooAbsPdf *model;
    RooRealVar *alpha,*beta,*gamma; // parameters of resolution function need to be defined here to be able to put an external constraint to it
    virtual void freeze_all_parameters(Bool_t freeze=kTRUE) = 0;
    virtual void Reset() = 0;
    virtual void Print() = 0;
};

class CUniform : public CBackgroundModel
{
public:
    CUniform(RooRealVar &m, const Bool_t pass, const Int_t ibin);
    ~CUniform(){};
    void freeze_all_parameters(Bool_t freeze=kTRUE){};
    void Reset(){};
    void Print(){};
};

class CExponential : public CBackgroundModel
{
public:
    CExponential(RooRealVar &m, const Bool_t pass, const Int_t ibin);
    ~CExponential();
    RooRealVar *t1;
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset();
    void Print();

};

class CLinear : public CBackgroundModel
{
public:
    CLinear(RooRealVar &m, const Bool_t pass, const Int_t ibin);
    ~CLinear();
    RooRealVar *a0;
    RooRealVar *a1;
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset();
    void Print();
};

class CQuadratic : public CBackgroundModel
{
public:
  CQuadratic(RooRealVar &m, const Bool_t pass, const Int_t ibin,
      const float p0=0, const float e0=0, const float p1=0, const float e1=0, const float p2=0, const float e2=0);
  ~CQuadratic();
  RooRealVar *a0;
  RooRealVar *a1;
  RooRealVar *a2;
  void freeze_all_parameters(Bool_t freeze=kTRUE);
  void Reset(){};
  void Print(){};

};

class CQuadPlusExp : public CBackgroundModel
{
public:
    CQuadPlusExp(RooRealVar &m, const Bool_t pass, const Int_t ibin,
      const float p0=0, const float e0=0, const float p1=0, const float e1=0, const float p2=0, const float e2=0);
    RooRealVar *a0, *a1, *a2, *t1, *frac;
    RooAbsPdf *exp, *quad;
    ~CQuadPlusExp();
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset(){};
    void Print(){};

};

class CDas: public CBackgroundModel
{
public:
    CDas(RooRealVar &m, const Bool_t pass, const Int_t ibin);
    RooRealVar *mean,*sigma,*kLo,*kHi;
    ~CDas();
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset();
    void Print();

};

class CDasPlusExp: public CBackgroundModel
{
public:
    CDasPlusExp(RooRealVar &m, const Bool_t pass, const Int_t ibin);
    RooRealVar *mean,*sigma,*kLo,*kHi,*t1, *frac;
    RooGaussDoubleSidedExp *dd;
    RooExponential *exp1;
    ~CDasPlusExp();
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset();
    void Print();

};

class CRooCMSShape: public CBackgroundModel
{
public:
    CRooCMSShape(RooRealVar &m, const Bool_t pass, const Int_t ibin, const Double_t massLo=20, const Double_t massHi=80);
    RooRealVar *peak;
    ~CRooCMSShape();
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset();
    void Print();

};

// class CRooChebyshev: public CBackgroundModel
// {
// public:
//     CRooChebyshev(RooRealVar &m, const Bool_t pass, const Int_t ibin);
//     RooRealVar *peak;
//     ~CRooChebyshev();
//     void freeze_all_parameters(Bool_t freeze=kTRUE);
//     void Reset();
//     void Print();

// };

class CQCD: public CBackgroundModel
{
public:
    CQCD(RooRealVar &m, TH1D* hist, const Bool_t pass, const Int_t ibin, Int_t intOrder=1);
    TH1D        *inHist;
    RooDataHist *dataHist;
    RooHistPdf  *histPdf;
    ~CQCD();
    void freeze_all_parameters(Bool_t freeze=kTRUE){};
    void Reset(){};
    void Print(){};

};

class CQCDPlusTT: public CBackgroundModel
{
public:
    CQCDPlusTT(RooRealVar &m, TH1D* histQCD, TH1D* histTT, const Bool_t pass, const Int_t ibin, Int_t intOrder=1);
    RooRealVar  *frac;
    TH1D        *inHistQCD;
    TH1D        *inHistTT;
    RooDataHist *dataHistQCD;
    RooHistPdf  *histPdfQCD;
    RooDataHist *dataHistTT;
    RooHistPdf  *histPdfTT;
    ~CQCDPlusTT();
    void freeze_all_parameters(Bool_t freeze=kTRUE);
    void Reset(){};
    void Print(){};

};

//--------------------------------------------------------------------------------------------------
CUniform::CUniform(RooRealVar &m, const Bool_t pass, const Int_t ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"background%s",name);
  model = new RooUniform(vname,"uniform",m);
}

//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass, const Int_t ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_t1%s",name);    t1 = new RooRealVar(vname,"bkg_t1",-0.1,-1.,0.);

  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,"exp",m,*t1);
}
void CExponential::Reset(){
    std::cout<<"CExponential::Reset() --- "<<std::endl;
    t1->setVal(-0.1);
}

void CExponential::Print(){
    std::cout<<"CExponential::Print() --- "<<std::endl;
    t1->Print();
}

void CExponential::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    t1->setConstant(freeze);
}

CExponential::~CExponential()
{
  delete t1;
}

//--------------------------------------------------------------------------------------------------
CLinear::CLinear(RooRealVar &m, const Bool_t pass, const Int_t ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_a0%s",name);  a0 = new RooRealVar(vname,"bkg_a0",0.,-100000.,100000.);
  sprintf(vname,"bkg_a1%s",name);  a1 = new RooRealVar(vname,"bkg_a1",0.,-1000.,1000.);

  sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,"linear","@1 + @2*@0",RooArgList(m,*a0,*a1));
}

void CLinear::Reset(){
    std::cout<<"CLinear::Reset() --- "<<std::endl;
    a0->setVal(0.);
    a1->setVal(0.);
}

void CLinear::Print(){
    std::cout<<"CLinear::Print() --- "<<std::endl;
    a0->Print();
    a1->Print();
}

void CLinear::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    a0->setConstant(freeze);
    a1->setConstant(freeze);
}

CLinear::~CLinear()
{
  delete a0;
  delete a1;
}


//--------------------------------------------------------------------------------------------------
CQuadratic::CQuadratic(RooRealVar &m, const Bool_t pass, const Int_t ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    sprintf(vname,"bkg_a0%s",name);  a0 = new RooRealVar(vname,"bkg_a0",p0,p0-e0,p0+e0);
    sprintf(vname,"bkg_a1%s",name);  a1 = new RooRealVar(vname,"bkg_a1",p1,p1-e1,p1+e1);
    float upper = p2+e2 > 0. ? 0. : p2+e2;
    std::cout<<p2<<", "<<p2-e2<<", "<<upper<<std::endl;
    sprintf(vname,"bkg_a2%s",name);  a2 = new RooRealVar(vname,"bkg_a2",p2,p2-e2,upper);
  }else{
    sprintf(vname,"bkg_a0%s",name);  a0 = new RooRealVar(vname,"bkg_a0",0.,-100000.,100000.);
    sprintf(vname,"bkg_a1%s",name);  a1 = new RooRealVar(vname,"bkg_a1",0.,-1000.,1000.);
    sprintf(vname,"bkg_a2%s",name);  a2 = new RooRealVar(vname,"bkg_a2",0.,-10.,10.);
  }

  sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,"quad","@1 + @2*@0 + @3*@0*@0",RooArgList(m,*a0,*a1,*a2));
}
void CQuadratic::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    a0->setConstant(freeze);
    a1->setConstant(freeze);
    a2->setConstant(freeze);
}

CQuadratic::~CQuadratic()
{
  delete a0;
  delete a1;
  delete a2;
}

//--------------------------------------------------------------------------------------------------
CQuadPlusExp::CQuadPlusExp(RooRealVar &m, const Bool_t pass, const Int_t ibin,
    const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{
    char name[10];
    if(pass) sprintf(name,"%s_%d","Pass",ibin);
    else     sprintf(name,"%s_%d","Fail",ibin);
    char vname[50];

    if((p0!=0.)||(p1!=0.)||(p2!=0.)){
        sprintf(vname,"bkg_a0%s",name);  a0 = new RooRealVar(vname, "bkg_a0",p0,p0-e0,p0+e0);
        sprintf(vname,"bkg_a1%s",name);  a1 = new RooRealVar(vname, "bkg_a1",p1,p1-e1,p1+e1);
        float upper = p2+e2 > 0. ? 0. : p2+e2;
        sprintf(vname,"bkg_a2%s",name);  a2 = new RooRealVar(vname, "bkg_a2",p2,p2-e2,upper);
    }else{
        sprintf(vname,"bkg_a0%s",name);  a0 = new RooRealVar(vname, "bkg_a0",0.,-100000.,100000.);
        sprintf(vname,"bkg_a1%s",name);  a1 = new RooRealVar(vname, "bkg_a1",0.,-1000.,1000.);
        sprintf(vname,"bkg_a2%s",name);  a2 = new RooRealVar(vname, "bkg_a2",0.,-10.,10.);
    }

    sprintf(vname,"bkg_pdf1%s",name);      quad = new RooGenericPdf(vname, "bkg_quad","@1 + @2*@0 + @3*@0*@0",RooArgList(m,*a0,*a1,*a2));
    sprintf(vname,"bkg_t1%s",name);        t1   = new RooRealVar(vname, "bkg_t1",-0.1,-1.,0.);
    sprintf(vname,"bkg_pdf2%s",name);      exp  = new RooExponential(vname, "bkg_exp",m,*t1);
    sprintf(vname,"bkg_frac%s",name);      frac   = new RooRealVar(vname, "bkg_frac", .95, 0.,1.);

    sprintf(vname,"background%s",name);
    model = new RooAddPdf(vname, "quad + exp", RooArgList(*quad,*exp), RooArgList(*frac));
}

void CQuadPlusExp::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    a0->setConstant(freeze);
    a1->setConstant(freeze);
    a2->setConstant(freeze);
    t1->setConstant(freeze);
    frac->setConstant(freeze);
}

CQuadPlusExp::~CQuadPlusExp()
{
  delete a0;
  delete a1;
  delete a2;
  delete t1;
  delete quad;
  delete exp;
  delete frac;
}

//--------------------------------------------------------------------------------------------------
CDas::CDas(RooRealVar &m, const Bool_t pass, const Int_t ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_mean%s",name);    mean   = new RooRealVar(vname, "bkg_mean", 90, 30, 200);
  sprintf(vname,"bkg_sigma%s",name);   sigma  = new RooRealVar(vname, "bkg_sigma", 30, 10, 60);
  sprintf(vname,"bkg_kLo%s",name);     kLo    = new RooRealVar(vname, "bkg_kLo", 1.5, .02, 10);
  sprintf(vname,"bkg_kHi%s",name);     kHi    = new RooRealVar(vname, "bkg_kHi", 1.5, .02, 10);

  sprintf(vname,"background%s",name);
  model  = new RooGaussDoubleSidedExp(vname,"das",m,*mean,*sigma,*kLo,*kHi);
}

void CDas::Reset(){
    std::cout<<"CDas::Reset() --- "<<std::endl;
    mean->setVal(90.);
    sigma->setVal(12);
    kLo->setVal(1.5);
    kHi->setVal(1.5);
}

void CDas::Print(){
    std::cout<<"CDas::Print() --- "<<std::endl;
    mean->Print();
    sigma->Print();
    kLo->Print();
    kHi->Print();
}

void CDas::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    mean->setConstant(freeze);
    sigma->setConstant(freeze);
    kLo->setConstant(freeze);
    kHi->setConstant(freeze);
}

CDas::~CDas()
{
  delete mean;
  delete sigma;
  delete kLo;
  delete kHi;
}

//--------------------------------------------------------------------------------------------------
CDasPlusExp::CDasPlusExp(RooRealVar &m, const Bool_t pass, const Int_t ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_mean%s",name);     mean   = new RooRealVar(vname, "bkg_mean" ,   90,    30,200);
  sprintf(vname,"bkg_sigma%s",name);    sigma  = new RooRealVar(vname, "bkg_sigma",   30,    10, 60);
  sprintf(vname,"bkg_kLo%s",name);      kLo    = new RooRealVar(vname, "bkg_kLo"  ,  1.5,   .02, 10);
  sprintf(vname,"bkg_kHi%s",name);      kHi    = new RooRealVar(vname, "bkg_kHi"  ,  1.5,   .02, 10);
  sprintf(vname,"bkg_pdf1%s",name);     dd     = new RooGaussDoubleSidedExp(vname, "bkg_das",m,*mean,*sigma,*kLo,*kHi);
  sprintf(vname,"bkg_t1%s",name);       t1     = new RooRealVar(vname, "bkg_t1"  ,-0.1,-1.,0.);
  sprintf(vname,"bkg_frac%s",name);     frac   = new RooRealVar(vname, "bkg_frac", .95, 0.,1.);
  sprintf(vname,"bkg_pdf2%s",name);     exp1   = new RooExponential(vname, "bkg_exp1",m,*t1);

  sprintf(vname,"background%s",name);
  model = new RooAddPdf(vname, "das + exp", RooArgList(*dd,*exp1), RooArgList(*frac));
}

void CDasPlusExp::Reset(){
    std::cout<<"CDasPlusExp::Reset() --- "<<std::endl;
    mean->setVal(90.);
    sigma->setVal(12);
    kLo->setVal(1.5);
    kHi->setVal(1.5);
    t1->setVal(-0.1);
    frac->setVal(0.95);
}

void CDasPlusExp::Print(){
    std::cout<<"CDasPlusExp::Print() --- "<<std::endl;
    mean->Print();
    sigma->Print();
    kLo->Print();
    kHi->Print();
    t1->Print();
    frac->Print();
}

void CDasPlusExp::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    mean->setConstant(freeze);
    sigma->setConstant(freeze);
    kLo->setConstant(freeze);
    kHi->setConstant(freeze);
    t1->setConstant(freeze);
    frac->setConstant(freeze);
}

CDasPlusExp::~CDasPlusExp()
{
  delete mean;
  delete sigma;
  delete kLo;
  delete kHi;
  delete t1;
  delete frac;
  delete exp1;
  delete dd;
}

//--------------------------------------------------------------------------------------------------
CRooCMSShape::CRooCMSShape(RooRealVar &m,
    const Bool_t pass, const Int_t ibin,
    const Double_t massLo, const Double_t massHi)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_alpha%s",name);   alpha = new RooRealVar(vname, "bkg_alpha", 50., massLo, massHi);
  sprintf(vname,"bkg_beta%s",name);    beta  = new RooRealVar(vname, "bkg_beta",  0.05, 0.0, 1.0);
  sprintf(vname,"bkg_gamma%s",name);   gamma = new RooRealVar(vname, "bkg_gamma", 0.05, 0.0, 1.0);
  sprintf(vname,"bkg_peak%s",name);    peak  = new RooRealVar(vname, "bkg_peak",  91.1876);
  sprintf(vname,"background%s",name);

  model  = new RooCMSShape(vname, "RooCMSShape", m, *alpha, *beta, *gamma, *peak);

}

void CRooCMSShape::Reset(){
    std::cout<<"CRooCMSShape::Reset() --- "<<std::endl;
    alpha->setVal(50.);
    beta->setVal(0.05);
    gamma->setVal(0.05);
}

void CRooCMSShape::Print(){
    std::cout<<"CRooCMSShape::Print() --- "<<std::endl;
    alpha->Print();
    beta->Print();
    gamma->Print();
    peak->Print();
}

void CRooCMSShape::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    alpha->setConstant(freeze);
    beta->setConstant(freeze);
    gamma->setConstant(freeze);
    peak->setConstant(freeze);
}

CRooCMSShape::~CRooCMSShape()
{
  delete alpha;
  delete beta;
  delete gamma;
  delete peak;
}

// //--------------------------------------------------------------------------------------------------
// CRooChebyshev::CRooChebyshev(RooRealVar &m,
//     const Bool_t pass, const Int_t ibin)
// {
//   char name[10];
//   if(pass) sprintf(name,"%s_%d","Pass",ibin);
//   else     sprintf(name,"%s_%d","Fail",ibin);
//   char vname[50];

//   sprintf(vname,"bkg_alpha%s",name);   p0 = new RooRealVar(vname, "bkg_p0", 0.0, -1.0, 1.0);
//   sprintf(vname,"bkg_beta%s",name);    p1  = new RooRealVar(vname, "bkg_p1", 0.0, -1.0, 1.0);
//   sprintf(vname,"bkg_gamma%s",name);   p2 = new RooRealVar(vname, "bkg_p2", 0.0, -1.0, 1.0);
//   sprintf(vname,"bkg_peak%s",name);    p3  = new RooRealVar(vname, "bkg_p3", 0.0, -1.0, 1.0);
//   sprintf(vname,"background%s",name);

//   model  = new CRooChebyshev(vname, "RooCMSShape", m, *alpha, *beta, *gamma, *peak);

// }

// void CRooChebyshev::Reset(){
//     std::cout<<"CRooChebyshev::Reset() --- "<<std::endl;
//     alpha->setVal(50.);
//     beta->setVal(0.05);
//     gamma->setVal(0.05);
// }

// void CRooChebyshev::Print(){
//     std::cout<<"CRooChebyshev::Print() --- "<<std::endl;
//     alpha->Print();
//     beta->Print();
//     gamma->Print();
//     peak->Print();
// }

// void CRooChebyshev::freeze_all_parameters(Bool_t freeze){
//     std::cout<<"freeze all parameters"<<std::endl;
//     alpha->setConstant(freeze);
//     beta->setConstant(freeze);
//     gamma->setConstant(freeze);
//     peak->setConstant(freeze);
// }

// CRooChebyshev::~CRooChebyshev()
// {
//   delete alpha;
//   delete beta;
//   delete gamma;
//   delete peak;
// }

//--------------------------------------------------------------------------------------------------
CQCD::CQCD(RooRealVar &m, TH1D* hist, const Bool_t pass, const Int_t ibin, Int_t intOrder)
{
  char name[10];
  if(pass) sprintf(name,"%s_%i","Pass",ibin);
  else     sprintf(name,"%s_%i","Fail",ibin);

  char vname[50];

  sprintf(vname,"bkg_inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);

  sprintf(vname,"bkg_hist%s",name);  dataHist = new RooDataHist(vname, "bkg_dataHist",RooArgSet(m),inHist);
  sprintf(vname,"bkg_pdf%s",name);   histPdf  = new RooHistPdf(vname, "bkg_histPdf",m,*dataHist,intOrder);

  sprintf(vname,"background%s",name);
  model = new RooHistPdf(vname,"QCD",m,*dataHist,intOrder);
}

CQCD::~CQCD()
{
  delete inHist;
  delete dataHist;
  delete histPdf;
}

//--------------------------------------------------------------------------------------------------
CQCDPlusTT::CQCDPlusTT(RooRealVar &m, TH1D* histQCD, TH1D* histTT, const Bool_t pass, const Int_t ibin, Int_t intOrder)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];

  sprintf(vname,"bkg_inHist_QCD_%s",histQCD->GetName());  inHistQCD = (TH1D*)histQCD->Clone(vname);
  sprintf(vname,"bkg_inHist_TT_%s",histTT->GetName());    inHistTT = (TH1D*)histTT->Clone(vname);

  sprintf(vname,"bkg_hist_QCD%s",name);   dataHistQCD = new RooDataHist(vname,"bkg_hist_QCD",RooArgSet(m),inHistQCD);
  sprintf(vname,"bkg_pdf1%s",name);       histPdfQCD  = new RooHistPdf(vname,"bkg_pdf_QCD",m,*dataHistQCD,intOrder);
  sprintf(vname,"bkg_hist_TT%s",name);    dataHistTT = new RooDataHist(vname,"bkg_hist_TT",RooArgSet(m),inHistTT);
  sprintf(vname,"bkg_pdf2%s",name);       histPdfTT  = new RooHistPdf(vname,"bkg_pdf_TT",m,*dataHistTT,intOrder);
  sprintf(vname,"bkg_frac%s",name);       frac = new RooRealVar(vname,"bkg_frac", .5, 0.,1.);

  sprintf(vname,"background%s",name);
  model = new RooAddPdf(vname,"QCD + ttbar MC",RooArgList(*histPdfQCD,*histPdfTT),RooArgList(*frac));

}

void CQCDPlusTT::freeze_all_parameters(Bool_t freeze){
    std::cout<<"freeze all parameters"<<std::endl;
    frac->setConstant(freeze);
}

CQCDPlusTT::~CQCDPlusTT()
{
  delete inHistQCD;
  delete dataHistQCD;
  delete histPdfQCD;
  delete inHistTT;
  delete dataHistTT;
  delete histPdfTT;
  delete frac;
}
