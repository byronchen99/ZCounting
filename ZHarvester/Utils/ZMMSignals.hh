#include "TROOT.h"
#include "TH1D.h"

#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"


class CSignalModel
{
protected:
    const double ZMASS = 91.1876;
    const double ZWIDTH = 2.4952;

    const int ipOrder = 3;
    
public:
    CSignalModel():model(0){}
    virtual ~CSignalModel(){ delete model; }
    virtual void Reset() = 0;

    RooRealVar  *mean, *sigma, *alpha; // parameters of resolution function need to be defined here to be able to put an external constraint to it
    RooAbsPdf *model;
};

class CBreitWigner : public CSignalModel
{
public:
    CBreitWigner(RooRealVar &m, const Bool_t pass, const int ibin);
    ~CBreitWigner();
    void Reset();
    RooAbsReal     *mass, *width;
};

class CBreitWignerConvCrystalBall : public CSignalModel
{
public:
    CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass, const int ibin);
    ~CBreitWignerConvCrystalBall();
    void Reset();
    RooAbsReal     *mass, *width;
    RooBreitWigner *bw;
    RooRealVar     *n;
    RooCBShape     *cb;
};

class CBreitWignerConvGaussian: public CSignalModel
{
public:
    CBreitWignerConvGaussian(RooRealVar &m, const Bool_t pass, const int ibin);
    ~CBreitWignerConvGaussian();
    void Reset();
    RooAbsReal     *mass, *width;
    RooBreitWigner *bw;
    RooGaussian *gaus;

};

class CMCTemplate : public CSignalModel
{
public:
    CMCTemplate(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder=3);
    ~CMCTemplate();
    void Reset();
    TH1D        *inHist;
    RooDataHist *dataHist;
    RooHistPdf  *histPdf;
};

class CMCTemplateConvGaussian : public CSignalModel
{
public:
    CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder=3);
    ~CMCTemplateConvGaussian();
    void Reset();
    RooGaussian *gaus;
    TH1D        *inHist;
    RooDataHist *dataHist;
    RooHistPdf  *histPdf;
};

class CMCTemplateConvCrystalBall : public CSignalModel
{
public:
    CMCTemplateConvCrystalBall(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder=3);
    ~CMCTemplateConvCrystalBall();
    void Reset();
    RooRealVar  *n;
    RooCBShape  *cb;
    TH1D        *inHist;
    RooDataHist *dataHist;
    RooHistPdf  *histPdf;
};

class CMCStackConvGaussian : public CSignalModel
{
public:
    CMCStackConvGaussian(RooRealVar &m, TH1D* hist_dy, TH1D* hist_tt, const Bool_t pass, const int ibin, int intOrder=3);
    ~CMCStackConvGaussian();
    void Reset(){};
    RooRealVar  *frac;
    RooGaussian *gaus;
    TH1D        *inHist_dy, *inHist_tt;
    RooDataHist *dataHist_dy, *dataHist_tt;
    RooHistPdf  *histPdf_dy, *histPdf_tt;
    RooAddPdf   *stack;
};

//--------------------------------------------------------------------------------------------------
CBreitWigner::CBreitWigner(RooRealVar &m, const Bool_t pass, const int ibin)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mass%s",name);
    mass = new RooConstVar(vname,vname,ZMASS);

    // mass = new RooRealVar(vname,vname,91,80,100);
    // mass->setVal(ZMASS);
    // mass->setConstant(kTRUE);

    sprintf(vname,"sid_width%s",name);
    width = new RooConstVar(vname,vname,ZWIDTH);

    // width = new RooRealVar(vname,vname,2.5,0.1,10);
    // width->setVal(ZWIDTH);
    // width->setConstant(kTRUE);

    sprintf(vname,"sig_bw%s",name);
    model = new RooBreitWigner(vname,"BW", m,*mass,*width);

}

void CBreitWigner::Reset(){
    std::cout<<"CBreitWigner::Reset() --- "<<std::endl;
}

CBreitWigner::~CBreitWigner()
{
    delete mass;
    delete width;
}

//--------------------------------------------------------------------------------------------------
CBreitWignerConvCrystalBall::CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass, const int ibin)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mass%s",name);
    mass = new RooConstVar(vname,vname,ZMASS);

    // mass = new RooRealVar(vname,vname,91,80,100);
    // mass->setVal(ZMASS);
    // mass->setConstant(kTRUE);

    sprintf(vname,"sid_width%s",name);
    width = new RooConstVar(vname,vname,ZWIDTH);

    // width = new RooRealVar(vname,vname,2.5,0.1,10);
    // width->setVal(ZWIDTH);
    // width->setConstant(kTRUE);

    sprintf(vname,"sig_bw%s",name);
    bw = new RooBreitWigner(vname,vname,m,*mass,*width);

    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0.1, 5);
    sprintf(vname,"sig_alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"sig_n%s",name);     n     = new RooRealVar(vname,vname,5,0.5,10);


    sprintf(vname,"sig_cb%s",name);
    cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);

    sprintf(vname,"signal%s",name);
    model = new RooFFTConvPdf(vname,"BW x CB",m,*bw,*cb, ipOrder);
}

void CBreitWignerConvCrystalBall::Reset(){
    std::cout<<"CBreitWignerConvCrystalBall::Reset() --- "<<std::endl;
    mean->setVal(0);
    sigma->setVal(1);
    alpha->setVal(5);
    n->setVal(1);
}

CBreitWignerConvCrystalBall::~CBreitWignerConvCrystalBall()
{
    delete mass;
    delete width;
    delete bw;
    delete mean;
    delete sigma;
    delete alpha;
    delete n;
    delete cb;
}

//--------------------------------------------------------------------------------------------------
CBreitWignerConvGaussian::CBreitWignerConvGaussian(RooRealVar &m, const Bool_t pass, const int ibin)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mass%s",name);
    mass = new RooConstVar(vname,vname,ZMASS);

    // mass = new RooRealVar(vname,vname,91,80,100);
    // mass->setVal(ZMASS);
    // mass->setConstant(kTRUE);

    sprintf(vname,"sid_width%s",name);
    width = new RooConstVar(vname,vname,ZWIDTH);

    // width = new RooRealVar(vname,vname,2.5,0.1,10);
    // width->setVal(ZWIDTH);
    // width->setConstant(kTRUE);

    sprintf(vname,"sig_bw%s",name);
    bw = new RooBreitWigner(vname,vname,m,*mass,*width);

    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0.1,5);
    sprintf(vname,"sig_gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);

    sprintf(vname,"signal%s",name);
    model = new RooFFTConvPdf(vname,"BW x Gaus",m,*bw,*gaus, ipOrder);
}

void CBreitWignerConvGaussian::Reset(){
    std::cout<<"CBreitWignerConvGaussian::Reset() --- "<<std::endl;
    mean->setVal(0);
    sigma->setVal(2);
}

CBreitWignerConvGaussian::~CBreitWignerConvGaussian()
{
    delete mass;
    delete width;
    delete bw;
    delete mean;
    delete sigma;
    delete gaus;
}

//--------------------------------------------------------------------------------------------------
CMCTemplate::CMCTemplate(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_inHist_%s",hist->GetName());
    inHist = (TH1D*)hist->Clone(vname);

    sprintf(vname,"sig_dataHist%s",name);
    dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
    
    sprintf(vname,"signal%s",name);
    model  = new RooHistPdf(vname,"MC",m,*dataHist,intOrder);

}

void CMCTemplate::Reset(){
    std::cout<<"CMCTemplate::Reset() --- "<<std::endl;

}

CMCTemplate::~CMCTemplate()
{
    delete inHist;
    delete dataHist;
}

//--------------------------------------------------------------------------------------------------
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0.1,5);
    sprintf(vname,"sig_gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);

    sprintf(vname,"sig_inHist_%s",hist->GetName());
    inHist = (TH1D*)hist->Clone(vname);

    sprintf(vname,"sig_dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
    sprintf(vname,"sig_histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);

    sprintf(vname,"signal%s",name);
    model = new RooFFTConvPdf(vname,"MC x Gaus",m,*histPdf,*gaus, ipOrder);
}

void CMCTemplateConvGaussian::Reset(){
    std::cout<<"CMCTemplateConvGaussian::Reset() --- "<<std::endl;
    mean->setVal(0);
    sigma->setVal(2);
}

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
    delete mean;
    delete sigma;
    delete gaus;
    delete inHist;
    delete dataHist;
    delete histPdf;
}

//--------------------------------------------------------------------------------------------------
CMCTemplateConvCrystalBall::CMCTemplateConvCrystalBall(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0.1, 5);
    sprintf(vname,"sig_alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"sig_n%s",name);     n     = new RooRealVar(vname,vname,5,0.5,10);

    sprintf(vname,"sig_cb%s",name);    
    cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);

    sprintf(vname,"sig_inHist_%s",hist->GetName());
    inHist = (TH1D*)hist->Clone(vname);

    sprintf(vname,"sig_dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
    sprintf(vname,"sig_histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);

    sprintf(vname,"signal%s",name);
    model = new RooFFTConvPdf(vname,"MC x CB",m,*histPdf,*cb, ipOrder);
}

void CMCTemplateConvCrystalBall::Reset(){
    std::cout<<"CMCTemplateConvCrystalBall::Reset() --- "<<std::endl;
    mean->setVal(0);
    sigma->setVal(1);
    alpha->setVal(5);
    n->setVal(1);
}

CMCTemplateConvCrystalBall::~CMCTemplateConvCrystalBall()
{
    delete mean;
    delete sigma;
    delete alpha;
    delete n;
    delete cb;
    delete inHist;
    delete dataHist;
    delete histPdf;
}

//--------------------------------------------------------------------------------------------------
CMCStackConvGaussian::CMCStackConvGaussian(RooRealVar &m, TH1D* hist_dy, TH1D* hist_tt, const Bool_t pass, const int ibin, int intOrder)
{
    char name[10];
    if(pass) sprintf(name,"%s_%i","Pass",ibin);
    else     sprintf(name,"%s_%i","Fail",ibin);
    char vname[50];

    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0.1,5);
    sprintf(vname,"sig_gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);

    sprintf(vname,"sig_inHist_dy_%s",hist_dy->GetName()); inHist_dy = (TH1D*)hist_dy->Clone(vname);
    sprintf(vname,"sig_inHist_tt_%s",hist_tt->GetName()); inHist_tt = (TH1D*)hist_tt->Clone(vname);

    sprintf(vname,"sig_dataHist_dy_%s",name); dataHist_dy = new RooDataHist(vname,vname,RooArgSet(m),inHist_dy);
    sprintf(vname,"sig_dataHist_tt_%s",name); dataHist_tt = new RooDataHist(vname,vname,RooArgSet(m),inHist_tt);

    sprintf(vname,"sig_histPdf_dy_%s",name);  histPdf_dy  = new RooHistPdf(vname,vname,m,*dataHist_dy,intOrder);
    sprintf(vname,"sig_histPdf_tt_%s",name);  histPdf_tt  = new RooHistPdf(vname,vname,m,*dataHist_tt,intOrder);

    sprintf(vname,"sig_frac_%s",name);      frac = new RooRealVar(vname, "sig_frac", .95, 0.,1.);
    sprintf(vname,"sig_stack_%s",name);
    stack = new RooAddPdf(vname, "dy + tt", RooArgList(*histPdf_dy,*histPdf_tt), RooArgList(*frac));

    sprintf(vname,"signal%s",name);
    model = new RooFFTConvPdf(vname,"MC x Gaus",m,*stack,*gaus, ipOrder);
}

CMCStackConvGaussian::~CMCStackConvGaussian()
{
    delete mean;
    delete sigma;
    delete gaus;
    delete inHist_dy;
    delete inHist_tt;
    delete dataHist_dy;
    delete dataHist_tt;
    delete histPdf_dy;
    delete histPdf_tt;
    delete frac;
    delete stack;
}
