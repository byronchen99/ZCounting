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



class CSignalModel
{
public:
  CSignalModel():model(0){}
  virtual ~CSignalModel(){ delete model; }
  virtual void Reset() = 0;

  RooAbsPdf *model;
};

class CBreitWignerConvCrystalBall : public CSignalModel
{
public:
  CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass, const int ibin);
  ~CBreitWignerConvCrystalBall();
  void Reset(){};
  RooRealVar     *mass, *width;
  RooBreitWigner *bw;
  RooRealVar     *mean, *sigma, *alpha, *n;
  RooCBShape     *cb;
};

class CMCTemplateConvGaussian : public CSignalModel
{
public:
  CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder=1);
  ~CMCTemplateConvGaussian();
  void Reset();
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
};

class CMCStackConvGaussian : public CSignalModel
{
public:
    CMCStackConvGaussian(RooRealVar &m, TH1D* hist_dy, TH1D* hist_tt, const Bool_t pass, const int ibin, int intOrder=1);
    ~CMCStackConvGaussian();
    void Reset(){};
    RooRealVar  *mean, *sigma, *frac;
    RooGaussian *gaus;
    TH1D        *inHist_dy, *inHist_tt;
    RooDataHist *dataHist_dy, *dataHist_tt;
    RooHistPdf  *histPdf_dy, *histPdf_tt;
    RooAddPdf   *stack;
};

//--------------------------------------------------------------------------------------------------
CBreitWignerConvCrystalBall::CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%i","Pass",ibin);
  else     sprintf(name,"%s_%i","Fail",ibin);
  char vname[50];

  sprintf(vname,"sig_mass%s",name);
  mass = new RooRealVar(vname,vname,91,80,100);
  mass->setVal(91.1876);
  mass->setConstant(kTRUE);

  sprintf(vname,"sid_width%s",name);
  width = new RooRealVar(vname,vname,2.5,0.1,10);
  width->setVal(2.4952);
  width->setConstant(kTRUE);

  sprintf(vname,"sig_bw%s",name);
  bw = new RooBreitWigner(vname,vname,m,*mass,*width);

  if(pass) {
    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"sig_alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"sig_n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  } else {
    sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-10,10);
    sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,1,0.1,5);
    sprintf(vname,"sig_alpha%s",name); alpha = new RooRealVar(vname,vname,5,0,20);
    sprintf(vname,"sig_n%s",name);     n     = new RooRealVar(vname,vname,1,0,10);
  }
//  n->setVal(1.0);
//  n->setConstant(kTRUE);

  sprintf(vname,"sig_cb%s",name);
  cb = new RooCBShape(vname,vname,m,*mean,*sigma,*alpha,*n);

  sprintf(vname,"signal%s",name);
  model = new RooFFTConvPdf(vname,"BW x CB",m,*bw,*cb);
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
CMCTemplateConvGaussian::CMCTemplateConvGaussian(RooRealVar &m, TH1D* hist, const Bool_t pass, const int ibin, int intOrder)
{
  char name[10];
  if(pass) sprintf(name,"%s_%i","Pass",ibin);
  else     sprintf(name,"%s_%i","Fail",ibin);
  char vname[50];

  sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
  sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5);
  sprintf(vname,"sig_gaus%s",name);  gaus  = new RooGaussian(vname,vname,m,*mean,*sigma);

  sprintf(vname,"sig_inHist_%s",hist->GetName());
  inHist = (TH1D*)hist->Clone(vname);

  sprintf(vname,"sig_dataHist%s",name); dataHist = new RooDataHist(vname,vname,RooArgSet(m),inHist);
  sprintf(vname,"sig_histPdf%s",name);  histPdf  = new RooHistPdf(vname,vname,m,*dataHist,intOrder);

  sprintf(vname,"signal%s",name);
  model = new RooFFTConvPdf(vname,"MC x Gaus",m,*histPdf,*gaus);
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
CMCStackConvGaussian::CMCStackConvGaussian(RooRealVar &m, TH1D* hist_dy, TH1D* hist_tt, const Bool_t pass, const int ibin, int intOrder)
{
  char name[10];
  if(pass) sprintf(name,"%s_%i","Pass",ibin);
  else     sprintf(name,"%s_%i","Fail",ibin);
  char vname[50];

  sprintf(vname,"sig_mean%s",name);  mean  = new RooRealVar(vname,vname,0,-2.5,2.5);
  sprintf(vname,"sig_sigma%s",name); sigma = new RooRealVar(vname,vname,2,0,5);
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
  model = new RooFFTConvPdf(vname,"MC x Gaus",m,*stack,*gaus);
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
