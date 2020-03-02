#include "TROOT.h"
#include "TH1D.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"

class CSignalModel
{
public:
  CSignalModel():model(0){}
  virtual ~CSignalModel(){ delete model; }
  RooAbsPdf *model;
};

class CBreitWignerConvCrystalBall : public CSignalModel
{
public:
  CBreitWignerConvCrystalBall(RooRealVar &m, const Bool_t pass, const int ibin);
  ~CBreitWignerConvCrystalBall();
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
  RooRealVar  *mean, *sigma;
  RooGaussian *gaus;
  TH1D        *inHist;
  RooDataHist *dataHist;
  RooHistPdf  *histPdf;
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

CMCTemplateConvGaussian::~CMCTemplateConvGaussian()
{
  delete mean;
  delete sigma;
  delete gaus;
  delete inHist;
  delete dataHist;
  delete histPdf;
}
