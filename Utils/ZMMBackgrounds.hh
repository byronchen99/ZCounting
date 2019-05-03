#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooGaussDoubleSidedExp.h"

class CBackgroundModel
{
public:
  CBackgroundModel():model(0){}
  virtual ~CBackgroundModel() { delete model; }
  RooAbsPdf *model;
};

class CExponential : public CBackgroundModel
{
public:
  CExponential(RooRealVar &m, const Bool_t pass, const int ibin);
  ~CExponential();
  RooRealVar *t;
};

class CQuadratic : public CBackgroundModel
{
public:
  CQuadratic(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2);
  ~CQuadratic();
  RooRealVar *a0;
  RooRealVar *a1;
  RooRealVar *a2;
};

class CQuadPlusExp : public CBackgroundModel
{
public:
  CQuadPlusExp(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2);
  RooRealVar *a0, *a1, *a2, *t1;
  RooAbsPdf *exp, *quad;
  ~CQuadPlusExp();
};

class CDas: public CBackgroundModel
{
public:
  CDas(RooRealVar &m, const Bool_t pass, const int ibin);
  RooRealVar *mean,*sigma,*kLo,*kHi;
  ~CDas();
};

class CDasPlusExp: public CBackgroundModel 
{
public:
  CDasPlusExp(RooRealVar &m, const Bool_t pass, const int ibin);
  RooRealVar *mean,*sigma,*kLo,*kHi,*t1, *frac;
  RooGaussDoubleSidedExp *dd;
  RooExponential *exp1;
  ~CDasPlusExp();
};


//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  
  char vname[50];
  
  sprintf(vname,"t%s",name);
  t = new RooRealVar(vname,vname,-0.1,-1.,0.);
      
  sprintf(vname,"background%s",name);
  model = new RooExponential(vname,vname,m,*t);
}

CExponential::~CExponential()
{
  delete t;
}

//--------------------------------------------------------------------------------------------------
CQuadratic::CQuadratic(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{ 
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);

  char a0name[50];
  sprintf(a0name,"a0%s",name); 
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    a0 = new RooRealVar(a0name,a0name,p0,p0-e0,p0+e0);
  }else{
    a0 = new RooRealVar(a0name,a0name,0.,-10.,10.);
  }

  char a1name[50]; 
  sprintf(a1name,"a1%s",name);
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    a1 = new RooRealVar(a1name,a1name,p1,p1-e1,p1+e1);
  }else{
    a1 = new RooRealVar(a1name,a1name,0.,-10.,10.);
  }

  char a2name[50]; 
  sprintf(a2name,"a2%s",name);
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
      float upper = p2+e2 > 0. ? 0. : p2+e2;
      a2 = new RooRealVar(a2name,a2name,p2,p2-e2,upper);
      std::cout<<p2<<", "<<p2-e2<<", "<<upper<<std::endl;
  }else{
    a2 = new RooRealVar(a2name,a2name,0.,0.,10.);
  }
  
  char formula[200];
  sprintf(formula,"(%s+%s*m+%s*m*m)",a0name, a1name,a2name);
  
  char vname[50]; sprintf(vname,"background%s",name);
  model = new RooGenericPdf(vname,vname,formula,RooArgList(m,*a0,*a1,*a2));
}

CQuadratic::~CQuadratic()
{
  delete a0;
  delete a1;
  delete a2;
}

//--------------------------------------------------------------------------------------------------
CQuadPlusExp::CQuadPlusExp(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{
  if((p0!=0.)||(p1!=0.)||(p2!=0.)){
    a0 = new RooRealVar("bkg_a0","bkg_a0",p0,p0-e0,p0+e0);
    a1 = new RooRealVar("bkg_a1","bkg_a1",p1,p1-e1,p1+e1);
    float upper = p2+e2 > 0. ? 0. : p2+e2;
    a2 = new RooRealVar("bkg_a2","bkg_a2",p2,p2-e2,upper);
  }else{
    a0 = new RooRealVar("bkg_a0","bkg_a0",0.,-10.,10.);
    a1 = new RooRealVar("bkg_a1","bkg_a1",0.,-10.,10.);
    a2 = new RooRealVar("bkg_a2","bkg_a2",0.,0.,10.);
  }

  quad = new RooGenericPdf("bkg_quad","bkg_quad","@1 + @2*@0 + @3*@0*@0",RooArgList(m,*a0,*a1,*a2));  
  t1   = new RooRealVar("bkg_t1","bkg_t1",-0.1,-1.,0.);
  exp  = new RooExponential("bkg_exp","bkg_exp",m,*t1);

  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];
  sprintf(vname,"background%s",name);

  model = new RooGenericPdf(vname,vname,"@0+@1",RooArgList(*quad,*exp));
}

CQuadPlusExp::~CQuadPlusExp()
{
  delete a0;
  delete a1;
  delete a2;
  delete t1;
  delete quad;
  delete exp;
}

//--------------------------------------------------------------------------------------------------
CDas::CDas(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];
  sprintf(vname,"background%s",name);

  mean   = new RooRealVar("bkg_mean" , "bkg_mean" ,   90,    30,200);
  sigma  = new RooRealVar("bkg_sigma", "bkg_sigma",   12,    10, 60);
  kLo    = new RooRealVar("bkg_kLo"  , "bkg_kLo"  ,  1.5,   .02, 10);
  kHi    = new RooRealVar("bkg_kHi"  , "bkg_kHi"  ,  1.5,   .02, 10);
  model  = new RooGaussDoubleSidedExp("bkgDas","bkgDas",m,*mean,*sigma,*kLo,*kHi);
}

CDas::~CDas()
{
  delete mean;
  delete sigma;
  delete kLo;
  delete kHi;
}

//--------------------------------------------------------------------------------------------------
CDasPlusExp::CDasPlusExp(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"%s_%d","Pass",ibin);
  else     sprintf(name,"%s_%d","Fail",ibin);
  char vname[50];
  sprintf(vname,"background%s",name);

  mean   = new RooRealVar("bkg_mean" , "bkg_mean" ,   90,    30,200);
  sigma  = new RooRealVar("bkg_sigma", "bkg_sigma",   12,    10, 60);
  kLo    = new RooRealVar("bkg_kLo"  , "bkg_kLo"  ,  1.5,   .02, 10);
  kHi    = new RooRealVar("bkg_kHi"  , "bkg_kHi"  ,  1.5,   .02, 10);
  dd     = new RooGaussDoubleSidedExp("bkgDas","bkgDas",m,*mean,*sigma,*kLo,*kHi);
  t1     = new RooRealVar("bkg_t1"  ,"bkg_t1"  ,-0.20,-.4,.4);
  frac   = new RooRealVar("bkg_frac","bkg_frac", 0.05, 0.,1.);
  exp1   = new RooExponential("bkg_exp1","bkg_exp1",m,*t1);
  model  = new RooAddPdf(vname,vname,RooArgList(*dd,*exp1),RooArgList(*frac));
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
