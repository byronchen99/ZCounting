#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"

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

class CExpQuad : public CBackgroundModel
{
public:
  CExpQuad(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2);
  CExponential *exp;
  CQuadratic *quad;
  ~CExpQuad();
};

//--------------------------------------------------------------------------------------------------
CExponential::CExponential(RooRealVar &m, const Bool_t pass, const int ibin)
{
  char name[10];
  if(pass) sprintf(name,"Ex%s_%d","Pass",ibin);
  else     sprintf(name,"Ex%s_%d","Fail",ibin);
  
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
  if(pass) sprintf(name,"Q%s_%d","Pass",ibin);
  else     sprintf(name,"Q%s_%d","Fail",ibin);

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
CExpQuad::CExpQuad(RooRealVar &m, const Bool_t pass, const int ibin, const float p0, const float e0, const float p1, const float e1, const float p2, const float e2)
{
  char name[10];
  if(pass) sprintf(name,"ExQ_%s_%d","Pass",ibin);
  else     sprintf(name,"ExQ_%s_%d","Fail",ibin);
  char vname[50];
  sprintf(vname,"background%s",name);
  
  quad = new CQuadratic(m, pass, ibin, p0, e0, p1, e1, p2, e2);
  exp = new CExponential(m, pass, ibin);
  
  std::cout<<quad->model<<std::endl;
  std::cout<<exp->model<<std::endl;

  model = new RooGenericPdf(vname,vname,"@0+@1",RooArgList(*(quad->model),*(exp->model)));
  std::cout<<model<<std::endl;
}

CExpQuad::~CExpQuad()
{
  delete quad;
  delete exp;
}

