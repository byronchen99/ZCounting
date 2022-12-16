//------------------------------------------------------------------------------
//
// BLUE: A ROOT class implementing the Best Linear Unbiased Estimate method.
//
// Copyright (C) 2012-2019, Richard.Nisius@mpp.mpg.de
// All rights reserved
//
// This file is part of BLUE - Version 2.2.0.
//
// BLUE is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// For the licensing terms see the file COPYING or http://www.gnu.org/licenses.
//
//------------------------------------------------------------------------------
#include "Blue.h"
#include <iostream>                 // standard I/O
#include "TROOT.h"

void LUMI_COMBINATION_1(Int_t Flag = 0){

    gROOT->SetBatch(1);

  //----------------------------------------------------------------------------
  // Flag steers which of the results should be calculated
  //  0: Z count to estimate 2017 lumi
  //  1: Z count to estimate 2018 lumi
  //----------------------------------------------------------------------------

  // The number of estimates, uncertainties and observables
  static const Int_t NumEst =  3;
  static const Int_t NumUnc =  4;
  static const Int_t MaxObs =  2;
  Int_t NumObs =  2;

  // The names
  TString NamEst[NumEst] = {"   P_17", "   P_18", "   Z_17"};
  TString NamUnc[NumUnc] = {" Uncorr", "   Corr", "      Z"};
  TString NamObs[MaxObs] = {"L(2017)", "L(2018)"};

  // Index for which estimates determines which observable
  Int_t IWhichObs[NumEst] = {0, 1, 0};

  // Preset according to Flag
  if(Flag == 0){
    printf("... Combination of 2017, 2018 PHYSICS luminosity with Z counts for 2017");
    printf(" Flag = %2i \n", Flag);
  }else if(Flag == 1){
    printf("... Combination of 2017, 2018 PHYSICS luminosity with Z counts for 2018");
    printf(" Flag = %2i \n", Flag);

    NamEst[2] = "   Z_18";
    IWhichObs[2] = 1;

  }else{
    printf("... : Not implemented Flag = %2i \n",Flag);
    return;
  }

  // Character to set the file names
  char Buffer[100];

  // Estimates 0-6
  // First the F0
  // 0 == PHYSICS luminosity for 2017
  // 1 == PHYSICS luminosity for 2018
  // 2 == Z Luminosity for 2017 from L_Z^2017 = L_PHYSICS^2018 * N_Z^2017 / N_Z^2018
  // OR: 2 == Z Luminosity for 2018 from L_Z^2018 = L_PHYSICS^2017 * N_Z^2018 / N_Z^2017

  // Some toy value for the ratio of Z count such that equal to ratio
  // assumed uncertainty on Z counts
  Double_t uncZ = 0.01;

  // uncertainties in absolute numbers
  static const Int_t LenXEst = NumEst * (NumUnc+1);
  Double_t XEst[LenXEst] = {
    //         0      1      2      3
    //     Uncorr  Corr1  Corr2     Z
    41.48, 0.8296, 0.3733, 0.2489, 0.0,
    59.83, 0.8975, 1.1966, 0.1197, 0.0,
    41.48, 0.6222, 0.8296, 0.0830, uncZ*41.48
  };

  if(Flag == 1){
    XEst[10] = 59.83;
    XEst[11] = 1.1966;
    XEst[12] = 0.5385;
    XEst[13] = 0.3590;
    XEst[14] =  uncZ*59.83;
  }

  // Uncertainty source 1: Uncorrelated PHYSICS luminosity -> but correlated with Z counts
  static const Int_t LenCor = NumEst * NumEst;
  // Z lumi 2017
  Double_t Cor00[LenCor] = {
    +1.00,  0.00,  0.00,
    +0.00,  1.00,  1.00,
    +0.00,  1.00,  1.00
  };
  if(Flag == 1){
      // Z lumi 2018
      Cor00[2] = 1.00;
      Cor00[5] = 0.00;
      Cor00[6] = 1.00;
      Cor00[7] = 0.00;

  }
  // Uncertainty source 2: correlated PHYSICS luminosity 16,17,18 -> correlated with Z counts
  Double_t Cor01[LenCor] = {
    +1.00,  1.00,  1.00,
    +1.00,  1.00,  1.00,
    +1.00,  1.00,  1.00
  };
  // Uncertainty source 2: correlated PHYSICS luminosity 17,18 -> correlated with Z counts
  Double_t Cor02[LenCor] = {
    +1.00,  1.00,  1.00,
    +1.00,  1.00,  1.00,
    +1.00,  1.00,  1.00
  };
  // Uncertainty source 3: Uncertainty on Z count ratio
  Double_t Cor03[LenCor] = {
    +1.00,  0.00,  0.00,
    +0.00,  1.00,  0.00,
    +0.00,  0.00,  1.00
  };

  //-- Local Structures for Blue output
  // TMatrices
  TMatrixD* LocRho    = new TMatrixD(NumEst,NumEst);
  TMatrixD* LocRhoRes = new TMatrixD(NumObs,NumObs);
  TMatrixD* LocWeight = new TMatrixD(NumEst,NumObs);
  //-- End

  // Define formats for Figures and Latex file
  const TString ForVal = "%5.3f";
  const TString ForUnc = ForVal;
  const TString ForWei = "%5.3f";
  const TString ForRho = ForWei;
  const TString ForPul = ForVal;
  const TString ForChi = "%5.3f";
  const TString ForUni = "";

  // Construct Object
  Blue *myBlue = new Blue(NumEst, NumUnc, NumObs, &IWhichObs[0]);
  myBlue->PrintStatus();
  myBlue->SetFormat(ForVal, ForUnc, ForWei, ForRho, ForPul, ForChi, ForUni);

  // Fill names
  myBlue->FillNamEst(&NamEst[0]);
  myBlue->FillNamUnc(&NamUnc[0]);
  myBlue->FillNamObs(&NamObs[0]);

  // Fill estimates
  Int_t ind = 0;
  for(Int_t i = 0; i<NumEst; i++){
    myBlue->FillEst(i,&XEst[ind]);
    ind = ind + NumUnc + 1;
  }

  // Fill correlations
  for(Int_t k = 0; k<NumUnc; k++){
    if(k == 0){myBlue->FillCor(k,&Cor00[0]);
    }else if(k ==  1){myBlue->FillCor(k,&Cor01[0]);
    }else if(k ==  2){myBlue->FillCor(k,&Cor02[0]);
    }else if(k ==  3){myBlue->FillCor(k,&Cor03[0]);
    }
  }


    printf("... LUMI_COMBINATION_1: Fix the input \n");

    myBlue->FixInp();

    printf("... LUMI_COMBINATION_1: Print estimators \n");
    for(Int_t i = 0; i<NumEst; i++){
      myBlue->PrintEst(i);
    }
    sprintf(Buffer,"LUMI_COMBINATION_1 %i",Flag);
    myBlue->PrintCompatEst(Buffer);
    printf("... LUMI_COMBINATION_1: The correlations");
    printf(" of the estimates in %%\n");
    myBlue->GetRho(LocRho);
    LocRho->operator*=(100);
    std::cout<<"LocRho[0][0] = "<<(*LocRho)[0][0]<<std::endl;
    myBlue->PrintMatrix(LocRho,"%+4.0f");

    myBlue->Solve();
    printf("... LUMI_COMBINATION_1: The Luminosity 17/18 Combination");
    printf(" Flag = %2i.\n",Flag);
    myBlue->PrintResult();
    myBlue->LatexResult("LUMI_COMBINATION_1");

    printf("... LUMI_COMBINATION_1: The correlations");
    printf(" of the observables in %%\n");
    myBlue->GetRhoRes(LocRhoRes);
    LocRhoRes->operator*=(100);
    myBlue->PrintMatrix(LocRhoRes,"%+4.0f");

    printf("... LUMI_COMBINATION_1: A difference is observed in the weights");
    printf(" that is under discussion with the authors \n");
    printf("... LUMI_COMBINATION_1: The observed weights\n");
    myBlue->GetWeight(LocWeight);
    myBlue->PrintMatrix(LocWeight," %+5.3f");
    // myBlue->PrintCompatObs();


  // Delete Object
  delete myBlue; myBlue = NULL;
  LocRho->Delete(); LocRho = NULL;
  LocRhoRes->Delete(); LocRhoRes = NULL;
  LocWeight->Delete(); LocWeight = NULL;
  return;
}
