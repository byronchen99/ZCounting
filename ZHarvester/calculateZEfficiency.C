#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <stdexcept>
#include <iomanip>                  // functions to format standard I/O
//#include <fstream>                  // functions for file I/O
//#include <sstream>                  // class for parsing strings

//(66,116), 27, 2.4: 2.0908e+06/6.11118e+06     0.342
//(66,116), 27, 0.9: 516500/6.11118e+06         0.0845
//(66,116), 30, 2.4: 1.84394e+06/6.11118e+06    0.302
//(66,116), 30, 0.9: 484740/6.11118e+06         0.0793

//76X
//6.11118e+06 1.84394e+06
//484740, 713373, 645823
//const std::vector<float> genAcc = {0.301732, 0.079320, 0.116732, 0.105679};//pre-FSR generator level acceptance
//92X
//8.72492e+06 2.6231e+06
//679707, 1.02256e+06, 920835
const std::vector<float> genAcc = {0.300644, 0.077904, 0.117200, 0.105541};//pre-FSR generator level acceptance

std::vector<float> getMCMuEff(int nPU){
  // for 92X
  if(nPU <= 1) {std::vector<float> v = {0.903776,0.879171,0.988622,0.99635,0.980652,0.998535};return v;}
  else if(nPU == 2) {std::vector<float> v = {0.91572,0.883317,0.987575,0.996304,0.981013,0.99819};return v;}
  else if(nPU == 3) {std::vector<float> v = {0.922173,0.88394,0.987494,0.996097,0.982268,0.998117};return v;}
  else if(nPU == 4) {std::vector<float> v = {0.925566,0.882926,0.987725,0.995858,0.983661,0.998132};return v;}
  else if(nPU == 5) {std::vector<float> v = {0.927116,0.881235,0.987902,0.995654,0.984852,0.998097};return v;}
  else if(nPU == 6) {std::vector<float> v = {0.92748,0.879438,0.987957,0.995509,0.985708,0.998003};return v;}
  else if(nPU == 7) {std::vector<float> v = {0.927107,0.87784,0.987941,0.995421,0.986244,0.997886};return v;}
  else if(nPU == 8) {std::vector<float> v = {0.926336,0.876547,0.987912,0.995375,0.986546,0.997775};return v;}
  else if(nPU == 9) {std::vector<float> v = {0.925414,0.875537,0.987898,0.995355,0.986707,0.997681};return v;}
  else if(nPU == 10) {std::vector<float> v = {0.924492,0.874734,0.987904,0.99535,0.986799,0.997603};return v;}
  else if(nPU == 11) {std::vector<float> v = {0.923649,0.874054,0.987924,0.995352,0.986859,0.997533};return v;}
  else if(nPU == 12) {std::vector<float> v = {0.922909,0.873427,0.987949,0.995356,0.986908,0.997466};return v;}
  else if(nPU == 13) {std::vector<float> v = {0.922268,0.872807,0.987973,0.995358,0.986951,0.997397};return v;}
  else if(nPU == 14) {std::vector<float> v = {0.921709,0.872168,0.987993,0.995358,0.986991,0.997325};return v;}
  else if(nPU == 15) {std::vector<float> v = {0.921211,0.8715,0.988006,0.995352,0.98703,0.997251};return v;}
  else if(nPU == 16) {std::vector<float> v = {0.920752,0.870802,0.988012,0.995341,0.987068,0.997175};return v;}
  else if(nPU == 17) {std::vector<float> v = {0.920317,0.870081,0.988012,0.995324,0.987107,0.997099};return v;}
  else if(nPU == 18) {std::vector<float> v = {0.919891,0.869344,0.988007,0.995302,0.987147,0.997021};return v;}
  else if(nPU == 19) {std::vector<float> v = {0.919465,0.8686,0.987998,0.995274,0.987188,0.996943};return v;}
  else if(nPU == 20) {std::vector<float> v = {0.919034,0.867853,0.987989,0.995242,0.98723,0.996862};return v;}
  else if(nPU == 21) {std::vector<float> v = {0.918596,0.867106,0.987981,0.995208,0.987275,0.996779};return v;}
  else if(nPU == 22) {std::vector<float> v = {0.918153,0.86636,0.987977,0.995173,0.987321,0.996692};return v;}
  else if(nPU == 23) {std::vector<float> v = {0.917706,0.865614,0.987978,0.995136,0.987369,0.996601};return v;}
  else if(nPU == 24) {std::vector<float> v = {0.917258,0.864864,0.987986,0.9951,0.987418,0.996504};return v;}
  else if(nPU == 25) {std::vector<float> v = {0.916812,0.864108,0.988003,0.995065,0.987469,0.996402};return v;}
  else if(nPU == 26) {std::vector<float> v = {0.916368,0.863343,0.988028,0.99503,0.98752,0.996296};return v;}
  else if(nPU == 27) {std::vector<float> v = {0.915927,0.862565,0.98806,0.994997,0.987571,0.996184};return v;}
  else if(nPU == 28) {std::vector<float> v = {0.915486,0.861772,0.9881,0.994966,0.987622,0.99607};return v;}
  else if(nPU == 29) {std::vector<float> v = {0.915043,0.860963,0.988147,0.994937,0.987674,0.995952};return v;}
  else if(nPU == 30) {std::vector<float> v = {0.914596,0.860137,0.988199,0.994909,0.987725,0.995833};return v;}
  else if(nPU == 31) {std::vector<float> v = {0.914143,0.859294,0.988256,0.994884,0.987777,0.995714};return v;}
  else if(nPU == 32) {std::vector<float> v = {0.91368,0.858434,0.988316,0.994861,0.987829,0.995595};return v;}
  else if(nPU == 33) {std::vector<float> v = {0.913207,0.857559,0.988379,0.99484,0.987883,0.995477};return v;}
  else if(nPU == 34) {std::vector<float> v = {0.912722,0.85667,0.988443,0.99482,0.987939,0.995361};return v;}
  else if(nPU == 35) {std::vector<float> v = {0.912226,0.855769,0.988507,0.994802,0.987997,0.995248};return v;}
  else if(nPU == 36) {std::vector<float> v = {0.911717,0.85486,0.988569,0.994785,0.98806,0.995137};return v;}
  else if(nPU == 37) {std::vector<float> v = {0.911198,0.853944,0.988628,0.994767,0.988127,0.99503};return v;}
  else if(nPU == 38) {std::vector<float> v = {0.91067,0.853026,0.988682,0.99475,0.9882,0.994927};return v;}
  else if(nPU == 39) {std::vector<float> v = {0.910135,0.852106,0.98873,0.994731,0.988279,0.994829};return v;}
  else if(nPU == 40) {std::vector<float> v = {0.909594,0.85119,0.988771,0.994711,0.988366,0.994735};return v;}
  else if(nPU == 41) {std::vector<float> v = {0.90905,0.850278,0.988805,0.994689,0.988462,0.994646};return v;}
  else if(nPU == 42) {std::vector<float> v = {0.908506,0.849374,0.988833,0.994665,0.988567,0.994563};return v;}
  else if(nPU == 43) {std::vector<float> v = {0.907962,0.848477,0.988855,0.994639,0.988681,0.994486};return v;}
  else if(nPU == 44) {std::vector<float> v = {0.907422,0.847587,0.988873,0.994612,0.988806,0.994416};return v;}
  else if(nPU == 45) {std::vector<float> v = {0.906886,0.846703,0.988889,0.994583,0.988939,0.994354};return v;}
  else if(nPU == 46) {std::vector<float> v = {0.906357,0.845819,0.988905,0.994553,0.989081,0.9943};return v;}
  else if(nPU == 47) {std::vector<float> v = {0.905835,0.844931,0.988924,0.994521,0.98923,0.994254};return v;}
  else if(nPU == 48) {std::vector<float> v = {0.90532,0.844031,0.988948,0.99449,0.989384,0.994218};return v;}
  else if(nPU == 49) {std::vector<float> v = {0.904813,0.843111,0.988978,0.994458,0.989542,0.994192};return v;}
  else if(nPU == 50) {std::vector<float> v = {0.904312,0.842162,0.989016,0.994428,0.989701,0.994176};return v;}
  else if(nPU == 51) {std::vector<float> v = {0.903818,0.841176,0.989063,0.9944,0.98986,0.994169};return v;}
  else if(nPU == 52) {std::vector<float> v = {0.903328,0.840146,0.989119,0.994374,0.990017,0.99417};return v;}
  else if(nPU == 53) {std::vector<float> v = {0.902841,0.839067,0.989181,0.994352,0.99017,0.994178};return v;}
  else if(nPU == 54) {std::vector<float> v = {0.902356,0.837939,0.989249,0.994334,0.990319,0.99419};return v;}
  else if(nPU == 55) {std::vector<float> v = {0.901873,0.836763,0.989318,0.994323,0.990464,0.994204};return v;}
  else if(nPU == 56) {std::vector<float> v = {0.901392,0.835548,0.989385,0.994319,0.990604,0.994216};return v;}
  else if(nPU == 57) {std::vector<float> v = {0.900914,0.834303,0.989446,0.994325,0.990739,0.994222};return v;}
  else if(nPU == 58) {std::vector<float> v = {0.900443,0.833044,0.989494,0.994341,0.990869,0.99422};return v;}
  else if(nPU == 59) {std::vector<float> v = {0.899984,0.831788,0.989526,0.99437,0.990996,0.994206};return v;}
  else if(nPU >= 60) {std::vector<float> v = {0.899547,0.830552,0.989536,0.994413,0.991119,0.994177};return v;}
  else std::invalid_argument("nPU");
}


std::vector<float> getMCZCount(int nPU){
  if(nPU <= 1)  {std::vector<float> v = {8.72492e+06,516096,823194,742601}; return v;}
  else if(nPU == 2)  {std::vector<float> v = {8.72492e+06,596739,911799,840891}; return v;}
  else if(nPU == 3)  {std::vector<float> v = {8.72492e+06,623743,935775,872603}; return v;}
  else if(nPU == 4)  {std::vector<float> v = {8.72492e+06,630984,943684,883225}; return v;}
  else if(nPU == 5)  {std::vector<float> v = {8.72492e+06,632190,948335,886750}; return v;}
  else if(nPU == 6)  {std::vector<float> v = {8.72492e+06,632165,950883,886036}; return v;}
  else if(nPU == 7)  {std::vector<float> v = {8.72492e+06,632178,953308,883848}; return v;}
  else if(nPU == 8)  {std::vector<float> v = {8.72492e+06,632436,955560,881127}; return v;}
  else if(nPU == 9)  {std::vector<float> v = {8.72492e+06,633894,958783,879818}; return v;}
  else if(nPU == 10) {std::vector<float> v = {8.72492e+06,633748,959266,876868}; return v;}
  else if(nPU == 11) {std::vector<float> v = {8.72492e+06,633678,960399,875444}; return v;}
  else if(nPU == 12) {std::vector<float> v = {8.72492e+06,634011,962047,875350}; return v;}
  else if(nPU == 13) {std::vector<float> v = {8.72492e+06,633505,962831,875055}; return v;}
  else if(nPU == 14) {std::vector<float> v = {8.72492e+06,632363,960282,872590}; return v;}
  else if(nPU == 15) {std::vector<float> v = {8.72492e+06,632952,962142,873900}; return v;}
  else if(nPU == 16) {std::vector<float> v = {8.72492e+06,632415,960248,872234}; return v;}
  else if(nPU == 17) {std::vector<float> v = {8.72492e+06,631510,958880,870788}; return v;}
  else if(nPU == 18) {std::vector<float> v = {8.72492e+06,631830,961379,872395}; return v;}
  else if(nPU == 19) {std::vector<float> v = {8.72492e+06,631472,959326,870427}; return v;}
  else if(nPU == 20) {std::vector<float> v = {8.72492e+06,631439,957891,869028}; return v;}
  else if(nPU == 21) {std::vector<float> v = {8.72492e+06,631784,960217,870427}; return v;}
  else if(nPU == 22) {std::vector<float> v = {8.72492e+06,631154,957535,867917}; return v;}
  else if(nPU == 23) {std::vector<float> v = {8.72492e+06,631341,958258,868172}; return v;}
  else if(nPU == 24) {std::vector<float> v = {8.72492e+06,631765,959393,868835}; return v;}
  else if(nPU == 25) {std::vector<float> v = {8.72492e+06,631754,956812,866734}; return v;}
  else if(nPU == 26) {std::vector<float> v = {8.72492e+06,631386,957147,866691}; return v;}
  else if(nPU == 27) {std::vector<float> v = {8.72492e+06,631763,954960,865095}; return v;}
  else if(nPU == 28) {std::vector<float> v = {8.72492e+06,631576,955618,865423}; return v;}
  else if(nPU == 29) {std::vector<float> v = {8.72492e+06,631508,955294,865093}; return v;}
  else if(nPU == 30) {std::vector<float> v = {8.72492e+06,631416,954842,864628}; return v;}
  else if(nPU == 31) {std::vector<float> v = {8.72492e+06,630523,952077,862289}; return v;}
  else if(nPU == 32) {std::vector<float> v = {8.72492e+06,631139,951926,862174}; return v;}
  else if(nPU == 33) {std::vector<float> v = {8.72492e+06,631082,950625,860996}; return v;}
  else if(nPU == 34) {std::vector<float> v = {8.72492e+06,631250,950797,860860}; return v;}
  else if(nPU == 35) {std::vector<float> v = {8.72492e+06,630948,949853,859699}; return v;}
  else if(nPU == 36) {std::vector<float> v = {8.72492e+06,630887,949820,859152}; return v;}
  else if(nPU == 37) {std::vector<float> v = {8.72492e+06,630400,948859,857697}; return v;}
  else if(nPU == 38) {std::vector<float> v = {8.72492e+06,630403,947874,856273}; return v;}
  else if(nPU == 39) {std::vector<float> v = {8.72492e+06,631131,946801,854931}; return v;}
  else if(nPU == 40) {std::vector<float> v = {8.72492e+06,631420,947889,854910}; return v;}
  else if(nPU == 41) {std::vector<float> v = {8.72492e+06,630610,945831,852258}; return v;}
  else if(nPU == 42) {std::vector<float> v = {8.72492e+06,631077,945210,851033}; return v;}
  else if(nPU == 43) {std::vector<float> v = {8.72492e+06,631527,944289,849561}; return v;}
  else if(nPU == 44) {std::vector<float> v = {8.72492e+06,631679,945173,849299}; return v;}
  else if(nPU == 45) {std::vector<float> v = {8.72492e+06,631560,944206,847719}; return v;}
  else if(nPU == 46) {std::vector<float> v = {8.72492e+06,631229,942332,845477}; return v;}
  else if(nPU == 47) {std::vector<float> v = {8.72492e+06,631644,942478,844979}; return v;}
  else if(nPU == 48) {std::vector<float> v = {8.72492e+06,630828,941190,843214}; return v;}
  else if(nPU == 49) {std::vector<float> v = {8.72492e+06,630447,940593,842165}; return v;}
  else if(nPU == 50) {std::vector<float> v = {8.72492e+06,630208,939802,841180} ; return v;}
  else if(nPU == 51) {std::vector<float> v = {8.72492e+06,629015,939066,840059} ; return v;}
  else if(nPU == 52) {std::vector<float> v = {8.72492e+06,627925,937917,838878} ; return v;}
  else if(nPU == 53) {std::vector<float> v = {8.72492e+06,627224,937067,838197} ; return v;}
  else if(nPU == 54) {std::vector<float> v = {8.72492e+06,625906,936442,837691} ; return v;}
  else if(nPU == 55) {std::vector<float> v = {8.72492e+06,624410,934948,836750} ; return v;}
  else if(nPU == 56) {std::vector<float> v = {8.72492e+06,623230,933600,836242} ; return v;}
  else if(nPU == 57) {std::vector<float> v = {8.72492e+06,621701,932435,836026} ; return v;}
  else if(nPU == 58) {std::vector<float> v = {8.72492e+06,619970,929849,835021} ; return v;}
  else if(nPU == 59) {std::vector<float> v = {8.72492e+06,618714,927656,834712} ; return v;}
  else if(nPU >= 60) {std::vector<float> v = {8.72492e+06,617394,925427,834627}; return v;}
  else std::invalid_argument("nPU");
}

float effFormula(float effHLT1, float effHLT2, float effSIT1, float effSIT2, float effSta1, float effSta2){
  return (1. - (1. - effHLT1)*(1. - effHLT2))* effSIT1*effSIT2 * effSta1*effSta2;
}

float calculateZEfficiency(
  int   dimuRegion,	// 0 for inclusive, 1 for barrel-barrel, 2 for barrel-endcap, 3 for endcap-endcap 
  float meanPileup,
  float dtHLTEffB, float dtHLTEffE, float dtSITEffB, float dtSITEffE, float dtStaEffB, float dtStaEffE
){

  auto vMCZCount = getMCZCount(int(meanPileup));
  auto vMCMuEff  = getMCMuEff(int(meanPileup));

  float nEvtsv = vMCZCount[0];
  float nSelCorrv = 0.;

  if(dimuRegion == 0){

    float dtZEffBB = effFormula( dtHLTEffB, dtHLTEffB, dtSITEffB, dtSITEffB, dtStaEffB, dtStaEffB );
    float dtZEffBE = effFormula( dtHLTEffB, dtHLTEffE, dtSITEffB, dtSITEffE, dtStaEffB, dtStaEffE );
    float dtZEffEE = effFormula( dtHLTEffE, dtHLTEffE, dtSITEffE, dtSITEffE, dtStaEffE, dtStaEffE );

    float mcZEffBB = effFormula( vMCMuEff.at(0), vMCMuEff.at(0), vMCMuEff.at(2), vMCMuEff.at(2), vMCMuEff.at(4), vMCMuEff.at(4) );
    float mcZEffBE = effFormula( vMCMuEff.at(0), vMCMuEff.at(1), vMCMuEff.at(2), vMCMuEff.at(3), vMCMuEff.at(4), vMCMuEff.at(5) );
    float mcZEffEE = effFormula( vMCMuEff.at(1), vMCMuEff.at(1), vMCMuEff.at(3), vMCMuEff.at(3), vMCMuEff.at(5), vMCMuEff.at(5) );

    float corrBB = dtZEffBB/mcZEffBB;
    float corrBE = dtZEffBE/mcZEffBE;
    float corrEE = dtZEffEE/mcZEffEE;

    nSelCorrv = vMCZCount[1] * corrBB + vMCZCount[2] * corrBE + vMCZCount[3] * corrEE;  

  }
  else if(dimuRegion == 1){

    float dtZEffBB = effFormula( dtHLTEffB, dtHLTEffB, dtSITEffB, dtSITEffB, dtStaEffB, dtStaEffB );
    float mcZEffBB = effFormula( vMCMuEff.at(0), vMCMuEff.at(0), vMCMuEff.at(2), vMCMuEff.at(2), vMCMuEff.at(4), vMCMuEff.at(4) );
    float corrBB = dtZEffBB/mcZEffBB;

    nSelCorrv = vMCZCount[1] * corrBB;

  }
  else if(dimuRegion == 2){

    float dtZEffBE = effFormula( dtHLTEffB, dtHLTEffE, dtSITEffB, dtSITEffE, dtStaEffB, dtStaEffE );
    float mcZEffBE = effFormula( vMCMuEff.at(0), vMCMuEff.at(1), vMCMuEff.at(2), vMCMuEff.at(3), vMCMuEff.at(4), vMCMuEff.at(5) );
    Float_t corrBE = dtZEffBE/mcZEffBE;

    nSelCorrv = vMCZCount[2] * corrBE;

  }
  else if(dimuRegion == 3){

    Float_t dtZEffEE = effFormula( dtHLTEffE, dtHLTEffE, dtSITEffE, dtSITEffE, dtStaEffE, dtStaEffE );
    Float_t mcZEffEE = effFormula( vMCMuEff.at(1), vMCMuEff.at(1), vMCMuEff.at(3), vMCMuEff.at(3), vMCMuEff.at(5), vMCMuEff.at(5) );
    Float_t corrEE = dtZEffEE/mcZEffEE;

    nSelCorrv = vMCZCount[3] * corrEE;
  }
  else std::invalid_argument("dimuRegion");


  // compute acceptances
  return (nSelCorrv/nEvtsv)/genAcc[dimuRegion];


  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================    
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << " Z -> mu mu" << endl;
  cout << endl;

  cout << "   ================================================" << endl;
  cout << endl;
  cout << "    *** Acceptance ***" << endl;
  cout << "   AllSFCorrected: " << setw(12) << nSelCorrv   << " / " << nEvtsv << " = " << nSelCorrv/nEvtsv << endl;
  cout << endl;
} 
