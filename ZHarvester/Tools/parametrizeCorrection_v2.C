#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <stdexcept>
#include <iomanip>                  // functions to format standard I/O
//#include <fstream>                  // functions for file I/O
//#include <sstream>                  // class for parsing strings
#include <TFile.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>

//(66,116), 27, 2.4: 2.0908e+06/6.11118e+06     0.342
//(66,116), 27, 0.9: 516500/6.11118e+06         0.0845
//(66,116), 30, 2.4: 1.84394e+06/6.11118e+06    0.302
//(66,116), 30, 0.9: 484740/6.11118e+06         0.0793

//76X: (66,116),30
//6.11118e+06 1.84394e+06
//484740, 713373, 645823
//const std::vector<float> genAcc = {0.301732, 0.079320, 0.116732, 0.105679};//pre-FSR generator level acceptance
//92X: (66,116),27
//8.72492e+06 2.96573e+06
//720702, 1.27156e+06, 973470
//const std::vector<float> genAcc = {0.339915, 0.082603, 0.145739, 0.111573};
//92X: (66,116),30
//8.72492e+06 2.6231e+06
//679707, 1.02256e+06, 920835
const std::vector<float> genAcc = {0.300644, 0.077904, 0.117200, 0.105541};//pre-FSR generator level acceptance

std::vector<float> getMCMuEffUnBiased(int nPU){
  if(nPU <= 1) {std::vector<float> v = {0.904623,0.879803,0.989084,0.997464,0.981463,0.998577};return v;}
  else if(nPU == 2) {std::vector<float> v = {0.916271,0.883756,0.98782,0.99578,0.981351,0.99833};return v;}
  else if(nPU == 3) {std::vector<float> v = {0.9225,0.884442,0.987531,0.994447,0.982572,0.998344};return v;}
  else if(nPU == 4) {std::vector<float> v = {0.925717,0.883565,0.987455,0.993697,0.984092,0.998422};return v;}
  else if(nPU == 5) {std::vector<float> v = {0.927181,0.881992,0.987363,0.993389,0.985411,0.998444};return v;}
  else if(nPU == 6) {std::vector<float> v = {0.92755,0.880271,0.987262,0.99333,0.986346,0.998404};return v;}
  else if(nPU == 7) {std::vector<float> v = {0.927241,0.878719,0.987198,0.993384,0.986913,0.998333};return v;}
  else if(nPU == 8) {std::vector<float> v = {0.926559,0.877454,0.987185,0.993477,0.98722,0.998256};return v;}
  else if(nPU == 9) {std::vector<float> v = {0.925726,0.876464,0.987212,0.993573,0.98738,0.998185};return v;}
  else if(nPU == 10) {std::vector<float> v = {0.924887,0.875681,0.987258,0.993659,0.987474,0.998122};return v;}
  else if(nPU == 11) {std::vector<float> v = {0.924118,0.875027,0.987305,0.993729,0.987546,0.998063};return v;}
  else if(nPU == 12) {std::vector<float> v = {0.923447,0.874439,0.987342,0.993785,0.987612,0.998007};return v;}
  else if(nPU == 13) {std::vector<float> v = {0.922874,0.873871,0.987365,0.993829,0.987678,0.99795};return v;}
  else if(nPU == 14) {std::vector<float> v = {0.922383,0.8733,0.987374,0.993863,0.98774,0.997893};return v;}
  else if(nPU == 15) {std::vector<float> v = {0.921955,0.872714,0.987372,0.993888,0.9878,0.997834};return v;}
  else if(nPU == 16) {std::vector<float> v = {0.921567,0.872112,0.987362,0.993903,0.987855,0.997775};return v;}
  else if(nPU == 17) {std::vector<float> v = {0.921203,0.871498,0.987348,0.993909,0.987906,0.997714};return v;}
  else if(nPU == 18) {std::vector<float> v = {0.92085,0.870877,0.987334,0.993905,0.987955,0.997654};return v;}
  else if(nPU == 19) {std::vector<float> v = {0.920499,0.870256,0.987321,0.993891,0.988002,0.997592};return v;}
  else if(nPU == 20) {std::vector<float> v = {0.920144,0.869638,0.987312,0.993869,0.988048,0.99753};return v;}
  else if(nPU == 21) {std::vector<float> v = {0.919786,0.869027,0.987307,0.993839,0.988093,0.997467};return v;}
  else if(nPU == 22) {std::vector<float> v = {0.919426,0.868422,0.98731,0.993804,0.988139,0.997402};return v;}
  else if(nPU == 23) {std::vector<float> v = {0.919068,0.867824,0.987319,0.993766,0.988185,0.997334};return v;}
  else if(nPU == 24) {std::vector<float> v = {0.918715,0.867231,0.987336,0.993726,0.988232,0.997263};return v;}
  else if(nPU == 25) {std::vector<float> v = {0.918371,0.866641,0.987359,0.993686,0.988279,0.997189};return v;}
  else if(nPU == 26) {std::vector<float> v = {0.918036,0.866051,0.987389,0.993647,0.988326,0.997111};return v;}
  else if(nPU == 27) {std::vector<float> v = {0.917713,0.86546,0.987424,0.993611,0.988374,0.997029};return v;}
  else if(nPU == 28) {std::vector<float> v = {0.917399,0.864865,0.987463,0.993579,0.988421,0.996944};return v;}
  else if(nPU == 29) {std::vector<float> v = {0.917093,0.864266,0.987506,0.993552,0.988467,0.996858};return v;}
  else if(nPU == 30) {std::vector<float> v = {0.916793,0.863662,0.987551,0.99353,0.988513,0.99677};return v;}
  else if(nPU == 31) {std::vector<float> v = {0.916496,0.863054,0.987599,0.993513,0.988559,0.996681};return v;}
  else if(nPU == 32) {std::vector<float> v = {0.916201,0.862442,0.987647,0.993501,0.988604,0.996593};return v;}
  else if(nPU == 33) {std::vector<float> v = {0.915906,0.861829,0.987696,0.993493,0.988648,0.996505};return v;}
  else if(nPU == 34) {std::vector<float> v = {0.915609,0.861217,0.987745,0.993488,0.988693,0.996419};return v;}
  else if(nPU == 35) {std::vector<float> v = {0.915312,0.860608,0.987793,0.993486,0.988738,0.996333};return v;}
  else if(nPU == 36) {std::vector<float> v = {0.915013,0.860005,0.98784,0.993485,0.988783,0.99625};return v;}
  else if(nPU == 37) {std::vector<float> v = {0.914713,0.859412,0.987883,0.993485,0.98883,0.996167};return v;}
  else if(nPU == 38) {std::vector<float> v = {0.914415,0.858831,0.987924,0.993484,0.98888,0.996087};return v;}
  else if(nPU == 39) {std::vector<float> v = {0.914119,0.858267,0.987961,0.993482,0.988934,0.996008};return v;}
  else if(nPU == 40) {std::vector<float> v = {0.913826,0.857722,0.987994,0.993478,0.988992,0.99593};return v;}
  else if(nPU == 41) {std::vector<float> v = {0.913539,0.8572,0.988024,0.993472,0.989057,0.995854};return v;}
  else if(nPU == 42) {std::vector<float> v = {0.91326,0.856704,0.988052,0.993464,0.989129,0.995779};return v;}
  else if(nPU == 43) {std::vector<float> v = {0.912989,0.856234,0.98808,0.993454,0.989208,0.995705};return v;}
  else if(nPU == 44) {std::vector<float> v = {0.912731,0.855792,0.988109,0.993442,0.989297,0.995632};return v;}
  else if(nPU == 45) {std::vector<float> v = {0.912486,0.855378,0.988141,0.993428,0.989393,0.995561};return v;}
  else if(nPU == 46) {std::vector<float> v = {0.912256,0.85499,0.98818,0.993413,0.989498,0.995492};return v;}
  else if(nPU == 47) {std::vector<float> v = {0.912043,0.854625,0.988226,0.993396,0.989609,0.995425};return v;}
  else if(nPU == 48) {std::vector<float> v = {0.911847,0.854278,0.988283,0.993379,0.989726,0.995361};return v;}
  else if(nPU == 49) {std::vector<float> v = {0.91167,0.853943,0.988352,0.993363,0.989848,0.995301};return v;}
  else if(nPU == 50) {std::vector<float> v = {0.911511,0.853612,0.988432,0.993348,0.989973,0.995245};return v;}
  else if(nPU == 51) {std::vector<float> v = {0.91137,0.853279,0.988526,0.993336,0.990101,0.995193};return v;}
  else if(nPU == 52) {std::vector<float> v = {0.911244,0.852935,0.98863,0.993328,0.990229,0.995147};return v;}
  else if(nPU == 53) {std::vector<float> v = {0.911132,0.852572,0.988745,0.993325,0.990358,0.995106};return v;}
  else if(nPU == 54) {std::vector<float> v = {0.911033,0.852186,0.988867,0.993328,0.990487,0.995069};return v;}
  else if(nPU == 55) {std::vector<float> v = {0.910944,0.851771,0.988994,0.993338,0.990617,0.995037};return v;}
  else if(nPU == 56) {std::vector<float> v = {0.910865,0.851326,0.989121,0.993357,0.990748,0.995007};return v;}
  else if(nPU == 57) {std::vector<float> v = {0.910799,0.850849,0.989245,0.993385,0.990882,0.994979};return v;}
  else if(nPU == 58) {std::vector<float> v = {0.910748,0.850343,0.989363,0.993423,0.991019,0.994951};return v;}
  else if(nPU == 59) {std::vector<float> v = {0.910719,0.84981,0.98947,0.993473,0.991162,0.994921};return v;}
  else if(nPU >= 60) {std::vector<float> v = {0.910723,0.849254,0.989563,0.993535,0.991313,0.994887};return v;}
  else std::invalid_argument("nPU");
}

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

void parametrizeCorrection_v2(){

	auto f = TFile::Open("ParametrizeCorrection_v2.root","RECREATE");

        TH2F h_BB("h_BB","data Z eff - MC Z eff, BB", 1000, 0., 100., 1000, 0., 0.1);
        TH2F h_BE("h_BE","data Z eff - MC Z eff, BE", 1000, 0., 100., 1000, 0., 0.1);
        TH2F h_EE("h_EE","data Z eff - MC Z eff, EE", 1000, 0., 100., 1000, 0., 0.1);

	for(int iPU = 0; iPU < 60; iPU++) {
		auto vMCZCount = getMCZCount(iPU);
		auto vMCMuEff  = getMCMuEff(iPU);

		//without correlation
		float mcZEffBB0 = effFormula( vMCMuEff.at(0), vMCMuEff.at(0), vMCMuEff.at(2), vMCMuEff.at(2), vMCMuEff.at(4), vMCMuEff.at(4) );
		float mcZEffBE0 = effFormula( vMCMuEff.at(0), vMCMuEff.at(1), vMCMuEff.at(2), vMCMuEff.at(3), vMCMuEff.at(4), vMCMuEff.at(5) );
		float mcZEffEE0 = effFormula( vMCMuEff.at(1), vMCMuEff.at(1), vMCMuEff.at(3), vMCMuEff.at(3), vMCMuEff.at(5), vMCMuEff.at(5) );

		//with correlation
		float mcZEffBB  = (vMCZCount[1]/vMCZCount[0])/0.077904;
		float mcZEffBE  = (vMCZCount[2]/vMCZCount[0])/0.117200;
		float mcZEffEE  = (vMCZCount[3]/vMCZCount[0])/0.105541;

		//cout<<"mcZEffBB0: "<<mcZEffBB0<<endl;
		//cout<<"mcZEffBE0: "<<mcZEffBE0<<endl;
		//cout<<"mcZEffEE0: "<<mcZEffEE0<<endl;

                //cout<<"mcZEffBB: "<<mcZEffBB<<endl;
                //cout<<"mcZEffBE: "<<mcZEffBE<<endl;
                //cout<<"mcZEffEE: "<<mcZEffEE<<endl;

		//cout<<mcZEffBB0 - mcZEffBB<<endl;
		//cout<<mcZEffBE0 - mcZEffBE<<endl;
		//cout<<mcZEffEE0 - mcZEffEE<<endl;

                h_BB.Fill(iPU, mcZEffBB0 - mcZEffBB);
                h_BE.Fill(iPU, mcZEffBE0 - mcZEffBE);
                h_EE.Fill(iPU, mcZEffEE0 - mcZEffEE);
	}

        //Fit three histograms to parametrize the correction as a function of pile-up
        TCanvas c("c");

        //TF1 f_BB("f_BB","[0]+[1]*x",15,60);
        TF1 f_BB_1("f_BB_1","[0]+[1]*x",15,50);
        TF1 f_BB_2("f_BB_2","[2]+[3]*x",50,60);
        //h_BB.Fit("f_BB", "R");
        h_BB.Fit("f_BB_1", "R");
        h_BB.Fit("f_BB_2", "R+");

        h_BB.GetXaxis()->SetTitle("mean pile-up");
        h_BB.GetYaxis()->SetTitle("Z efficiency correction");
        h_BB.GetXaxis()->SetTitleSize(0.05);
        h_BB.GetYaxis()->SetTitleSize(0.05);
        //h_BB.GetXaxis()->SetLimits(xmin,xmax);
        //h_BB.GetYaxis()->SetRangeUser(ymin,ymax);
        h_BB.SetMarkerStyle(20);
        h_BB.SetMarkerColor(kBlue);
        h_BB.SetMarkerSize(1);
        h_BB.Draw("P");

        c.SaveAs("Corr_vs_PU_BB.png");
        c.Clear();

        //TF1 f_BE("f_BE","[0]+[1]*x",15,60);
        TF1 f_BE_1("f_BE_1","[0]+[1]*x",15,55);
        TF1 f_BE_2("f_BE_2","[2]+[3]*x",55,60);
        //h_BE.Fit("f_BE", "R");
        h_BE.Fit("f_BE_1", "R");
        h_BE.Fit("f_BE_2", "R+");

        h_BE.GetXaxis()->SetTitle("mean pile-up");
        h_BE.GetYaxis()->SetTitle("Z efficiency correction");
        h_BE.GetXaxis()->SetTitleSize(0.05);
        h_BE.GetYaxis()->SetTitleSize(0.05);
        //h_BE.GetXaxis()->SetLimits(xmin,xmax);
        //h_BE.GetYaxis()->SetRangeUser(ymin,ymax);
        h_BE.SetMarkerStyle(20);
        h_BE.SetMarkerColor(kBlue);
        h_BE.SetMarkerSize(1);
        h_BE.Draw("P");

        c.SaveAs("Corr_vs_PU_BE.png");
        c.Clear();

        //TF1 f_EE("f_EE","[0]+[1]*x",15,60);
        TF1 f_EE_1("f_EE_1","[0]+[1]*x",15,40);
        TF1 f_EE_2("f_EE_2","[2]+[3]*x",40,60);
        //h_EE.Fit("f_EE", "R");
        h_EE.Fit("f_EE_1", "R");
        h_EE.Fit("f_EE_2", "R+");

        h_EE.GetXaxis()->SetTitle("mean pile-up");
        h_EE.GetYaxis()->SetTitle("Z efficiency correction");
        h_EE.GetXaxis()->SetTitleSize(0.05);
        h_EE.GetYaxis()->SetTitleSize(0.05);
        //h_EE.GetXaxis()->SetLimits(xmin,xmax);
        //h_EE.GetYaxis()->SetRangeUser(ymin,ymax);
        h_EE.SetMarkerStyle(20);
        h_EE.SetMarkerColor(kBlue);
        h_EE.SetMarkerSize(1);
        h_EE.Draw("P");

        c.SaveAs("Corr_vs_PU_EE.png");
        c.Clear();

        f->Write();
}
