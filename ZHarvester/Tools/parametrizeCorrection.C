#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <TFile.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
using namespace std;

void parametrizeCorrection() {
	// read file pileupfunccsvfile.csv, and fill correction = (ZeffData - ZeffMC) into histograms as a function of pile-up 
	ifstream in;
	//in.open("pileupfunccsvfile.csv");
	in.open("/afs/cern.ch/work/x/xniu/public/CMSSW_9_2_8/src/ZCountHarvest/CloneJob/TestResult/TestApr13/50-10invpb/pileupfunccsvfile.csv");

	auto f = TFile::Open("ParametrizeCorrection.root","RECREATE");	

   	TH2F h_BB("h_BB","data Z eff - MC Z eff, BB", 1000, 0., 100., 1000, 0., 0.1);
   	TH2F h_BE("h_BE","data Z eff - MC Z eff, BE", 1000, 0., 100., 1000, 0., 0.1);
   	TH2F h_EE("h_EE","data Z eff - MC Z eff, EE", 1000, 0., 100., 1000, 0., 0.1);

	string line;
	while(getline(in,line)) {
		stringstream    linestream(line);
		string          value;
		vector <string> pileup_eff;
		bool badMarker = 0;

		while(getline(linestream,value,',')) {
			if(value.compare("0.0") == 0) {
				badMarker = 1;
				break;
			}
			pileup_eff.push_back(value);
		}

		if(badMarker == 1) continue;

	        float pileup = strtof((pileup_eff[0]).c_str(),0);

	        float ZdataEffBB = strtof((pileup_eff[1]).c_str(),0);
        	float ZdataEffBE = strtof((pileup_eff[2]).c_str(),0);
	        float ZdataEffEE = strtof((pileup_eff[3]).c_str(),0);

	        float ZMCEffBB = strtof((pileup_eff[4]).c_str(),0);
	        float ZMCEffBE = strtof((pileup_eff[5]).c_str(),0);
	        float ZMCEffEE = strtof((pileup_eff[6]).c_str(),0);

                //cout<<"pileup :"<<pileup<<endl;
                //cout<<"EffdataBB :"<<ZdataEffBB<<endl;
                //cout<<"EffdataBE :"<<ZdataEffBE<<endl;
                //cout<<"EffdataEE :"<<ZdataEffEE<<endl;
                //cout<<"EffMCBB :"<<ZMCEffBB<<endl;
                //cout<<"EffMCBE :"<<ZMCEffBE<<endl;
                //cout<<"EffMCEE :"<<ZMCEffEE<<endl;

                h_BB.Fill(pileup, ZdataEffBB - ZMCEffBB);
                h_BE.Fill(pileup, ZdataEffBE - ZMCEffBE);
                h_EE.Fill(pileup, ZdataEffEE - ZMCEffEE);

	}
	in.close();

	//Fit three histograms to parametrize the correction as a function of pile-up 
	TCanvas c("c");

        TF1 f_BB("f_BB","[0]+[1]*x",0,60);
        h_BB.Fit("f_BB");

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

	TF1 f_BE("f_BE","[0]+[1]*x",0,60);
	h_BE.Fit("f_BE");

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

        TF1 f_EE("f_EE","[0]+[1]*x",0,60);
        h_EE.Fit("f_EE");

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
