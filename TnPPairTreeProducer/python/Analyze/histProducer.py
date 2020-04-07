from __future__ import division, print_function

# production of histograms from TnP trees
#   the histograms have the same format as the onse produced by DQM
#   and can therefore used with the same tools, e.g. ZCounting.py

import numpy as np
import pandas as pd
import os
import root_numpy as rn
import pdb
from Utils.Utils import tree_to_df
import ROOT
import argparse
import glob

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)

args = parser.parse_args()

output = args.output


if os.path.isdir(output):
    print("output dir already exists, please remove or specify another name")
    exit()
os.mkdir(output)


MassBin_ = 50
MassMin_ = 66.
MassMax_ = 116.

LumiBin_ = 5000
LumiMin_ = 0.5
LumiMax_ = 5000.5

ptCut = 30.
tkIsoCut = 0.05

# acceptance selection
selection = 'pt1 > {0} & pt2 > {0} & dilepMass > {1} &  dilepMass < {2} & q1 != q2'.format(ptCut, MassMin_, MassMax_)

if tkIsoCut is not None:
    selection += ' & (tkIso1 < {0} | (is2HLT & tkIso2 < {0}))'.format(tkIsoCut)

# specify which branches to load
branches = ['dilepMass', 'nPV', 'run', 'ls',
            'eta1', 'eta2',
            'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'tkIso1', 'tkIso2',
            ]

storage="/pnfs/desy.de/cms/tier2/store/user/dwalter/SingleMuon/TnPPairTrees_V04/"

inputs = {
    'B':glob.glob(storage+"200324_212241/0000/output_0_*"),
    'C': glob.glob(storage+"200324_212252/0000/output_0_*"),
    'D': glob.glob(storage+"200406_123834/0000/output_0_*"),
    'E': glob.glob(storage+"200324_212318/0000/output_0_*"),
    'F': glob.glob(storage+"200324_212331/0000/output_0_*"),
    'H':glob.glob(storage+"200324_212345/0000/output_0_*")
}



for era, input in inputs.iteritems():
    print("##### era {0}".format(era))

    treeName = rn.list_trees(input[0])[0]
    print(">>> Load Events")
    _df = [tree_to_df(rn.root2array(i, treeName, selection=selection, branches=branches), 5) for i in input]
    print(">>> Concatenate")
    df = pd.concat(_df)

    ### additional requirement at Sel:
    # case for is2HLT and only tag fails additional cut: tag and probe has to be chanched
    selected = df['is2HLT'] & (df['tkIso1'] > tkIsoCut) & (df['tkIso2'] < tkIsoCut)
    to_switch = df[selected]
    df = df[selected==False]
    to_switch = to_switch.rename(index=str, columns={
        'tkIso1': 'tkIso2', 'tkIso2': 'tkIso1',
        'pt1': 'pt2', 'pt2': 'pt1',
        'eta1': 'eta2', 'eta2': 'eta1'
        })
    df = pd.concat([df,to_switch],sort=True)
    # downgrade is2HLT and isSel to isGlo
    df['isGlo'] = df['isGlo'] + df['is2HLT'] * (df['tkIso2'] > tkIsoCut) + df['isSel'] * (df['tkIso2'] > tkIsoCut)
    df['is2HLT'] = df['is2HLT'] * (df['tkIso2'] < tkIsoCut)
    df['isSel'] = df['isSel'] * (df['tkIso2'] < tkIsoCut)



    for run, data_run in df.groupby('run'):
        print(">>> run {0}".format(run))
        if max(data_run['ls']) > LumiMax_:
            print("WARNING: this run has {0} lumisections, where the maximum bin is {1}".format(max(data_run['ls']), LumiMax_))

        subdir = output+"/000"+str(run)[:-2]+"xx"
        if not os.path.isdir(subdir):
            os.mkdir(subdir)

        tfile = ROOT.TFile.Open(subdir+"/DQM_V0001_R000{0}__SingleMuon__Run2017{1}-09Aug2019_UL2017-v1__DQMIO.root".format(run, era),"RECREATE")

        h_mass_HLT_pass_central = ROOT.TH2D("h_mass_HLT_pass_central", "Muon HLT passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_HLT_pass_forward = ROOT.TH2D("h_mass_HLT_pass_forward", "Muon HLT passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_HLT_fail_central = ROOT.TH2D("h_mass_HLT_fail_central", "Muon HLT failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_HLT_fail_forward = ROOT.TH2D("h_mass_HLT_fail_forward", "Muon HLT failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        h_mass_SEL_pass_central = ROOT.TH2D("h_mass_SIT_pass_central", "Muon SIT passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_pass_forward = ROOT.TH2D("h_mass_SIT_pass_forward", "Muon SIT passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_fail_central = ROOT.TH2D("h_mass_SIT_fail_central", "Muon SIT failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_fail_forward = ROOT.TH2D("h_mass_SIT_fail_forward", "Muon SIT failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        h_mass_GLO_pass_central = ROOT.TH2D("h_mass_Glo_pass_central", "Muon Glo passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_GLO_pass_forward = ROOT.TH2D("h_mass_Glo_pass_forward", "Muon Glo passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_GLO_fail_central = ROOT.TH2D("h_mass_Glo_fail_central", "Muon Glo failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_GLO_fail_forward = ROOT.TH2D("h_mass_Glo_fail_forward", "Muon Glo failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        h_mass_yield_Z = ROOT.TH2D("h_mass_yield_Z", "reconstructed Z bosons", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_yieldBB_Z = ROOT.TH1D("h_yieldBB_Z", "reconstructed Z bosons in barrel", LumiBin_, LumiMin_, LumiMax_);
        h_yieldEE_Z = ROOT.TH1D("h_yieldEE_Z", "reconstructed Z bosons in endcap", LumiBin_, LumiMin_, LumiMax_);

        #fill tag and probe histograms
        #--- barrel
        # fill tags if both muons passed HLT (2HLT category)
        for m,ls in data_run.query("is2HLT==1 & abs(eta1) < 0.9")[['dilepMass','ls']].values:
            h_mass_HLT_pass_central.Fill(ls,m)
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_GLO_pass_central.Fill(ls,m)
        # fill probes
        for m,ls,eta1 in data_run.query("is2HLT==1 & abs(eta2) < 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_HLT_pass_central.Fill(ls,m)
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_GLO_pass_central.Fill(ls,m)
            h_mass_yield_Z.Fill(ls,m)
            if(abs(eta1) < 0.9):
                h_yieldBB_Z.Fill(ls)

        for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) < 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_HLT_fail_central.Fill(ls,m)
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_GLO_pass_central.Fill(ls,m)
            h_mass_yield_Z.Fill(ls,m)
            if(abs(eta1) < 0.9):
                h_yieldBB_Z.Fill(ls)

        for m,ls in data_run.query("isGlo==1 & abs(eta2) < 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_fail_central.Fill(ls,m)
            h_mass_GLO_pass_central.Fill(ls,m)

        for m,ls in data_run.query("(isSta==1 | isTrk==1) & abs(eta2)  < 0.9")[['dilepMass','ls']].values:
            h_mass_GLO_fail_central.Fill(ls,m)

        #--- endcap
        # fill tags if both muons passed HLT (2HLT category)
        for m,ls in data_run.query("is2HLT==1 & abs(eta1) > 0.9")[['dilepMass','ls']].values:
            h_mass_HLT_pass_forward.Fill(ls,m)
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_GLO_pass_forward.Fill(ls,m)

        # fill probes
        for m,ls,eta1 in data_run.query("is2HLT==1 & abs(eta2) > 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_HLT_pass_forward.Fill(ls,m)
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_GLO_pass_forward.Fill(ls,m)
            h_mass_yield_Z.Fill(ls,m)
            if(abs(eta1) > 0.9):
                h_yieldEE_Z.Fill(ls)

        for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) > 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_HLT_fail_forward.Fill(ls,m)
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_GLO_pass_forward.Fill(ls,m)
            h_mass_yield_Z.Fill(ls,m)
            if(abs(eta1) > 0.9):
                h_yieldEE_Z.Fill(ls)

        for m,ls in data_run.query("isGlo==1 & abs(eta2) > 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_fail_forward.Fill(ls,m)
            h_mass_GLO_pass_forward.Fill(ls,m)

        for m,ls in data_run.query("(isSta==1 | isTrk==1) & abs(eta2)  > 0.9")[['dilepMass','ls']].values:
            h_mass_GLO_fail_forward.Fill(ls,m)

        tfile.Write()
        tfile.Close()
        tfile.Delete()