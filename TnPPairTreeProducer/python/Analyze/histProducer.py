from __future__ import division, print_function

# production of histograms from TnP trees
#   the histograms have the same format as the onse produced by DQM
#   and can therefore used with the same tools, e.g. ZCounting.py

import numpy as np
import pandas as pd
import os
import root_numpy as rn
import pdb
import ROOT
import argparse
import glob

#local imports
os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import tree_to_df


parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '--pt', default=30., type=float,
    help='specify pt cut'
)
parser.add_argument(
    '--tkIso', default=None, type=float,
    help='specify tracker isolation cut (0.05/0.1 for tight/loose Muon Iso)'
)
parser.add_argument(
    '--pfIso', default=None, type=float,
    help='specify particle flow based isolation cut (0.12/0.2 for tight/loose Muon Iso)'
)
parser.add_argument(
    '--dxy', default=None, type=float,
    help='specify dxy cut (0.2 for tight Muon ID)'
)
parser.add_argument(
    '--dz', default=None, type=float,
    help='specify dz cut (0.5 for tight Muon ID)'
)
parser.add_argument(
    '-o', '--output', default='./',
    help='specify output dir'
)

args = parser.parse_args()

output = args.output


if not os.path.isdir(output):
    os.mkdir(output)

MassMin_ = 30
MassMax_ = 301
MassBin_ = 4*int(MassMax_ - MassMin_)

LumiBin_ = 5000
LumiMin_ = 0.5
LumiMax_ = 5000.5

nPVMin_ = -0.5
nPVMax_ = 74.5
nPVBin_ = int(nPVMax_ - nPVMin_)

ptCut = args.pt if args.pt else 0.
etaCut = 2.4

dxyDistCut = None# 0.2
dzDistCut = None#0.5


# acceptance selection
selection = "pt1 >= {0} & pt2 >= {0} & eta1 < {1} & eta2 < {1} & dilepMass >= {2} & dilepMass <= {3} & delR > 0.8 & q1 != q2".format(
    ptCut, etaCut, MassMin_, MassMax_
    )

if dxyDistCut is not None:
    selection += ' & abs(dxy1 - dxy2) < {0}'.format(dxyDistCut)
if dzDistCut is not None:
    selection += ' & abs(dz1 - dz2) < {0}'.format(dzDistCut)

if args.tkIso is not None:
    selection += ' & (tkIso1 < {0} | (is2HLT & tkIso2 < {0}))'.format(tkIsoCut)
if args.pfIso is not None:
    selection += ' & (pfIso1 < {0} | (is2HLT & pfIso2 < {0}))'.format(pfIsoCut)
if args.dxy is not None:
    selection += ' & (dxy1 < {0} | (is2HLT & dxy2 < {0}))'.format(dxyCut)
if args.dz is not None:
    selection += ' & (dz1 < {0} | (is2HLT & dz2 < {0}))'.format(dzCut)

tkIsoCut = args.tkIso if args.tkIso else 99999
pfIsoCut = args.pfIso if args.pfIso else 99999
dxyCut = args.dxy if args.dxy else 99999
dzCut = args.dz if args.dz else 99999

# specify which branches to load
branches = ['dilepMass', 'nPV', 'eventNumber',
            'run', 'ls',
            'eta1', 'eta2',
            'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            # 'tkIso1', 'tkIso2', 'pfIso1', 'pfIso2',
            # 'dxy1', 'dxy2', 'dz1', 'dz2',
            # 'q1','q2',
            'nTrackerLayers2',
            'nValidPixelHits2'
            ]

storage="/pnfs/desy.de/cms/tier2/store/user/dwalter/SingleMuon/TnPPairTrees"

inputs = {

    '16B': glob.glob(storage+"_V13_UL2016B/210712_141723/0000/output_0_*"),
    '16C': glob.glob(storage+"_V13_UL2016C/210711_084521/0000/output_0_*"),
    '16D': glob.glob(storage+"_V13_UL2016D/210711_084540/0000/output_0_*"),
    '16E': glob.glob(storage+"_V13_UL2016E/210711_084604/0000/output_0_*"),
    '16F': glob.glob(storage+"_V13_ULpreVFP2016F/210710_081713/0000/output_0_*")
       +glob.glob(storage+"_V13_ULpostVFP2016F/210710_081732/0000/output_0_*"),
    '16G': glob.glob(storage+"_V13_UL2016G/210711_084625/0000/output_0_*"),
    '16H': glob.glob(storage+"_V13_UL2016H/210730_134940/0000/output_0_*"),
    
    '17B': glob.glob(storage+"_V13_UL2017B/210712_141903/0000/output_0_*"),
    '17C': glob.glob(storage+"_V13_UL2017C/210711_084818/0000/output_0_*"),
    '17D': glob.glob(storage+"_V13_UL2017D/210711_084856/0000/output_0_*"),
    '17E': glob.glob(storage+"_V13_UL2017E/210711_084911/0000/output_0_*"),
    '17F': glob.glob(storage+"_V13_UL2017F/210711_084923/0000/output_0_*"),
    '17H': glob.glob(storage+"_V14_UL2017H/211126_134142/0000/output_0_*"),

    '18A': glob.glob(storage+"_V13_UL2018A/210822_193705/0000/output_0_*"),
    '18B': glob.glob(storage+"_V13_UL2018B/211014_083356/0000/output_0_*"),
    '18C': glob.glob(storage+"_V13_UL2018C/210711_085057/0000/output_0_*"),
    '18D': glob.glob(storage+"_V13_UL2018D/210711_085109/000?/output_0_*"),
    # # Missing runs from 2018

}



for era, input in inputs.iteritems():
    print("##### era {0}".format(era))

    treeName = rn.list_trees(input[0])[0]
    print(">>> Load Events from {0} files".format(len(input)))
    _df = []
    for i in input:
        print("> file "+i)
        tfile = ROOT.TFile.Open(i)
        if(tfile.Get(treeName).GetEntries(selection) == 0):
            print("> no events in this file found! continue with next file")
            continue
        _df.append(tree_to_df(rn.root2array(i, treeName, selection=selection, branches=branches)))
    print(">>> Concatenate")
    df = pd.concat(_df)

    # df['q'] = (df['q1'] + df['q2']) / 2.

    df['isTrk'] = df['isTrk'] * (df['nTrackerLayers2'] >= 6) * (df['nValidPixelHits2'] >= 1)


    # for q, data in df.groupby('q'):
    #
    #     if q == -1:
    #         continue
    #         suffix="_Minus"
    #     elif q == 1:
    #         continue
    #         suffix="_Plus"
    #     else:
    #         suffix=""
    suffix = ""
    for run, data_run in df.groupby('run'):
        print(">>> run {0}".format(run))
        if max(data_run['ls']) > LumiMax_:
            print("WARNING: this run has {0} lumisections, where the maximum bin is {1}".format(max(data_run['ls']), LumiMax_))

        subdir = output+"/000"+str(run)[:-2]+"xx"
        if not os.path.isdir(subdir):
            os.mkdir(subdir)

        tfile = ROOT.TFile.Open(subdir+"/DQM_V0001_R000{0}__SingleMuon__Run20{1}.root".format(run, era),"UPDATE")
        
        # histogram with number of PV
        h_nPV = ROOT.TH2D("h_nPV", "number of PV",
            LumiBin_, LumiMin_, LumiMax_, nPVBin_, nPVMin_, nPVMax_);

        # Sel to Glo
        h_mass_SEL_pass_central = ROOT.TH2D("h_mass_SIT_pass_central"+suffix, "Muon SIT passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_pass_forward = ROOT.TH2D("h_mass_SIT_pass_forward"+suffix, "Muon SIT passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_fail_central = ROOT.TH2D("h_mass_SIT_fail_central"+suffix, "Muon SIT failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_SEL_fail_forward = ROOT.TH2D("h_mass_SIT_fail_forward"+suffix, "Muon SIT failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        # Trk to Sta
        h_mass_TRK_pass_central = ROOT.TH2D("h_mass_Trk_pass_central"+suffix, "Muon track passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_TRK_pass_forward = ROOT.TH2D("h_mass_Trk_pass_forward"+suffix, "Muon track passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_TRK_fail_central = ROOT.TH2D("h_mass_Trk_fail_central"+suffix, "Muon track failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_TRK_fail_forward = ROOT.TH2D("h_mass_Trk_fail_forward"+suffix, "Muon track failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        # Sta to Trk
        h_mass_STA_pass_central = ROOT.TH2D("h_mass_Sta_pass_central"+suffix, "Muon standalone passing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_STA_pass_forward = ROOT.TH2D("h_mass_Sta_pass_forward"+suffix, "Muon standalone passing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_STA_fail_central = ROOT.TH2D("h_mass_Sta_fail_central"+suffix, "Muon standalone failing probes central",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_STA_fail_forward = ROOT.TH2D("h_mass_Sta_fail_forward"+suffix, "Muon standalone failing probes forward",
            LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        h_mass_2HLTBB_Z = ROOT.TH2D("h_mass_2HLTBB_Z"+suffix, "Events where 2 muons pass HLT, both muons in barrel",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_2HLTEE_Z = ROOT.TH2D("h_mass_2HLTEE_Z"+suffix, "Events where 2 muons pass HLT, both muons in endcap",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_2HLTBE_Z = ROOT.TH2D("h_mass_2HLTBE_Z"+suffix, "Events where 2 muons pass HLT, one muon in barrel and one in endcap",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

        h_mass_1HLTBB_Z = ROOT.TH2D("h_mass_1HLTBB_Z"+suffix, "Events where 1 muon pass HLT, both muons in barrel",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_1HLTEE_Z = ROOT.TH2D("h_mass_1HLTEE_Z"+suffix, "Events where 1 muon pass HLT, both muons in endcap",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        h_mass_1HLTBE_Z = ROOT.TH2D("h_mass_1HLTBE_Z"+suffix, "Events where 1 muon pass HLT, one barrel one endcap",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
        
        data_run = data_run.sort_values(['ls', 'eventNumber'])

        # fill PV histogram
        last_evt = 0
        for ls, pv, evt in data_run[['ls', 'nPV', 'eventNumber']].values:
            # fill nPV one per event
            if evt == last_evt: 
                continue
            h_nPV.Fill(ls, pv)
            last_evt = evt

        # fill tag and probe histograms
        #--- barrel
        # fill tags if both muons passed HLT (2HLT category)
        for m,ls,eta2 in data_run.query("is2HLT==1 & abs(eta1) < 0.9")[['dilepMass','ls','eta2']].values:
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_TRK_pass_central.Fill(ls,m)
            h_mass_STA_pass_central.Fill(ls,m)
            if(abs(eta2) < 0.9):
                h_mass_2HLTBB_Z.Fill(ls,m)
            else:
                h_mass_2HLTBE_Z.Fill(ls,m)

        # fill probes
        for m,ls in data_run.query("is2HLT==1 & abs(eta2) < 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_TRK_pass_central.Fill(ls,m)
            h_mass_STA_pass_central.Fill(ls,m)


        for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) < 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_SEL_pass_central.Fill(ls,m)
            h_mass_TRK_pass_central.Fill(ls,m)
            h_mass_STA_pass_central.Fill(ls,m)
            if(abs(eta1) < 0.9):
                h_mass_1HLTBB_Z.Fill(ls,m)
            else:
                h_mass_1HLTBE_Z.Fill(ls,m)

        for m,ls in data_run.query("isGlo==1 & isTrk==1 & abs(eta2) < 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_fail_central.Fill(ls,m)
            h_mass_STA_pass_central.Fill(ls,m)
            h_mass_TRK_pass_central.Fill(ls,m)

        for m,ls in data_run.query("isTrk==1 & ((isGlo==0 & isSel==0 & is2HLT==0) | isSta==0) & abs(eta2)  < 0.9")[['dilepMass','ls']].values:
            h_mass_STA_fail_central.Fill(ls,m)

        for m,ls in data_run.query("isSta==1 & ((isGlo==0 & isSel==0 & is2HLT==0) | isTrk==0)  & abs(eta2)  < 0.9")[['dilepMass','ls']].values:
            h_mass_TRK_fail_central.Fill(ls,m)

        #--- endcap
        # fill tags if both muons passed HLT (2HLT category)
        for m,ls,eta2 in data_run.query("is2HLT==1 & abs(eta1) >= 0.9")[['dilepMass','ls','eta2']].values:
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_TRK_pass_forward.Fill(ls,m)
            h_mass_STA_pass_forward.Fill(ls,m)
            if(abs(eta2) >= 0.9):
                h_mass_2HLTEE_Z.Fill(ls,m)
            else:
                h_mass_2HLTBE_Z.Fill(ls,m)

        # fill probes
        for m,ls in data_run.query("is2HLT==1 & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_TRK_pass_forward.Fill(ls,m)
            h_mass_STA_pass_forward.Fill(ls,m)

        for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) >= 0.9")[['dilepMass','ls','eta1']].values:
            h_mass_SEL_pass_forward.Fill(ls,m)
            h_mass_TRK_pass_forward.Fill(ls,m)
            h_mass_STA_pass_forward.Fill(ls,m)
            if(abs(eta1) >= 0.9):
                h_mass_1HLTEE_Z.Fill(ls,m)
            else:
                h_mass_1HLTBE_Z.Fill(ls,m)

        for m,ls in data_run.query("isGlo==1 & isTrk==1 & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
            h_mass_SEL_fail_forward.Fill(ls,m)
            h_mass_TRK_pass_forward.Fill(ls,m)
            h_mass_STA_pass_forward.Fill(ls,m)

        for m,ls in data_run.query("isTrk==1 & ((isGlo==0 & isSel==0 & is2HLT==0) | isSta==0) & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
            h_mass_STA_fail_forward.Fill(ls,m)

        for m,ls in data_run.query("isSta==1 & ((isGlo==0 & isSel==0 & is2HLT==0) | isTrk==0) & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
            h_mass_TRK_fail_forward.Fill(ls,m)

        tfile.Write()
        tfile.Close()
        tfile.Delete()
