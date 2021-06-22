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
MassMax_ = 300
MassBin_ = int(MassMax_ - MassMin_)

LumiBin_ = 5000
LumiMin_ = 0.5
LumiMax_ = 5000.5

ptCut = args.pt if args.pt else 0.

dxyDistCut = None# 0.2
dzDistCut = None#0.5


# acceptance selection
selection = 'pt1 >= {0} & pt2 >= {0} & dilepMass >= {1} &  dilepMass <= {2}'.format(ptCut, MassMin_, MassMax_)

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
branches = ['dilepMass', 'nPV',
            'run', 'ls',
            'eta1', 'eta2',
            'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            'tkIso1', 'tkIso2', 'pfIso1', 'pfIso2',
            'dxy1', 'dxy2', 'dz1', 'dz2',
            'q1','q2'
            ]

storage="/pnfs/desy.de/cms/tier2/store/user/dwalter/SingleMuon/TnPPairTrees"

inputs = {
    # '16B': glob.glob(storage+"_V11_UL2016B/210417_163807/0000/output_0_*"),
    # '16C': glob.glob(storage+"_V11_UL2016C_v2/210419_184043/0000/output_0_*"),
    # '16D': glob.glob(storage+"_V11_UL2016D_v2/210419_183904/0000/output_0_*"),
    # '16E': glob.glob(storage+"_V11_UL2016E/210417_164907/0000/output_0_*"),
    # '16F': glob.glob(storage+"_V11_ULpreVFP2016F_v2/210420_065401/0000/output_0_*")
    #    +glob.glob(storage+"_V11_ULpostVFP2016F/210417_164957/0000/output_0_*"),
    # '16G': glob.glob(storage+"_V11_UL2016G_v3/210419_183936/0000/output_0_*"),
    # '16H': glob.glob(storage+"_V11_UL2016H/210417_165023/0000/output_0_*"),
    #
    # '17B': glob.glob(storage+"_V11_UL2017B/210420_072555/0000/output_0_*"),
    # '17C': glob.glob(storage+"_V11_UL2017C/210420_072607/0000/output_0_*"),
    # '17D': glob.glob(storage+"_V09_UL2017D/201214_134621/0000/output_0_*"),
    # '17E': glob.glob(storage+"_V11_UL2017E/210526_151826/0000/output_0_*"),
    # '17F': glob.glob(storage+"_V11_UL2017F/210526_151850/0000/output_0_*"),
    # '17H': glob.glob(storage+"_V11_UL2017H/210420_072656/0000/output_0_*"),
    #
    # '18A': glob.glob(storage+"_V11_UL2018A/210420_072218/0000/output_0_*"),
    # '18B': glob.glob(storage+"_V11_UL2018B/210420_072153/0000/output_0_*"),
    # '18C': glob.glob(storage+"_V11_UL2018C/210420_072235/0000/output_0_*"),
    '18D': glob.glob(storage+"_V11_UL2018D/210420_072249/000?/output_0_*"),
    # Missing runs from 2018
    '18D_2': glob.glob(storage+"_V11_UL2018D_v3/210506_184017/0000/output_0_*")
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

    df['q'] = (df['q1'] + df['q2']) / 2.

    # ### additional requirement at Sel:
    # # case for is2HLT and only tag fails additional cut: tag and probe has to be chanched
    # selectedOr = df['is2HLT']
    #
    # if tkIsoCut is not None:
    #     selectedOr = selectedOr | (df['tkIso1'] >= tkIsoCut) & (df['tkIso2'] < tkIsoCut)
    # if pfIsoCut is not None:
    #     selectedOr = selectedOr | ((df['pfIso1'] >= pfIsoCut) & (df['pfIso2'] < pfIsoCut))
    # if dxyCut is not None:
    #     selectedOr = selectedOr | ((df['dxy1'] >= dxyCut) & (df['dxy2'] < dxyCut))
    # if dzCut is not None:
    #     selectedOr = selectedOr | ((df['dz1'] >= dzCut) & (df['dz2'] < dzCut))
    #
    # selected = df['is2HLT'] & selectedOr
    #
    # to_switch = df[selected]
    # df = df[selected==False]
    # to_switch = to_switch.rename(index=str, columns={
    #     'tkIso1': 'tkIso2', 'tkIso2': 'tkIso1',
    #     'pfIso1': 'pfIso2', 'pfIso2': 'pfIso1',
    #     'pt1': 'pt2', 'pt2': 'pt1',
    #     'eta1': 'eta2', 'eta2': 'eta1',
    #     'dxy1': 'dxy2', 'dxy2': 'dxy1',
    #     'dz1': 'dz2', 'dz2': 'dz1'
    #     })
    # df = pd.concat([df,to_switch],sort=True)
    # # downgrade is2HLT and isSel to isGlo
    # df['isGlo'] = df['isGlo'] \
    #     + df['is2HLT'] * (df['tkIso2'] >= tkIsoCut) + df['isSel'] * (df['tkIso2'] >= tkIsoCut) \
    #     + df['is2HLT'] * (df['pfIso2'] >= pfIsoCut) + df['isSel'] * (df['pfIso2'] >= pfIsoCut) \
    #     + df['is2HLT'] * (df['dxy2'] >= dxyCut) + df['isSel'] * (df['dxy2'] >= dxyCut) \
    #     + df['is2HLT'] * (df['dz2'] >= dzCut) + df['isSel'] * (df['dz2'] >= dzCut)
    # df['is2HLT'] = df['is2HLT'] \
    #     * (df['tkIso2'] < tkIsoCut) \
    #     * (df['pfIso2'] < pfIsoCut) \
    #     * (df['dxy2'] < dxyCut) \
    #     * (df['dz2'] < dzCut)
    # df['isSel'] = df['isSel'] \
    #     * (df['tkIso2'] < tkIsoCut) \
    #     * (df['pfIso2'] < pfIsoCut) \
    #     * (df['dxy2'] < dxyCut) \
    #     * (df['dz2'] < dzCut)

    for q, data in df.groupby('q'):

        if q == -1:
            suffix="_Minus"
        elif q == 1:
            suffix="_Plus"
        else:
            suffix=""

        for run, data_run in data.groupby('run'):
            print(">>> run {0}".format(run))
            if max(data_run['ls']) > LumiMax_:
                print("WARNING: this run has {0} lumisections, where the maximum bin is {1}".format(max(data_run['ls']), LumiMax_))

            subdir = output+"/000"+str(run)[:-2]+"xx"
            if not os.path.isdir(subdir):
                os.mkdir(subdir)

            tfile = ROOT.TFile.Open(subdir+"/DQM_V0001_R000{0}__SingleMuon__Run20{1}.root".format(run, era),"UPDATE")

            # HLT to Sel
            h_mass_HLT_pass_central = ROOT.TH2D("h_mass_HLT_pass_central"+suffix, "Muon HLT passing probes central",
                LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
            h_mass_HLT_pass_forward = ROOT.TH2D("h_mass_HLT_pass_forward"+suffix, "Muon HLT passing probes forward",
                LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
            h_mass_HLT_fail_central = ROOT.TH2D("h_mass_HLT_fail_central"+suffix, "Muon HLT failing probes central",
                LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
            h_mass_HLT_fail_forward = ROOT.TH2D("h_mass_HLT_fail_forward"+suffix, "Muon HLT failing probes forward",
                LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

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

            h_mass_yieldBB_Z = ROOT.TH2D("h_mass_yieldBB_Z"+suffix, "reconstructed Z bosons, both muons in barrel",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
            h_mass_yieldBE_Z = ROOT.TH2D("h_mass_yieldBE_Z"+suffix, "reconstructed Z bosons, muons in barrel and endcap", LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);
            h_mass_yieldEE_Z = ROOT.TH2D("h_mass_yieldEE_Z"+suffix, "reconstructed Z bosons, both muons in endcap",       LumiBin_, LumiMin_, LumiMax_, MassBin_, MassMin_, MassMax_);

            #fill tag and probe histograms
            #--- barrel
            # fill tags if both muons passed HLT (2HLT category)
            for m,ls in data_run.query("is2HLT==1 & abs(eta1) < 0.9")[['dilepMass','ls']].values:
                h_mass_HLT_pass_central.Fill(ls,m)
                h_mass_SEL_pass_central.Fill(ls,m)
                # h_mass_GLO_pass_central.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_central.Fill(ls,m)
                h_mass_TRK_pass_central.Fill(ls,m)
                h_mass_STA_pass_central.Fill(ls,m)

            # fill probes
            for m,ls,eta1 in data_run.query("is2HLT==1 & abs(eta2) < 0.9")[['dilepMass','ls','eta1']].values:
                h_mass_HLT_pass_central.Fill(ls,m)
                h_mass_SEL_pass_central.Fill(ls,m)
                # h_mass_GLO_pass_central.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_central.Fill(ls,m)
                h_mass_TRK_pass_central.Fill(ls,m)
                h_mass_STA_pass_central.Fill(ls,m)
                if(abs(eta1) < 0.9):
                    h_mass_yieldBB_Z.Fill(ls,m)
                else:
                    h_mass_yieldBE_Z.Fill(ls,m)

            for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) < 0.9")[['dilepMass','ls','eta1']].values:
                h_mass_HLT_fail_central.Fill(ls,m)
                h_mass_SEL_pass_central.Fill(ls,m)
                # h_mass_GLO_pass_central.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_central.Fill(ls,m)
                h_mass_TRK_pass_central.Fill(ls,m)
                h_mass_STA_pass_central.Fill(ls,m)
                if(abs(eta1) < 0.9):
                    h_mass_yieldBB_Z.Fill(ls,m)
                else:
                    h_mass_yieldBE_Z.Fill(ls,m)

            for m,ls in data_run.query("isGlo==1 & abs(eta2) < 0.9")[['dilepMass','ls']].values:
                h_mass_SEL_fail_central.Fill(ls,m)
                # h_mass_GLO_pass_central.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_central.Fill(ls,m)
                h_mass_STA_pass_central.Fill(ls,m)
                h_mass_TRK_pass_central.Fill(ls,m)

            for m,ls in data_run.query("isTrk==1 & abs(eta2)  < 0.9")[['dilepMass','ls']].values:
                # h_mass_GLO_fail_central.Fill(ls,m)
                h_mass_STA_fail_central.Fill(ls,m)
                h_mass_TRK_pass_central.Fill(ls,m)

            for m,ls in data_run.query("isSta==1 & abs(eta2)  < 0.9")[['dilepMass','ls']].values:
                h_mass_TRK_fail_central.Fill(ls,m)
                h_mass_STA_pass_central.Fill(ls,m)
                # h_mass_GLOtoSTA_fail_central.Fill(ls,m)

            #--- endcap
            # fill tags if both muons passed HLT (2HLT category)
            for m,ls in data_run.query("is2HLT==1 & abs(eta1) >= 0.9")[['dilepMass','ls']].values:
                h_mass_HLT_pass_forward.Fill(ls,m)
                h_mass_SEL_pass_forward.Fill(ls,m)
                # h_mass_GLO_pass_forward.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_forward.Fill(ls,m)
                h_mass_TRK_pass_forward.Fill(ls,m)
                h_mass_STA_pass_forward.Fill(ls,m)

            # fill probes
            for m,ls,eta1 in data_run.query("is2HLT==1 & abs(eta2) >= 0.9")[['dilepMass','ls','eta1']].values:
                h_mass_HLT_pass_forward.Fill(ls,m)
                h_mass_SEL_pass_forward.Fill(ls,m)
                # h_mass_GLO_pass_forward.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_forward.Fill(ls,m)
                h_mass_TRK_pass_forward.Fill(ls,m)
                h_mass_STA_pass_forward.Fill(ls,m)
                if(abs(eta1) > 0.9):
                    h_mass_yieldEE_Z.Fill(ls,m)
                else:
                    h_mass_yieldBE_Z.Fill(ls,m)

            for m,ls,eta1 in data_run.query("isSel==1 & abs(eta2) >= 0.9")[['dilepMass','ls','eta1']].values:
                h_mass_HLT_fail_forward.Fill(ls,m)
                h_mass_SEL_pass_forward.Fill(ls,m)
                # h_mass_GLO_pass_forward.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_forward.Fill(ls,m)
                h_mass_TRK_pass_forward.Fill(ls,m)
                h_mass_STA_pass_forward.Fill(ls,m)
                if(abs(eta1) > 0.9):
                    h_mass_yieldEE_Z.Fill(ls,m)
                else:
                    h_mass_yieldBE_Z.Fill(ls,m)


            for m,ls in data_run.query("isGlo==1 & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
                h_mass_SEL_fail_forward.Fill(ls,m)
                # h_mass_GLO_pass_forward.Fill(ls,m)
                # h_mass_GLOtoSTA_pass_forward.Fill(ls,m)
                h_mass_TRK_pass_forward.Fill(ls,m)
                h_mass_STA_pass_forward.Fill(ls,m)

            for m,ls in data_run.query("isTrk==1 & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
                # h_mass_GLOtoTRK_fail_forward.Fill(ls,m)
                h_mass_TRK_pass_forward.Fill(ls,m)
                h_mass_STA_fail_forward.Fill(ls,m)

            for m,ls in data_run.query("isSta==1 & abs(eta2) >= 0.9")[['dilepMass','ls']].values:
                h_mass_TRK_fail_forward.Fill(ls,m)
                h_mass_STA_pass_forward.Fill(ls,m)
                # h_mass_GLOtoSTA_fail_forward.Fill(ls,m)

            tfile.Write()
            tfile.Close()
            tfile.Delete()
