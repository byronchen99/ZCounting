from __future__ import division, print_function

# production of histograms from TnP trees
#   the histograms are for data driven background shapes in nPV and
#   and can therefore used with the same toopu, e.g. ZCounting.py

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

PUMin_ = -0.5
PUMax_ = 100.5
PUBin_ = int(PUMax_ - PUMin_)

ptCut = args.pt if args.pt else 0.

dxyDistCut = None# 0.2
dzDistCut = None#0.5


# acceptance selection
selection = 'q1==q2 & pt1 >= {0} & pt2 >= {0} & dilepMass >= {1} &  dilepMass <= {2}'.format(ptCut, MassMin_, MassMax_)

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
            'run',
            'eta1', 'eta2',
            # 'pt1', 'pt2',
            'is2HLT', 'isSel', 'isGlo', 'isSta', 'isTrk',
            # 'tkIso1', 'tkIso2', 'pfIso1', 'pfIso2',
            # 'dxy1', 'dxy2', 'dz1', 'dz2',
            # 'q1','q2'
            ]

storage="/pnfs/desy.de/cms/tier2/store/user/dwalter/SingleMuon/TnPPairTrees"

inPVts = {

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
    # '17E': glob.glob(storage+"_V09_UL2017E/201214_134634/0000/output_0_*"),
    # '17F': glob.glob(storage+"_V09_UL2017F/201214_141544/0000/output_0_*"),
    '17H': glob.glob(storage+"_V11_UL2017H/210420_072656/0000/output_0_*"),

    # '18A': glob.glob(storage+"_V11_UL2018A/210420_072218/0000/output_0_*"),
    # '18B': glob.glob(storage+"_V11_UL2018B/210420_072153/0000/output_0_*"),
    # '18C': glob.glob(storage+"_V11_UL2018C/210420_072235/0000/output_0_*"),
    # '18D': glob.glob(storage+"_V11_UL2018D/210420_072249/000?/output_0_*"),
    # # Missing runs from 2018
    # '18D': glob.glob(storage+"_V11_UL2018D_v3/210506_184017/0000/output_0_*")
}



for era, inPVt in inPVts.iteritems():
    print("##### era {0}".format(era))

    treeName = rn.list_trees(inPVt[0])[0]
    print(">>> Load Events from {0} files".format(len(inPVt)))
    _df = []
    for i in inPVt:
        print("> file "+i)
        tfile = ROOT.TFile.Open(i)
        if(tfile.Get(treeName).GetEntries(selection) == 0):
            print("> no events in this file found! continue with next file")
            continue
        _df.append(tree_to_df(rn.root2array(i, treeName, selection=selection, branches=branches)))
    print(">>> Concatenate")
    df = pd.concat(_df)

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

    tfile = ROOT.TFile.Open(output+"/Run20{0}.root".format(era),"RECREATE")

    # HLT to Sel
    h_mass_HLT_pass_central = ROOT.TH2D("h_mass_HLT_pass_central", "Muon HLT passing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_HLT_pass_forward = ROOT.TH2D("h_mass_HLT_pass_forward", "Muon HLT passing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_HLT_fail_central = ROOT.TH2D("h_mass_HLT_fail_central", "Muon HLT failing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_HLT_fail_forward = ROOT.TH2D("h_mass_HLT_fail_forward", "Muon HLT failing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);

    # Sel to Glo
    h_mass_SEL_pass_central = ROOT.TH2D("h_mass_Sel_pass_central", "Muon Sel passing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_SEL_pass_forward = ROOT.TH2D("h_mass_Sel_pass_forward", "Muon Sel passing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_SEL_fail_central = ROOT.TH2D("h_mass_Sel_fail_central", "Muon Sel failing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_SEL_fail_forward = ROOT.TH2D("h_mass_Sel_fail_forward", "Muon Sel failing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);

    # Trk to Sta
    h_mass_TRK_pass_central = ROOT.TH2D("h_mass_Trk_pass_central", "Muon track passing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_TRK_pass_forward = ROOT.TH2D("h_mass_Trk_pass_forward", "Muon track passing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_TRK_fail_central = ROOT.TH2D("h_mass_Trk_fail_central", "Muon track failing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_TRK_fail_forward = ROOT.TH2D("h_mass_Trk_fail_forward", "Muon track failing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);

    # Sta to Trk
    h_mass_STA_pass_central = ROOT.TH2D("h_mass_Sta_pass_central", "Muon standalone passing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_STA_pass_forward = ROOT.TH2D("h_mass_Sta_pass_forward", "Muon standalone passing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_STA_fail_central = ROOT.TH2D("h_mass_Sta_fail_central", "Muon standalone failing probes central",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_STA_fail_forward = ROOT.TH2D("h_mass_Sta_fail_forward", "Muon standalone failing probes forward",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);

    h_mass_yieldBB_Z = ROOT.TH2D("h_mass_yieldBB_Z", "reconstructed Z bosons, both muons in barrel",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_yieldBE_Z = ROOT.TH2D("h_mass_yieldBE_Z", "reconstructed Z bosons, muons in barrel and endcap",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);
    h_mass_yieldEE_Z = ROOT.TH2D("h_mass_yieldEE_Z", "reconstructed Z bosons, both muons in endcap",
        PUBin_, PUMin_, PUMax_, MassBin_, MassMin_, MassMax_);


    print(">>> Fill histograms")
    #fill tag and probe histograms
    #--- barrel
    # fill tags if both muons passed HLT (2HLT category)
    for m,pu in df.query("is2HLT==1 & abs(eta1) < 0.9")[['dilepMass','nPV']].values:
        h_mass_HLT_pass_central.Fill(pu,m)
        h_mass_SEL_pass_central.Fill(pu,m)
        h_mass_TRK_pass_central.Fill(pu,m)
        h_mass_STA_pass_central.Fill(pu,m)

    # fill probes
    for m,pu,eta1 in df.query("is2HLT==1 & abs(eta2) < 0.9")[['dilepMass','nPV','eta1']].values:
        h_mass_HLT_pass_central.Fill(pu,m)
        h_mass_SEL_pass_central.Fill(pu,m)
        h_mass_TRK_pass_central.Fill(pu,m)
        h_mass_STA_pass_central.Fill(pu,m)
        if(abs(eta1) < 0.9):
            h_mass_yieldBB_Z.Fill(pu,m)
        else:
            h_mass_yieldBE_Z.Fill(pu,m)

    for m,pu,eta1 in df.query("isSel==1 & abs(eta2) < 0.9")[['dilepMass','nPV','eta1']].values:
        h_mass_HLT_fail_central.Fill(pu,m)
        h_mass_SEL_pass_central.Fill(pu,m)
        h_mass_TRK_pass_central.Fill(pu,m)
        h_mass_STA_pass_central.Fill(pu,m)
        if(abs(eta1) < 0.9):
            h_mass_yieldBB_Z.Fill(pu,m)
        else:
            h_mass_yieldBE_Z.Fill(pu,m)

    for m,pu in df.query("isGlo==1 & abs(eta2) < 0.9")[['dilepMass','nPV']].values:
        h_mass_SEL_fail_central.Fill(pu,m)
        h_mass_STA_pass_central.Fill(pu,m)
        h_mass_TRK_pass_central.Fill(pu,m)

    for m,pu in df.query("isTrk==1 & abs(eta2)  < 0.9")[['dilepMass','nPV']].values:
        h_mass_STA_fail_central.Fill(pu,m)
        h_mass_TRK_pass_central.Fill(pu,m)

    for m,pu in df.query("isSta==1 & abs(eta2)  < 0.9")[['dilepMass','nPV']].values:
        h_mass_TRK_fail_central.Fill(pu,m)
        h_mass_STA_pass_central.Fill(pu,m)

    #--- endcap
    # fill tags if both muons passed HLT (2HLT category)
    for m,pu in df.query("is2HLT==1 & abs(eta1) >= 0.9")[['dilepMass','nPV']].values:
        h_mass_HLT_pass_forward.Fill(pu,m)
        h_mass_SEL_pass_forward.Fill(pu,m)
        h_mass_TRK_pass_forward.Fill(pu,m)
        h_mass_STA_pass_forward.Fill(pu,m)

    # fill probes
    for m,pu,eta1 in df.query("is2HLT==1 & abs(eta2) >= 0.9")[['dilepMass','nPV','eta1']].values:
        h_mass_HLT_pass_forward.Fill(pu,m)
        h_mass_SEL_pass_forward.Fill(pu,m)
        h_mass_TRK_pass_forward.Fill(pu,m)
        h_mass_STA_pass_forward.Fill(pu,m)
        if(abs(eta1) > 0.9):
            h_mass_yieldEE_Z.Fill(pu,m)
        else:
            h_mass_yieldBE_Z.Fill(pu,m)

    for m,pu,eta1 in df.query("isSel==1 & abs(eta2) >= 0.9")[['dilepMass','nPV','eta1']].values:
        h_mass_HLT_fail_forward.Fill(pu,m)
        h_mass_SEL_pass_forward.Fill(pu,m)
        h_mass_TRK_pass_forward.Fill(pu,m)
        h_mass_STA_pass_forward.Fill(pu,m)
        if(abs(eta1) > 0.9):
            h_mass_yieldEE_Z.Fill(pu,m)
        else:
            h_mass_yieldBE_Z.Fill(pu,m)


    for m,pu in df.query("isGlo==1 & abs(eta2) >= 0.9")[['dilepMass','nPV']].values:
        h_mass_SEL_fail_forward.Fill(pu,m)
        h_mass_TRK_pass_forward.Fill(pu,m)
        h_mass_STA_pass_forward.Fill(pu,m)

    for m,pu in df.query("isTrk==1 & abs(eta2) >= 0.9")[['dilepMass','nPV']].values:
        h_mass_TRK_pass_forward.Fill(pu,m)
        h_mass_STA_fail_forward.Fill(pu,m)

    for m,pu in df.query("isSta==1 & abs(eta2) >= 0.9")[['dilepMass','nPV']].values:
        h_mass_TRK_fail_forward.Fill(pu,m)
        h_mass_STA_pass_forward.Fill(pu,m)

    tfile.Write()
    tfile.Close()
    tfile.Delete()
