import argparse
from python.logging import child_logger
log = child_logger(__name__)

# basic signal models (index for root)
sigModels = {
    "BWxCB":1,
    "MCxGauss": 2,
    "BW":3,
    "MC":4,
    "BWxGauss":5,
    "MCxBW":6,
    "Gen":16
}

# basic background models (index for root)
bkgModels = {
    "Const":0,
    "Exp":1,
    "Quad": 2,
    "QuadPlusExp":3,
    "Das":4,
    "CMSShape":6,
    "default":None,
    "alt":None
}

# setting default and alternative choices (differnt for passing and failing categories)
bkgModelsPass = bkgModels
bkgModelsPass.update({
    "default": bkgModels["Exp"],
    "alt": bkgModels["Const"]
})

bkgModelsFail = bkgModels
bkgModelsFail.update({
    "default": bkgModels["CMSShape"],
    "alt": bkgModels["Das"]
})

def parser():
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4],
                        help="Set verbosity level with logging, the larger the more verbose")

    return parser

def parser_zharvest(parser):

    parser.add_argument("-b", "--beginRun", type=int, default=272007, 
                        help="first run to analyze [%(default)s]")
    parser.add_argument("-e", "--endRun", type=int, default=1000000, 
                        help="analyze stops when comes to this run [%(default)s]")
    parser.add_argument('--mcCorrections', default="default", type=str,
                        help='specify .json file with MC corrections for muon correlations')
    parser.add_argument("-c", "--writeSummaryCSV", default=False, action="store_true",
                        help="produce merged CSV with all runs")
    parser.add_argument("-i", "--input", default="/eos/cms/store/group/comm_luminosity/ZCounting/2022/DQMFiles/cmsweb.cern.ch/dqm/offline/data/browse/ROOT/OfflineData/Run2022/", 
                        help="Directory to the input root files from the DQM Offline module")
    parser.add_argument("--byLsCSV", default="/eos/cms/store/group/comm_luminosity/ZCounting/2022/brilcalcByLS/byLS_Collisions22_355100_362760_Muon_20230210.csv", 
                        help="ByLs csv input generated by testBril.sh")
    parser.add_argument("--sigModel", default="MCxGauss", type=str, choices=sigModels.keys(),
                        help="Choose one of the options for signal model.")
    parser.add_argument("--bkgModel", default="default", type=str, choices=bkgModels.keys(),
                        help="Choose one of the options for background model.")
    parser.add_argument('--ptCut', type=float, default=25.,
                        help='specify lower pt cut on tag and probe muons')
    parser.add_argument('--etaCut', type=float, default=2.4,
                        help='specify upper |eta| cut on tag and probe muons')
    parser.add_argument('--mass', nargs=3, metavar=('LOW', 'HIGH', 'NUMBER'), default=(60,120,120), type=int,
                        help='specify mass range for tag and probe muon pairs')
    parser.add_argument('--LumiPerMeasurement', default=20, type=float,
                        help='specify amount of luminosity per measurement in pb-1')
    parser.add_argument('--collect', default=False, action="store_true",
                        help='specify whether or not to run the fits or just collect the results')
    parser.add_argument("-o", "--output", default="./", 
                        help="where to store the output files")

    return parser

def set_parser_default(parser, argument, newDefault):
    # change the default argument of the parser, must be called before parse_arguments
    f = next((x for x in parser._actions if x.dest ==argument), None)
    if f:
        log.info(f" Modifying default of {f.dest} from {f.default} to {newDefault}")
        f.default = newDefault
    else:
        log.warning(f" Parser argument {argument} not found!")

