from __future__ import division, print_function

import argparse
import glob
import pdb
import os

parser = argparse.ArgumentParser(prog='./histMerge')
parser.add_argument(
    '-i','--input', required=True, type=str,
    help='specify input directory'
)

args = parser.parse_args()

for runDir in glob.glob(args.input+"/000*"):


    filelist = glob.glob(runDir+"/DQM_V0001_R000*")
    print("Found "+str(len(filelist))+" files in "+runDir+" to hadd")

    os.system("hadd {0}/merged.root ".format(runDir)+" ".join(filelist))

filelist = glob.glob(args.input+"/000*/merged.root")
print("Found "+str(len(filelist))+" merged files to hadd")

os.system("hadd {0}/merged.root ".format(args.input)+" ".join(filelist))

print("remove temporary merged files")
for file in filelist:
    os.system("rm {0}".format(file))

print("Done!")
