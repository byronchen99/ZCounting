import os,sys
import StringIO
import ROOT
from ROOT import TGraph
from ROOT import TGraphErrors
from array import array
from ROOT import *
from operator import truediv
import random
import argparse
import scipy.integrate as integrate
ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetCanvasPreferGL(1)

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cms", default="nothing", type=string, help="give the CMS csv as input")
args = parser.parse_args()
if args.cms=="nothing":
	print "please provide cms input files"
	sys.exit()
print args.cms

cmsfile=open(str(args.cms))

linescms=cmsfile.readlines()

with open('CMSZrates_timestampSwap_Inclusive.csv', 'w') as f:
    for linecms in range (0,len(linescms)):
	elements=linescms[linecms].split(",")
	interElOne=elements[1].split(" ")
	newElOne=" "+str(interElOne[0].split("/")[2])+"/"+str(interElOne[0].split("/")[0])+"/"+interElOne[0].split("/")[1]+" "+interElOne[1]
	interElTwo=elements[2].split(" ")
	newElTwo=" "+str(interElTwo[0].split("/")[2])+"/"+str(interElTwo[0].split("/")[0])+"/"+interElTwo[0].split("/")[1]+" "+interElTwo[1]
	
	print elements[2]
	print newElTwo
	line=elements[0]+","+newElOne+","+newElTwo+","+elements[3]+","+elements[4]+","+elements[5]+","+elements[6]
        f.write(line)
#	f.write('\n')


