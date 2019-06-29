import os,sys
#sys.argv.append( '-b-' )
#from scipy.interpolate import UnivariateSpline
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
#linesatlas=atlasfile.readlines()
import StringIO

with open('CMSZrates_Inclusive.csv', 'w') as f:
    for linecms in range(0,len(linescms)):
	###print linescms[linecms]
	elements=linescms[linecms].split(",")
	line=elements[0]+","+elements[1]+","+elements[2]+","+str(1.13*float(elements[3]))+","+elements[4]+","+elements[5]+","+str(1.13*float(elements[6]))
	#print elements[0]+"  "+str(elements[3])
	
        f.write(line)
	f.write('\n')


