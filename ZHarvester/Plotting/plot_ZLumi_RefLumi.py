import argparse
import pandas as pd
import uncertainties as unc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from datetime import datetime
import os, sys
import pdb

sys.path.append(os.getcwd())
print(os.getcwd())

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from python.corrections import apply_muon_prefire, apply_ECAL_prefire
from python.utils import to_DateTime, load_input_csv
from ZUtils.python.utils import to_RootTime

# disable panda warnings when assigning a new column in the dataframe
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--rates", required=True, type=str, help="csv file with z rates per measurement")
parser.add_argument("-l", "--refLumi", default="", type=str, help="give a ByLs.csv as input for additional reference Luminosity")
parser.add_argument("-x", "--xsec", default="", type=str, help="csv file with z rates per measurement for absolute scale")
parser.add_argument("-f", "--fill", nargs="*",  type=int, default=[], help="specify a single fill to plot")
parser.add_argument("-y", "--year",  default=2017, type=int, help="give a year for calculation of time")
parser.add_argument("-s", "--saveDir",  default='./',  type=str, help="give output dir")
args = parser.parse_args()
outDir = args.saveDir
if not os.path.isdir(outDir):
    os.mkdir(outDir)

########## Settings ##########

currentYear = args.year
secPerLS=float(23.3)

fmt = "png"

# plotting options

textsize = 16
markersize = 4.0
lefttitle = "$\sqrt{s}=13\,\mathrm{TeV}$"
xlabel = "LHC runtime [h]"
ylabelLumi = "Inst. luminosity [nb$^{-1}$s$^{-1}$]"
ylabelEff = "Efficiency"

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize,
    'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})

    
########## Data Acquisition ##########

if args.refLumi != "":
    # --- aditional external reference luminosity
    extLumi = True 

    data_ref = load_input_csv(str(args.refLumi))

    # # convert reference data from ByLs.csv file
    # lumiFile=open(str(args.refLumi))
    # lumiLines=lumiFile.readlines()
    # data_ref = pd.read_csv(str(args.refLumi), sep=',',low_memory=False, skiprows=lambda x: lumiLines[x].startswith('#') and not lumiLines[x].startswith('#run'))
    # if 'recorded(/ub)' in data_ref.columns.tolist():      #convert to /pb
    #     data_ref['recorded(/ub)'] = data_ref['recorded(/ub)'].apply(lambda x:x / 1000000.)
    #     data_ref = data_ref.rename(index=str, columns={'recorded(/ub)':'recorded(/pb)' })
    # elif 'recorded(/fb)' in data_ref.columns.tolist():      #convert to /pb
    #     data_ref['recorded(/fb)'] = data_ref['recorded(/fb)'].apply(lambda x:x * 1000.)
    #     data_ref = data_ref.rename(index=str, columns={'recorded(/fb)':'recorded(/pb)' })

    # data_ref['fill'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[1])
    # data_ref['run'] = pd.to_numeric(data_ref['#run:fill'].str.split(':',expand=True)[0])
    # data_ref['ls'] = pd.to_numeric(data_ref['ls'].str.split(':',expand=True)[0])
    # data_ref = data_ref.drop(['#run:fill','hltpath','source'],axis=1)
    # data_ref = data_ref.sort_values(['fill','run','ls','recorded(/pb)'])
    # data_ref = data_ref.drop_duplicates(['fill','run','ls'])

    if args.fill != []:
        data_ref = data_ref.loc[data_ref['fill'].isin(args.fill)]
    else:
        print("Plot all fills")

    data_ref['time'] = data_ref['time'].apply(lambda x: to_DateTime(x))
    # data_ref['timeE'] = data_ref['time'] - data_ref['time']

    data_ref['recorded(/nb)'] = data_ref['recorded(/pb)'] * 1000  # convert into /nb 
    data_ref['dLRec(/nb)'] = data_ref['recorded(/nb)']/secPerLS
    data_ref = data_ref[['fill','dLRec(/nb)','time', 'recorded(/nb)', 'avgpu']]		#Keep only what you need

else:
    extLumi = False 
  

# --- absolute scale for z luminosity
if args.xsec != "":
    print("get Z cross section")
    data_xsec = pd.read_csv(str(args.xsec), sep=',',low_memory=False)#, skiprows=[1,2,3,4,5])
    data_xsec['zDelI'] = data_xsec['zDel'].apply(lambda x: unc.ufloat_fromstr(x).n)

    apply_muon_prefire(data_xsec)
    apply_ECAL_prefire(data_xsec)

    print("apply prefire corrections - done")

    xsecI = sum(data_xsec['zDelI'])/sum(data_xsec['lumiRec'])
else:
    xsecI = 1


# --- z luminosity
data = pd.read_csv(str(args.rates), sep=',',low_memory=False) #, skiprows=[1,2,3,4,5])
if args.fill != []:
    data = data.loc[data['fill'].isin(args.fill)]

data['tdate_begin'] = data['beginTime'].apply(lambda x: to_RootTime(x,currentYear)).astype(float)
data['tdate_end'] = data['endTime'].apply(lambda x: to_RootTime(x,currentYear)).astype(float)

data = data.sort_values(['fill','tdate_begin','tdate_end'])

data['tdate_begin'] = data['tdate_begin']# - 31536000
data['tdate_end'] = data['tdate_end']# - 31536000

# origin of ROOT TDatime
data['time'] = data['tdate_begin'] + (data['tdate_end'] - data['tdate_begin'])/2

data['time'] = data['time'].apply(datetime.utcfromtimestamp)
data['timeUp'] = data['tdate_end'].apply(datetime.utcfromtimestamp)
data['timeDown'] = data['tdate_begin'].apply(datetime.utcfromtimestamp)

def unorm(x):
    # for counting experiments: define ufloat with poisson uncertainty
    return unc.ufloat(x, np.sqrt(abs(x)))

# sort out nan values
nan = np.isnan(data['ZRate'])
print(f"Sort out {sum(nan)} nan values")
data = data.loc[nan==False]

# efficiency corrected number of Z bosons in timewindow, w/o deadtime corrections
data['zDelI_mc'] = data['ZRate'] * data['timewindow'] * data['deadtime']

data['zDelI_mc'] = data['zDelI_mc'].apply(lambda x: unorm(x))

data['zLumi_mc'] = data['zDelI_mc'] / xsecI
data['zLumiInst_mc'] = data['zLumi_mc'] / data['timewindow'] * 1000  # convert into /nb
# reference lumi during one measurement
data['dLRec(/nb)'] = data['recLumi'] / data['timewindow'] * 1000  # convert into /nb 

# # convert inputs to uncertainties
# for key in ('HLTeff', 'Seleff', 'Gloeff', 'Staeff'):
#     data[key] = data[key].apply(lambda x: unc.ufloat_fromstr(x))


do_ratio=False

for fill, data_fill in data.groupby("fill"):
    print(f"Now at fill {fill}")

    
    if len(data_fill) == 1:
        print("Only one measurement, next fill ")
        continue

    if extLumi:
        # group external reference luminosity
        # # only take range where we measured Zs
        # ref_fill = ref_fill.loc[(ref_fill['time'] >= data_fill['timeDown'].values[0]) & (ref_fill['time'] < data_fill['timeUp'].values[-1])]
        # 

        ref_fill = data_ref.loc[data_ref["fill"] == fill]

        if len(ref_fill) == 0:
            print("WARNING! Fill is not present in external luminosity source!")

        ref_fill_lumi = []
        for tUp, tDown in data_fill[['timeUp','timeDown']].values:
            ref_fill_lumi.append(ref_fill.loc[(ref_fill['time'] >= tDown) & (ref_fill['time'] < tUp)]['recorded(/nb)'].sum() / ((mpl.dates.date2num(tUp) - mpl.dates.date2num(tDown)) *24 * 3600) )
        
        yRefExt = np.array(ref_fill_lumi)

    x = mpl.dates.date2num(data_fill['time'].values)
    xUp = mpl.dates.date2num(data_fill['timeUp'].values)
    xDown = mpl.dates.date2num(data_fill['timeDown'].values)

    # convert into hours
    xUp = (xUp - x) * 24 
    xDown = (x - xDown) * 24
    x = x * 24
    x = (x - x[0] + xDown[0])
    
    # x axis range
    xMin = min(x-xDown)
    xMax = max(x+xUp)
    xRange = xMax - xMin
    xMin = xMin - xRange * 0.025
    xMax = xMax + xRange * 0.025
    xRange = xMax - xMin

    if int(max(x)) > 8:
        ticksteps = 2
    else:
        ticksteps = 1
        
    xTicks = np.arange(0, int(max(x))+ticksteps, ticksteps)

    # average pileup as x axis
    xPU = data_fill['pileUp'].values
    xMinPU = min(xPU)
    xMaxPU = max(xPU)   
    xRangePU = xMaxPU - xMinPU 
    xMinPU = xMinPU - xRangePU * 0.025
    xMaxPU = xMaxPU + xRangePU * 0.025
    xRangePU = xMaxPU - xMinPU

    if int(max(xPU)) > 8:
        ticksteps = 2
    else:
        ticksteps = 1
        
    xTicksPU = np.arange(0, int(max(xPU))+ticksteps, ticksteps)
    
    ### Make plot for efficiency vs time
    plt.clf()
    fig = plt.figure()
    
    # ax1 = fig.add_subplot(111)
    
    gs = gridspec.GridSpec(2, 1)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.125, hspace=0.0)

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabelEff)
        
    maxY = 0
    minY = 1
    for eff, name, col in (
        # ("effGlo", "$ \epsilon_\mathrm{Glo|Sta}^\mu $", "green"), 
        # ("effSta", "$ \epsilon_\mathrm{Sta|Trk}^\mu $", "blue"), 
        ("effSel", "$ \epsilon_\mathrm{ID|Glo}^\mu $", "red"), 
    ):
        
        # y = [y.n for y in data_fill[eff].values]
        # yErr = [y.s for y in data_fill[eff].values]

        y = data_fill[eff].values

        ax1.errorbar(x, y, xerr=(xDown, xUp), #yerr=yErr, 
            label=name,
            marker="o", linewidth=0, color=col, ecolor=col, elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)

        # ax1.plot(x, y, label=name,
        #     marker="o", linewidth=0, color=col, markersize=markersize)

        # yMax = [y.n + y.s for y in data_fill[eff].values]
        # yMin = [y.n - y.s for y in data_fill[eff].values]
                    
        maxY = max(maxY, max(y))
        minY = min(minY, min(y))

    leg = ax1.legend(loc="upper left", ncol=3,
        frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    leg.get_frame().set_linewidth(0.8)
    
    yRange = maxY - minY
    ax1.set_ylim([minY-yRange*0.05, maxY + yRange*0.5])
    ax1.set_xlim([xMin, xMax])
    ax1.set_xticks(xTicks)

    ax1.xaxis.set_major_locator(ticker.NullLocator())

    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabelEff)

    ax2.text(0.54, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax2.transAxes, fontweight="bold")
    ax2.text(0.65, 0.97, "\\emph{Work in progress}", verticalalignment='top', transform=ax2.transAxes, style='italic')    
    ax2.text(0.54, 0.86, f"Fill {fill}", verticalalignment='top', transform=ax2.transAxes)    

    maxY = 0
    minY = 1
    for eff, name, col in (
        ("effHLT", "$ \epsilon_\mathrm{HLT}^\mu $", "k"), 
        # ("ZIeff", "$ ( \epsilon_\mathrm{ID}^\mu ) ^2 $", "r") 
    ):    
        # y = [y.n for y in data_fill[eff].values]
        # yErr = [y.s for y in data_fill[eff].values]

        y = data_fill[eff].values

        ax2.errorbar(x, y, xerr=(xDown, xUp), #yerr=yErr, 
            label=name,
            marker="o", linewidth=0, color=col, ecolor=col, elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)

        # ax2.plot(x, y, label=name,
        #     marker="o", linewidth=0, color=col, markersize=markersize)
        # yMax = [y.n + y.s for y in data_fill[eff].values]
        # yMin = [y.n - y.s for y in data_fill[eff].values]
                    
        maxY = max(maxY, max(y))
        minY = min(minY, min(y))

    leg = ax2.legend(loc="upper left", ncol=1, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black",
        # borderpad=1, labelspacing=1
        )    
    leg.get_frame().set_linewidth(0.8)
        
    yRange = maxY - minY
    ax2.set_ylim([minY-yRange*0.05, maxY + yRange*0.4])
    ax2.set_xlim([xMin, xMax])
    ax2.set_xticks(xTicks)

    
    # align y labels
    ax1.yaxis.set_label_coords(-0.12, 0.5)
    ax2.yaxis.set_label_coords(-0.12, 0.5)

    plt.savefig(outDir+f"/fill{fill}_eff.{fmt}")
    plt.close()    
    
    
    ### Make plot for lumi vs time
    plt.clf()
    fig = plt.figure()
    if do_ratio:
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
    else:
        ax1 = fig.add_subplot(111)
        
    fig.subplots_adjust(hspace=0.0, left=0.15, right=0.99, top=0.99, bottom=0.125)
        
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabelLumi)
    ax1.text(0.54, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax1.transAxes, weight="bold")
    ax1.text(0.65, 0.97, "\\emph{Work in progress}", verticalalignment='top', transform=ax1.transAxes,style='italic')    
    ax1.text(0.54, 0.89, f"Fill {fill}", verticalalignment='top', transform=ax1.transAxes)    
        
    if args.xsec == "":
        # normalize Z luminosity to reference luminosity    
        # data_fill['zLumiInst_mc'] *= ref_fill['recorded(/pb)'].sum() / data_fill['zLumi_mc'].sum()
        # data_fill['zLumi_mc'] *= ref_fill['recorded(/pb)'].sum() / data_fill['zLumi_mc'].sum()

        data_fill['zLumiInst_mc'] *= data_fill['recLumi'].sum() / data_fill['zLumi_mc'].sum()
        data_fill['zLumi_mc'] *= data_fill['recLumi'].sum() / data_fill['zLumi_mc'].sum()


    y = np.array([yy.n for yy in data_fill['zLumiInst_mc'].values])
    yErr = np.array([y.s for y in data_fill['zLumiInst_mc'].values])

    yRef = data_fill['dLRec(/nb)'].values   
    
    ax1.errorbar(x, yRef, xerr=(xDown, xUp), label="Reference luminosity", color="blue", 
        linestyle='', zorder=0)

    if extLumi:
        ax1.errorbar(x, yRefExt, xerr=(xDown, xUp), label="RAMSES", color="red", linestyle='',
            zorder=0)

    ax1.errorbar(x, y, xerr=(xDown, xUp), yerr=yErr, 
        label="Z rate",# measurement",
        fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)

    leg = ax1.legend(loc="lower left", ncol=2, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    leg.get_frame().set_linewidth(0.8)

    yMin = min(min(y),min(yRef))
    yMax = max(max(y),max(yRef))

    yRange = yMax - yMin 
    ax1.set_ylim([yMin - yRange*0.5, yMax + yRange*0.2])
    ax1.set_xlim([xMin, xMax])
    ax1.set_xticks(xTicks)
    
    if do_ratio:
        ax1.xaxis.set_major_locator(ticker.NullLocator())
    
        ratios = []
        for zlumi, lumi in data_fill[['zLumi_mc','recLumi']].values:
            ratios.append(zlumi/lumi)

        yRatio = np.array([yy.n for yy in ratios])
        yRatioErr = np.array([y.s for y in ratios])
    
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel("Ratio")
        
        ax2.errorbar(x, yRatio, xerr=(xDown, xUp), yerr=yRatioErr, 
            label="Z rate measurement",
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
            
        ax2.plot(np.array([xMin, xMax]), np.array([1.0, 1.0]), color="black",linestyle="-", linewidth=1)
        
        ax2.set_ylim([0.951,1.049])
        ax2.set_xlim([xMin, xMax])
        ax2.set_xticks(xTicks)

    # align y labels
    ax1.yaxis.set_label_coords(-0.12, 0.5)
    ax2.yaxis.set_label_coords(-0.12, 0.5)

    plt.savefig(outDir+f"/fill{fill}_lumi.{fmt}")
    plt.close()


    ### Make plot for lumi vs pileup
    plt.clf()
    fig = plt.figure()
    if do_ratio:
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
    else:
        ax1 = fig.add_subplot(111)
        
    fig.subplots_adjust(hspace=0.0, left=0.15, right=0.99, top=0.99, bottom=0.125)
        
    ax1.set_xlabel("average pileup")
    ax1.set_ylabel(ylabelLumi)
    ax1.text(0.54, 0.97, "\\bf{CMS}", verticalalignment='top', transform=ax1.transAxes, weight="bold")
    ax1.text(0.65, 0.97, "\\emph{Work in progress}", verticalalignment='top', transform=ax1.transAxes,style='italic')    
    ax1.text(0.54, 0.89, f"Fill {fill}", verticalalignment='top', transform=ax1.transAxes)    

    y = np.array([yy.n for yy in data_fill['zLumiInst_mc'].values])
    yErr = np.array([y.s for y in data_fill['zLumiInst_mc'].values])

    yRef = data_fill['dLRec(/nb)'].values   
    
    ax1.plot(xPU, yRef, #xerr=(xDown, xUp), 
        label="Reference luminosity", color="blue", 
        marker = "_",
        linestyle='', zorder=0)

    ax1.errorbar(xPU, y, #xerr=(xDown, xUp), 
        yerr=yErr, 
        label="Z rate",
        fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
        zorder=1)

    leg = ax1.legend(loc="lower left", ncol=2, frameon=True, framealpha=1.0, fancybox=False, edgecolor="black")
    leg.get_frame().set_linewidth(0.8)

    yMin = min(min(y),min(yRef))
    yMax = max(max(y),max(yRef))

    yRange = yMax - yMin 
    ax1.set_ylim([yMin - yRange*0.5, yMax + yRange*0.2])
    ax1.set_xlim([xMinPU, xMaxPU])
    # ax1.set_xticks(xTicksPU)
    
    if do_ratio:
        ax1.xaxis.set_major_locator(ticker.NullLocator())
    
        ratios = []
        for zlumi, lumi in data_fill[['zLumi_mc','recLumi']].values:
            ratios.append(zlumi/lumi)

        yRatio = np.array([yy.n for yy in ratios])
        yRatioErr = np.array([y.s for y in ratios])
    
        ax2.set_xlabel("average pileup")
        ax2.set_ylabel("Ratio")
        
        ax2.errorbar(xPU, yRatio, xerr=(xDown, xUp), yerr=yRatioErr, 
            label="Z rate measurement",
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize,
            zorder=1)
            
        ax2.plot(np.array([xMinPU, xMaxPU]), np.array([1.0, 1.0]), color="black",linestyle="-", linewidth=1)
        
        ax2.set_ylim([0.951,1.049])
        ax2.set_xlim([xMinPU, xMaxPU])
        # ax2.set_xticks(xTicksPU)

    # align y labels
    ax1.yaxis.set_label_coords(-0.12, 0.5)
    ax2.yaxis.set_label_coords(-0.12, 0.5)

    plt.savefig(outDir+f"/fill{fill}_lumi_pileup.{fmt}")
    plt.close()