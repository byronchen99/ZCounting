import os
os.environ['ZFIT_DISABLE_TF_WARNINGS'] = "1."

import hist
from hist import Hist
import numpy as np
import pdb
import uncertainties as unc

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--label",  default='Work in progress',  type=str, help="specify label ('Work in progress', 'Preliminary', )")
parser.add_argument("--unblind",  default=False, action="store_true", help="Fit on data")
parser.add_argument("-o","--output",  default='Test',  type=str, help="give output dir")
args = parser.parse_args()

outDir = args.output
if not os.path.isdir(outDir):
    os.mkdir(outDir)

# fiducial Z cross section at 13TeV
xsec = 734
xsec_uncertainty = None # no uncertainty on cross section - free floating parameter
# xsec_uncertainty = 0.03   # uncertainty on cross section - gaussian constraint

# ratio of fiducial cross sections at 13.6 and 13TeV
ratio13p6to13 = 1.05


# --- input
# efficiency corrected number of Z bosons for: 2016preVFP, 2016postVFP, 2017, 2017H, 2018,
z_yields = {
    "2016preVFP": 13148175.0,
    "2016postVFP": 11658670.0,
    "2017": 26799689.0,
    "2017H": 144360.0,
    "2018": 41651170.0,
}

# # reference lumi in pb using normtags PHYSICS
# ref_lumi = {
#     "2016preVFP": 18910.381126,
#     "2016postVFP": 16968.409375000003,
#     "2017": 37921.727600000006,
#     "2017H": 199.269742,
#     "2018": 59742.764247,
# }

# reference lumi in pb using normtags PHYSICS (2016), hfet17New_v0 (2017), hfoc18NEW_v0 (2018)
ref_lumi = {
    "2016preVFP": 18910.381126,
    "2016postVFP": 16968.409375000003,
    "2017": 37568.897576999996,
    "2017H": 202.692817,
    "2018": 59629.669634,
    "2022": 38040,   # Just a dummy number for now
}

z_yields["2022"] = ref_lumi["2022"] * xsec * ratio13p6to13

if not args.unblind:
    # Use asymov data
    z_yields = {key: value * xsec for key, value in ref_lumi.items()}

# statistical uncertainties of Z yields
z_yields_stat = {
    "2016preVFP": 0.00030886129248932546,
    "2016postVFP": 0.0003240928650176417,
    "2017": 0.00021328236230086072,
    "2017H": 0.002884547355729706,
    "2018": 0.0001702104409933952,
    "2022": 0.0002,   # Just a dummy number for now    
}

# uncertainties
uncertainties_lumi = {
    "ref_2016": {
        "2016preVFP": 1.04,
        "2016postVFP": 1.04
        },
    "ref_2017": {
        "2017": 1.96,
        "2017H": 1.28
        },
    "ref_2018": {
        "2018": 1.46
    },
    # "ref_2022": {
    #     "2022": 3.0
    # },
    "ref_20172018": {
        "2017": 0.6,
        "2017H": 0.6,
        "2018": 0.2
    },
    "ref_correlated": {
        "2016preVFP": 0.8,
        "2016postVFP": 0.8,
        "2017": 1.04,
        "2017H": 1.04,
        "2018": 2.11
    }
}
uncertainties_z = {
    "z_2016preVFP": {
        "2016preVFP": 0.48
    },
    "z_2016postVFP": {
        "2016postVFP": 0.37
    },
    "z_2017": {
        "2017": 0.3
    },
    "z_2017H": {
        "2017H": 0.3
    },
    "z_2018": {
        "2018": 0.31
    },
    # "z_2022": {
    #     "2022": 2.0
    # },
    "z_2016preVFP2016postVFP20172017H": {
        "2016preVFP": 0.03,
        "2016postVFP": 0.04,
        "2017": 0.03,
        "2017H": 0.14
    },
    "z_20172017H": {
        "2017": 0.08,
        "2017H": 0.08
    },
    "z_2016preVFP2016postVFP20172018": {
        "2016preVFP": 0.39,
        "2016postVFP": 0.06,
        "2017": 0.01,
        "2018": 0.01
    },
    "z_correlated": {
        "2016preVFP": 1.3,
        "2016postVFP": 1.45,
        "2017": 1.37,
        "2017H": 1.53,
        "2018": 1.38
    },
}
 

# uncertainties = {
#     "2016preVFP": 0.1,
#     "2016postVFP": 0.1,
#     "2017": 0.1,
#     "2017H": 0.1,
#     "2018": 0.1
# }

# define eras used in combination
eras = ["2016preVFP", "2016postVFP", "2017", "2017H", "2018"]

# 
nBins = len(eras)
binLo = 0
binHi = nBins

# set the Z counts
hZ = Hist(
    hist.axis.Regular(bins=nBins, start=binLo, stop=binHi, name="x"))

# set yields and variance (statistical uncertainty)
z_yields_uncorrected = {era: (1./(z_yields_stat[era])**2) for era in eras}  # the uncorrected z yields is rederived from the statistical uncertainty
z_weights = {era: z_yields[era]/z_yields_uncorrected[era] for era in eras}  # the per event weights to get from uncorrected to corrected yield
# hZ[:] = [(z_yields[era], z_yields_uncorrected[era] * z_weights[era]**2) for era in eras]      # the variance is given by square of weights times uncorrected yields

# set uncorrected yields for data
hZ[:] = [int(z_yields_uncorrected[era]) for era in eras]

# efficiencies to get from corrected to uncorrected number
z_efficiencies = {era: z_yields_uncorrected[era] / z_yields[era] for era in eras}


# set the reference lumi counts
hists_ref = {}
for i, era in enumerate(eras):
    hist_ref = Hist(
        hist.axis.Regular(bins=nBins, start=binLo, stop=binHi, name="x"))

    # set lumi
    hist_ref[i] = ref_lumi[era]

    hists_ref[era] = hist_ref

print("import zfit")
import zfit
print("imported zfit")


# obs_nobin = zfit.Space('x', (binLo, binHi))

binning = zfit.binned.RegularBinning(nBins, binLo, binHi, name="x")
obs = zfit.Space("x", binning=binning)

data = zfit.data.BinnedData.from_hist(hZ)

# make extended pdfs, each bin is scaled by a separate histogram
# cross section as a common normalization
rate_xsec = zfit.Parameter('r_xsec', 0, -1., 1.)

# nuisance parameters on luminosity
nuisances_lumi = {key: zfit.Parameter('n_{0}'.format(key), 0, -5., 5.) for key in uncertainties_lumi.keys()}

# nuisance parameters on Z yield
nuisances_z = {key: zfit.Parameter('n_{0}'.format(key), 0, -5., 5.) for key in uncertainties_z.keys()}

# all nuisance parameters
nuisances = {**nuisances_lumi, **nuisances_z}

# put gaussian constraints on all nuisances
constraints = {key: zfit.constraint.GaussianConstraint(param, observation=0.0, uncertainty=1.0) for key, param in nuisances.items()}

if xsec_uncertainty != None:
    # put a gaussian constraint on the cross section
    constraints["r_xsec"] = zfit.constraint.GaussianConstraint(rate_xsec, observation=0.0, uncertainty=xsec_uncertainty)

# rearange dictionaries: for each era a list of uncertainties, nuisances
uncertainties_lumi_era = {}
nuisances_lumi_era = {}
uncertainties_z_era = {}
nuisances_z_era = {}
for era in eras:
    uncertainties_lumi_era[era] = []
    nuisances_lumi_era[era] = []
    for key, uncertainty in uncertainties_lumi.items():
        if era in uncertainty.keys():
            uncertainties_lumi_era[era].append(uncertainty[era]/100.+1)
            nuisances_lumi_era[era].append(nuisances_lumi[key])

    uncertainties_z_era[era] = []
    nuisances_z_era[era] = []
    for key, uncertainty in uncertainties_z.items():
        if era in uncertainty.keys():
            uncertainties_z_era[era].append(uncertainty[era]/100.+1)
            nuisances_z_era[era].append(nuisances_z[key])


def get_luminosity_function(era):
    # return function to calculate luminosity for a year

    # central value of luminosity
    central = ref_lumi[era]

    # dictionary with uncertainties for the considered era
    def function(params):
        l = central
        for i, p in enumerate(params):
            l *= uncertainties_lumi_era[era][i]**p
        return l
    
    return function


def get_scale_function(era):

    z_eff = z_efficiencies[era]

    def function(rate_xsec, lumi, params):#, parameters=[]):
        s = xsec * z_eff * (1+rate_xsec) * lumi
        # apply nuisance parameters
        # for p in parameters:
        #     s *= p
        for i, p in enumerate(params):
            s *= uncertainties_z_era[era][i]**p
        return s
    
    return function


models = {}
lumis = {}
scales = {}
p_lumis = {}
for era in eras:
    print("create model for {0}".format(era))

    # p_lumi = zfit.Parameter('l_{0}'.format(era), 1, 0.5, 1.5, floating=True)

    # lumi parameter including uncetainties
    l = zfit.ComposedParameter('lumi_{0}'.format(era),
        get_luminosity_function(era),
        # params=[p_lumi, nuisances_lumi_era[era]]
        params=[nuisances_lumi_era[era]]
        )

    # absolute scale including uncertainties on cross section, acceptance, and efficiencies
    s = zfit.ComposedParameter('scale_{0}'.format(era),
        get_scale_function(era),
        params=[rate_xsec, l, nuisances_z_era[era]] 
        )

    m = zfit.pdf.HistogramPDF(hists_ref[era], extended=s, name="PDF_Bin{0}".format(era))

    lumis[era] = l
    models[era] = m
    scales[era] = s
    # p_lumis[era] = p_lumi


# build composite model
model = zfit.pdf.BinnedSumPDF([m for m in models.values()])

# if not args.unblind:
#     # azimov_hist = model.to_hist()
#     data = model.to_binneddata()

loss = zfit.loss.ExtendedBinnedChi2(model, data, constraints=[c for c in constraints.values()])
# loss = zfit.loss.ExtendedBinnedNLL(model, data, constraints=[c for c in constraints.values()])

# minimization
minimizer = zfit.minimize.Minuit(mode=2)#, gradient=True)
result = minimizer.minimize(loss)

# calculate hesse
result.hesse()

# calculate minuit error by doing likelihood scan
errors, new_result = result.errors(name='errors') # name='minuit_minos'

print("Function minimum:", result.fmin)
print("Converged:", result.converged)
print("Full minimizer information:", result)
# print("Correlation matrix:", result.correlation())

if new_result:
    print("New result was found!")
    print(result)


# error propagation to get uncertainty on luminosity
for era in eras:
    # uncertainties from hesse
    values = unc.correlated_values(
        [result.params[p]["value"] for p in nuisances_lumi_era[era]], 
        result.covariance(nuisances_lumi_era[era]))

    # uncertainties from minos
    err_lo = unc.correlated_values_norm([v for v in zip(
        [result.params[p]["value"] for p in nuisances_lumi_era[era]], 
        [result.params[p]["errors"]["lower"] for p in nuisances_lumi_era[era]])], 
        result.correlation(nuisances_lumi_era[era]))

    err_hi = unc.correlated_values_norm([v for v in zip(
        [result.params[p]["value"] for p in nuisances_lumi_era[era]], 
        [result.params[p]["errors"]["upper"] for p in nuisances_lumi_era[era]])], 
        result.correlation(nuisances_lumi_era[era]))

    lumi = get_luminosity_function(era)(values)
    lumi_lo = get_luminosity_function(era)(err_lo)
    lumi_hi = get_luminosity_function(era)(err_hi)

    print("Lumi[{0}] = {1} (+{2} / -{3})".format(era, lumi, lumi_hi.s, lumi_lo.s))
    print("Del Lumi[{0}] = {1} (+{2} / -{3})".format(era, lumi/lumi.n, lumi_hi.s/lumi.n, lumi_lo.s/lumi.n))


# correlation matrix on luminosities
correlated_values = unc.correlated_values([result.params[p]["value"] for p in nuisances_lumi.values()], result.covariance(nuisances_lumi.values()))

for k,v in zip(nuisances_lumi.values(), correlated_values):
    result.params[k]["correlated_value"] = v

lumi_function = [get_luminosity_function(era) for era in eras]
lumi_values = [lumi_function[i]([result.params[iv]["correlated_value"] for iv in v]) for i, v in enumerate(nuisances_lumi_era.values())]

corr_matrix_lumi = unc.correlation_matrix(lumi_values)

print(corr_matrix_lumi)


all_params = [rate_xsec,] + [n for n in nuisances.values()]

### --- plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

textsize = 15

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino",],
    "font.size": textsize
    # 'text.latex.preamble': [r"""\usepackage{bm}"""]
})

mpl.rcParams.update({
    "legend.fontsize" : "medium",
    "axes.labelsize" : "medium",
    "axes.titlesize" : "medium",
    "xtick.labelsize" : "medium",
    "ytick.labelsize" : "medium",
})

### --- plot correlation matrix
def plot_matrix(matrix, labels=[], correlation=True, name=""):
    fig, ax = plt.subplots()
    im = ax.imshow(matrix)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(labels)):
        for j in range(len(labels)):
            # if correlation and i == j:
            #     # correlation matrix is always 1 on diagonal, so we don't plot it
            #     continue
            val = round(matrix[i, j],3)
            text = ax.text(j, i, val, ha="center", va="center", color="w")

    # ax.set_title("Harvest of local farmers (in tons/year)")
    fig.tight_layout()
    suffix = "_correlation" if correlation else ""
    name = name if name.startswith("_") else "_"+name
    plt.savefig("{0}/matrix{1}{2}.png".format(outDir, name, suffix))


plot_matrix(corr_matrix_lumi, labels = [e for e in eras], name="lumi")

exit()


### --- plot likelihood scan
def plot_scan(param, name="param", profile=True, limits=2.):
    print("Make likelihood scan for {0}".format(name))

    val = result.params[name]["value"]
    err = result.params[name]["hesse"]["error"]

    # scan around +/- 2 sigma (prefit) interval
    xLo = val-limits
    xHi = val+limits

    x = np.linspace(xLo, xHi, num=100)
    y = []
    param.floating = False
    for val in x:
        param.set_value(val)
        if profile:
            minimizer.minimize(loss)
        y.append(loss.value())

    y = (np.array(y) - result.fmin)*2

    param.floating = True
    zfit.param.set_values(loss.get_params(), result)

    # ymin = min(y)

    yargmin = np.argmin(y)
    xmin = x[np.argmin(y)]
    # left and right 68% intervals
    xL1 = x[np.argmin(abs(y[:yargmin]-1))]
    xR1 = x[yargmin+np.argmin(abs(y[yargmin:]-1))]
    # left and right 95% intervals
    xL2 = x[np.argmin(abs(y[:yargmin]-4))]
    xR2 = x[yargmin+np.argmin(abs(y[yargmin:]-4))]

    plt.figure()
    
    if max(y) >= 1:
        plt.plot([xLo, xL1], [1,1], linestyle="dashed", color="gray")
        plt.plot([xR1, xHi], [1,1], linestyle="dashed", color="gray")
        plt.plot([xL1, xL1], [0,1], linestyle="dashed", color="gray")
        plt.plot([xR1, xR1], [0,1], linestyle="dashed", color="gray")

    if max(y) >= 4:
        plt.plot([xLo, xL2], [4,4], linestyle="dashed", color="gray")
        plt.plot([xR2, xHi], [4,4], linestyle="dashed", color="gray")
        plt.plot([xL2, xL2], [0,4], linestyle="dashed", color="gray")
        plt.plot([xR2, xR2], [0,4], linestyle="dashed", color="gray")

    plt.plot(x, y, color="red")
    plt.xlabel(name)
    plt.ylabel("-2 $\\Delta$ ln(L)")

    plt.xlim(xLo,xHi)
    plt.ylim(0,5)

    fig.tight_layout()

    suffix = "_profiling" if profile else "_scan"
    plt.savefig("{0}/loss_{1}{2}.png".format(outDir, name, suffix))

plot_scan(rate_xsec, "r_xsec", limits=0.03, profile=False)
plot_scan(rate_xsec, "r_xsec", limits=0.03)

# for era, p in p_lumis.items():
#     plot_scan(p, "l_"+era, limits=0.1)

for n, p in nuisances.items():
    plot_scan(p, "n_"+n, profile=False)
    plot_scan(p, "n_"+n)


### --- plot pulls scan
