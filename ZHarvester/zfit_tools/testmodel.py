import uproot
import matplotlib.pyplot as plt
import vector
from hist import Hist
import numpy as np
# from templates import generate_template,unfolding

import zfit

# import zfit_physics as zphysics


import pdb


def return_mc_data_reco(filename_ntuples):

    f1 = uproot.open(filename_ntuples)
    treename = "zcounting/tree"
    data = uproot.open(filename_ntuples+":"+treename)

    data_reco = data.arrays(["muon_genPt","muon_genEta", "muon_genPhi", "antiMuon_genPt","antiMuon_genEta","antiMuon_genPhi", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_ID", "Muon_triggerBits", "nPV"], "(decayMode==13) & (nMuon>=2)")

    return data_reco


# filename_ntuples = "C:/Users/byron/ZCounting/DYJetsToLL_M_50_LO_FlatPU0to75_Autumn18.root"
# data_reco = return_mc_data_reco(filename_ntuples)


# mc_data_PV = generate_template(data_reco, "HLT2")

# test PV data
# PV_data = [3.40000e+01, 4.50000e+01, 2.28000e+02, 1.00900e+03, 3.57900e+03,
#        9.42500e+03, 2.13230e+04, 4.27480e+04, 7.67070e+04, 1.23301e+05,
#        1.82771e+05, 2.51702e+05, 3.23823e+05, 3.95676e+05, 4.60453e+05,
#        5.12120e+05, 5.46302e+05, 5.66748e+05, 5.70929e+05, 5.59420e+05,
#        5.35546e+05, 5.02174e+05, 4.64715e+05, 4.19288e+05, 3.74842e+05,
#        3.29376e+05, 2.84960e+05, 2.42599e+05, 2.03105e+05, 1.68695e+05,
#        1.36796e+05, 1.09307e+05, 8.63000e+04, 6.73620e+04, 5.13320e+04,
#        3.86330e+04, 2.87350e+04, 2.12400e+04, 1.51040e+04, 1.06400e+04,
#        7.54900e+03, 5.15700e+03, 3.50500e+03, 2.32400e+03, 1.50700e+03,
#        9.84000e+02, 7.28000e+02, 4.19000e+02, 2.54000e+02, 1.76000e+02,
#        1.03000e+02, 5.90000e+01, 4.30000e+01, 2.70000e+01, 1.40000e+01,
#        6.00000e+00, 2.00000e+00, 3.00000e+00, 0.00000e+00, 1.00000e+00,
#        1.00000e+00, 0.00000e+00, 0.00000e+00, 1.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,
#        0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00]

# PV_data = PV_data / np.sum(PV_data)

# hPV = Hist.new.Regular(100, 0.5, 100.5, name = 'PV').Double()
# for i in range(len(PV_data)):
#     hPV[i] = PV_data[i]



# mc_data = unfolding(mc_data_PV, hPV)

# test data set
# data = [  234.,   259.,   235.,   330.,   339.,   353.,   384.,   432.,
#          493.,   542.,   619.,   691.,   786.,   891.,  1148.,  1375.,
#         1837.,  2377.,  3031.,  4240.,  6148.,  9462., 15088., 23601.,
#        31513., 30315., 21900., 13332.,  7677.,  4784.,  3087.,  2036.,
#         1605.,  1232.,   996.,   845.,   671.,   619.,   490.,   444.,
#          389.,   345.,   310.,   299.,   248.,   236.,   221.,   207.,
#          194.,   186.]

# data_hist = Hist.new.Regular(50, 66, 116, name="Mass").Double()

# for i in range(len(data)):
#     data_hist[i] = data[i]

# binned_data = zfit.data.BinnedData.from_hist(data_hist)

# pdb.set_trace()


# troubleshoot with constant histogram
# constant = np.full(50,1)
# mc_data_hist = Hist.new.Regular(50, 66, 116, name="x").Double()
# for i in range(len(constant)):
#     mc_data_hist[i] = constant[i]*3
# mc_data = zfit.data.BinnedData.from_hist(mc_data_hist)

# data_hist = Hist.new.Regular(50, 66, 116, name="x").Double()
# for i in range(len(constant)):
#     data_hist[i] = constant[i]
# binned_data = zfit.data.BinnedData.from_hist(data_hist)




# # just a histogram fit
def template_fit(data_hist, template, category):
    scale = zfit.Parameter('sig_yield_{0}'.format(category), sum(data_hist), 0, sum(data_hist)*1.5, step_size=1)
    func = zfit.pdf.HistogramPDF(template, scale)
    
    # loss = zfit.loss.ExtendedBinnedNLL(sig, binned_data, options = { "numhess" : False })
    loss = zfit.loss.ExtendedBinnedNLL(func, data_hist)
    # minimizer = zfit.minimize.ScipyTrustConstrV1(hessian = "zfit")
    minimizer = zfit.minimize.Minuit()
    result = minimizer.minimize(loss)
    
    return result.params['sig_yield_{0}'.format(category)]["value"]

# pdb.set_trace()


# binned Gaussian fit
# sig = zfit.pdf.Gauss(obs=zfit.Space('x', (-10., 10.)), mu=zfit.Parameter("mu", 90, 88, 92), sigma=zfit.Parameter("sigma", 5, 2, 8))
# sig_yield = zfit.Parameter("yield", sum(data_hist), 0, sum(data_hist)*1.5, step_size=1)
# sig.set_yield(sig_yield)

# binning = zfit.binned.RegularBinning(50, 66, 116, name="x")
# obs = zfit.Space("x", binning=binning)
# binned_sig = zfit.pdf.BinnedFromUnbinnedPDF(sig, obs)


# loss = zfit.loss.ExtendedBinnedNLL(binned_sig, binned_data)

# minimizer = zfit.minimize.Minuit()



# convolution fit (doesn't work_)
# scale = zfit.Parameter("scale", sum(data_hist), 0, sum(data_hist)*1.5, step_size=1)
# func = zfit.pdf.HistogramPDF(mc_data)

# resolution = zfit.pdf.Gauss(obs=zfit.Space('x', (-10., 10.)), mu=zfit.Parameter("mu", 0, -2.5, 2.5), sigma=zfit.Parameter("sigma", 2, 0.1, 5))

# binning = zfit.binned.RegularBinning(50, 66, 116, name="x")
# obs = zfit.Space("x", binning=binning)
# binned_resolution = zfit.pdf.BinnedFromUnbinnedPDF(resolution, obs)

# sig = zfit.pdf.FFTConvPDFV1(func, binned_resolution, interpolation="spline:3")

# sig_yield = zfit.Parameter("yield", sum(data_hist), 0, sum(data_hist)*1.5, step_size=1)
# sig.set_yield(sig_yield)

# binned_sig = zfit.pdf.BinnedFromUnbinnedPDF(sig, obs)


# loss = zfit.loss.ExtendedBinnedNLL(binned_sig, binned_data)



