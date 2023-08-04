import zfit
import numpy as np
import uncertainties as unc
from math import sqrt
from zfit_tools.models import get_signal, get_background
from zfit_tools.plot import plot_comp_model

import pdb

# expected signal fraction for initial
signal_fraction = {
    "HLT_2": 0.99,
    "HLT_1": 0.98,
    "ID_pass": 0.98,
    "ID_fail": 0.5,
    "Glo_pass": 0.95,
    "Glo_fail": 0.5,
    "Sta_pass": 0.95,
    "Sta_fail": 0.5
}


# def fit(hist, nBins, binLo, binHi, category, implementation="roofit", **kwargs):

#     if implementation=="zfit":
#         zfit(hist, nBins, binLo, binHi, category, kwargs)
#     elif implementation=="roofit":
#         roofit(hist, nBins, binLo, binHi, category, kwargs)

# def roofit(hist, nBins, binLo, binHi, category):

#     h2HLT = utils.np_to_hist(h2HLT, nBins, binLo, binHi, "TH1", category)

#     ROOT.getZyield(h2HLT, m, "HLT", etaRegion, sigModel, bkgModelPass, 2, sigTemplates, 0)


# new fit (NOT DONE)
def fit(data, template_1d, category):
    alpha_scale = zfit.Parameter('alpha_scale', 0, -10, 10, step_size=0.01)
    alpha_spread = zfit.Parameter('alpha_spread', 0, -10, 10, step_size=0.01)
    sig_yield_scale = zfit.Parameter("sig_yield_scale", sum(data), 0, sum(data)*1.5, step_size=1)
    sig_yield_spread = zfit.Parameter("sig_yield_spread", sum(data), 0, sum(data)*1.5, step_size=1)
    morphing_pdf_scale = zfit.pdf.SplineMorphingPDF(alpha_scale, templates1, sig_yield_scale)
    morphing_pdf_spread = zfit.pdf.SplineMorphingPDF(alpha_spread, templates2, sig_yield_spread)


# old fit
# def fit(hist, nBins, binLo, binHi, category):


    # category = category.replace(" ","_")

    # # --- unbinned
    # # create observable space
    # obs_nobin = zfit.Space('x', (binLo, binHi))

    # # --- binned
    # binning = zfit.binned.RegularBinning(nBins, binLo, binHi, name="x")
    # obs = zfit.Space("x", binning=binning)

    # data = zfit.data.BinnedData.from_hist(hist)

    # # # -- create model

    # # signal component
    # sig = get_signal(obs_nobin, "bw", "gauss", category=category)


    # sig_yield = zfit.Parameter('sig_yield_{0}'.format(category), sum(hist)*signal_fraction[category], 0, sum(hist)*1.5, step_size=1)
    # sig.set_yield(sig_yield)

    # # combinatorial background
    # bkg = get_background(obs_nobin, "chebyshev", category=category)

    # bkg_yield = zfit.Parameter('bkg_yield_{0}'.format(category), sum(hist)*(1-signal_fraction[category]), 0, sum(hist)*1.5, step_size=1)
    # bkg.set_yield(bkg_yield)

    # model_nobin = zfit.pdf.SumPDF([sig, bkg])

    # # plot_comp_model(model_nobin, data, name="prefit", bkg=bkg)x

    # model = zfit.pdf.BinnedFromUnbinnedPDF(model_nobin, obs)
    # loss = zfit.loss.ExtendedBinnedNLL(model, data, options = { "numhess" : False })


    # # minimization
    # # minimizer = zfit.minimize.Minuit()
    # minimizer = zfit.minimize.ScipyTrustConstrV1(hessian = "zfit")
    # result = minimizer.minimize(loss)

    # status = result.valid

    # print(f"status: {status}")

    # try:
    #     hessval = result.loss.hessian(list(result.params)).numpy()
    #     cov = np.linalg.inv(hessval)
    #     eigvals = np.linalg.eigvalsh(hessval)
    #     covstatus = eigvals[0] > 0.
    # except:
    #     cov = None
    #     covstatus = False

    # print(f"covariance status: {covstatus}")

    # print(result)

    # values = [result.params[p]["value"] for p in result.params]

    # pdb.set_trace()

    # plot_comp_model(model_nobin, data, bkg=bkg, name=category)

    # # return parameter values in a dictionary
    # # parameter_value = {key:result.params[key] for key in ["sig_yield_{0}".format(category), "bkg_yield_{0}".format(category)]}
    # # parameter_value = {key.name: result.params[key]["value"] for key in result.params}


    # chi2 = 1

    # Nsig = unc.ufloat(result.params[f"sig_yield_{category}"]["value"],sqrt(cov[0][0]) ) 


    # return Nsig, result.params[f"bkg_yield_{category}"]["value"], chi2, result.params[f"sig_mu_{category}"]["value"], result.params[f"sig_sigma_{category}"]["value"]