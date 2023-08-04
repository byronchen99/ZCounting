import zfit
import numpy as np
import uncertainties as unc
from hist import Hist
from zfit_tools.plot import plot_model
from zfit_tools.models import get_background

import pdb

def fit(templates_1d, data_hist, category):
    # obtain normalized 1d arrays
    correction_types = ["nominal", "plus_scale", "minus_scale", "plus_spread", "minus_spread"]

    normalized_array = {}
    for correction in correction_types:
        normalized_array[correction] = templates_1d["{0}_{1}".format(category, correction)] / templates_1d["{0}_nominal".format(category)]
        normalized_array[correction] = np.nan_to_num(normalized_array[correction], nan=1)


    # create binned templates for zfit
    if category == "Glopass" or category == "Glofail":
        nbins, bin_min, bin_max = 70, 56, 126
    else:
        nbins, bin_min, bin_max = 50, 66, 116
    
    hist = {}
    for template in normalized_array:
        hist[template] = Hist.new.Regular(nbins, bin_min, bin_max, name="x").Double()
        for i in range(len(normalized_array[template])):
            hist[template][i] = normalized_array[template][i]

    hist["original"] = Hist.new.Regular(nbins, bin_min, bin_max, name="x").Double()
    for i in range(len(templates_1d["{0}_nominal".format(category)])):
        hist["original"][i] = templates_1d["{0}_nominal".format(category)][i]

    binned_template = {}
    for template in hist:
        binned_template[template] = zfit.data.BinnedData.from_hist(hist[template])


    
    templates_scale = (binned_template["minus_scale"], binned_template["nominal"], binned_template["plus_scale"])
    templates_spread = (binned_template["minus_spread"], binned_template["nominal"], binned_template["plus_spread"])

    # define parameters
    alpha = zfit.Parameter('alpha_{0}'.format(category), 0, -10, 10, step_size=0.01)
    beta = zfit.Parameter('beta_{0}'.format(category), 0, -10, 10, step_size=0.01)
    sig_yield = zfit.Parameter('sig_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)

    # define PDFs
    morphing_pdf_scale = zfit.pdf.SplineMorphingPDF(alpha, templates_scale)
    morphing_pdf_spread = zfit.pdf.SplineMorphingPDF(beta, templates_spread)
    nominal = zfit.pdf.HistogramPDF(binned_template["original"])

    pdfs = (morphing_pdf_scale, morphing_pdf_spread, nominal)
    productPDF = zfit.pdf.ProductPDF(pdfs)
    productPDF.set_yield(sig_yield)

    binning = zfit.binned.RegularBinning(nbins, bin_min, bin_max, name="x")
    obs_binned = zfit.Space("x", binning=binning)
    productPDF_binned = zfit.pdf.BinnedFromUnbinnedPDF(productPDF, obs_binned)

    # impose constraints
    constraint_alpha = zfit.constraint.GaussianConstraint(alpha, 0, 1)
    constraint_beta = zfit.constraint.GaussianConstraint(beta, 0, 1)
    constraints = (constraint_alpha, constraint_beta)


    # combinatorial background
    obs_nobin = zfit.Space('x', (bin_min, bin_max))
    bkg1 = get_background(obs_nobin, "chebyshev", category=category)
    bkg2 = get_background(obs_nobin, "exp", category=category)

    bkg1_yield = zfit.Parameter('bkg1_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)
    bkg2_yield = zfit.Parameter('bkg2_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)
    bkg1.set_yield(bkg1_yield)
    bkg2.set_yield(bkg2_yield)

    bkg = zfit.pdf.SumPDF((bkg1, bkg2))
    bkg_binned = zfit.pdf.BinnedFromUnbinnedPDF(bkg, obs_binned)

    combinedPDF = zfit.pdf.BinnedSumPDF((productPDF_binned, bkg_binned))
   

    # zfit binning of data
    binned_data = zfit.data.BinnedData.from_hist(data_hist)

    # fitting
    # loss = zfit.loss.ExtendedBinnedNLL(productPDF_binned, binned_data, constraints = constraints)
    loss = zfit.loss.ExtendedBinnedNLL(combinedPDF, binned_data, constraints = constraints)
    minimizer = zfit.minimize.Minuit()
    result = minimizer.minimize(loss)

    
    status = result.valid
    print(f"status: {status}")

    # best fit params
    sig_yield_nominal = result.params['sig_yield_{0}'.format(category)]["value"]
    alpha_nominal = result.params['alpha_{0}'.format(category)]["value"]
    beta_nominal = result.params['beta_{0}'.format(category)]["value"]
    bkg1_yield_nominal = result.params['bkg1_yield_{0}'.format(category)]["value"]
    bkg2_yield_nominal = result.params['bkg2_yield_{0}'.format(category)]["value"]

    # get uncertainties
    # try:
    #     hessval = result.loss.hessian(list(result.params)).numpy()
    #     cov = np.linalg.inv(hessval)
    #     eigvals = np.linalg.eigvalsh(hessval)
    #     covstatus = eigvals[0] > 0.
    # except:
    #     cov = None
    #     covstatus = False

    # print(f"covariance status: {covstatus}")

    # get uncertainties
    param_hesse = result.hesse()
    sig_yield_unc = param_hesse[[p for p in result.params if p.name ==  "sig_yield_{0}".format(category)][0]]["error"]
    alpha_unc = param_hesse[[p for p in result.params if p.name ==  "alpha_{0}".format(category)][0]]["error"]
    beta_unc = param_hesse[[p for p in result.params if p.name ==  "beta_{0}".format(category)][0]]["error"]
    bkg1_yield_unc = param_hesse[[p for p in result.params if p.name ==  "bkg1_yield_{0}".format(category)][0]]["error"]
    bkg2_yield_unc = param_hesse[[p for p in result.params if p.name ==  "bkg2_yield_{0}".format(category)][0]]["error"]

    
    print(result)
    
    sig_yield = unc.ufloat(sig_yield_nominal, sig_yield_unc)
    alpha = unc.ufloat(alpha_nominal, alpha_unc) 
    beta = unc.ufloat(beta_nominal, beta_unc)
    bkg_yield = unc.ufloat(bkg1_yield_nominal, bkg1_yield_unc) + unc.ufloat(bkg2_yield_nominal, bkg2_yield_unc)


    # plot_comp_model(model_nobin, data, bkg=bkg, name=category)

    # plotting
    model_total = combinedPDF.values()
    model_bkg = bkg_binned.values()
    data = data_hist.view()
    plot_model(model_total, model_bkg, data, nbins, bin_min, bin_max, name=category)


    return sig_yield, alpha, beta, bkg_yield
    # return result.params['sig_yield_{0}'.format(category)]["value"], result.params['alpha_{0}'.format(category)]["value"], result.params['beta_{0}'.format(category)]["value"]


def fit_2(templates_1d, data_hist, category): # temporary solution for templates without corrections yet, only gives yield
    template = templates_1d["{0}_none".format(category)]

    if category == "Glopass" or category == "Glofail":
        nbins, bin_min, bin_max = 70, 56, 126
    else:
        nbins, bin_min, bin_max = 50, 66, 116
    
    hist = Hist.new.Regular(nbins, bin_min, bin_max, name="x").Double()
    for i in range(len(template)):
        hist[i] = template[i]
    
    binned_template = zfit.data.BinnedData.from_hist(hist)

    sig_yield = zfit.Parameter('sig_yield_{0}'.format(category), sum(data_hist), 0, sum(data_hist)*1.5, step_size=1)
    func = zfit.pdf.HistogramPDF(binned_template, sig_yield)

    # combinatorial background
    obs_nobin = zfit.Space('x', (bin_min, bin_max))
    bkg1 = get_background(obs_nobin, "chebyshev", category=category)
    bkg2 = get_background(obs_nobin, "exp", category=category)

    bkg1_yield = zfit.Parameter('bkg1_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)
    bkg2_yield = zfit.Parameter('bkg2_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)
    bkg1.set_yield(bkg1_yield)
    bkg2.set_yield(bkg2_yield)

    bkg = zfit.pdf.SumPDF((bkg1, bkg2))

    binning = zfit.binned.RegularBinning(nbins, bin_min, bin_max, name="x")
    obs_binned = zfit.Space("x", binning=binning)
    bkg_binned = zfit.pdf.BinnedFromUnbinnedPDF(bkg, obs_binned)

    combinedPDF = zfit.pdf.BinnedSumPDF((func, bkg_binned))

    # zfit binning of data
    binned_data = zfit.data.BinnedData.from_hist(data_hist)
    
    # loss = zfit.loss.ExtendedBinnedNLL(sig, binned_data, options = { "numhess" : False })
    loss = zfit.loss.ExtendedBinnedNLL(combinedPDF, data_hist)
    # minimizer = zfit.minimize.ScipyTrustConstrV1(hessian = "zfit")
    minimizer = zfit.minimize.Minuit()
    result = minimizer.minimize(loss)

    
    status = result.valid
    print(f"status: {status}")

    
    # best fit params
    sig_yield_nominal = result.params['sig_yield_{0}'.format(category)]["value"]
    bkg1_yield_nominal = result.params['bkg1_yield_{0}'.format(category)]["value"]
    bkg2_yield_nominal = result.params['bkg2_yield_{0}'.format(category)]["value"]

    # get uncertainties
    param_hesse = result.hesse()
    sig_yield_unc = param_hesse[[p for p in result.params if p.name ==  "sig_yield_{0}".format(category)][0]]["error"]
    bkg1_yield_unc = param_hesse[[p for p in result.params if p.name ==  "bkg1_yield_{0}".format(category)][0]]["error"]
    bkg2_yield_unc = param_hesse[[p for p in result.params if p.name ==  "bkg2_yield_{0}".format(category)][0]]["error"]

    print(result)

    sig_yield = unc.ufloat(sig_yield_nominal, sig_yield_unc)
    bkg_yield = unc.ufloat(bkg1_yield_nominal, bkg1_yield_unc) + unc.ufloat(bkg2_yield_nominal, bkg2_yield_unc)
    

    # plotting
    model_total = combinedPDF.values()
    model_bkg = bkg_binned.values()
    data = data_hist.view()
    plot_model(model_total, model_bkg, data, nbins, bin_min, bin_max, name=category)
    
    return sig_yield, bkg_yield