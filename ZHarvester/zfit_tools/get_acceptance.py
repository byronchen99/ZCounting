import zfit
import numpy as np
import uncertainties as unc
from hist import Hist
import h5py
import matplotlib.pyplot as plt


def get_acceptance(morphing_params, acceptance_templates, indexing = "", plot = False):
    # obtain normalized 1d arrays
    normalized_array = {}
    for template in acceptance_templates:
        normalized_array[template] = acceptance_templates[template] / acceptance_templates["nominal"]
        normalized_array[template] = np.nan_to_num(normalized_array[template], nan=1)
    
    acceptance_templates["nominal"] = np.nan_to_num(acceptance_templates["nominal"], nan=1)
    

    # create binned templates for zfit
    nbins, bin_min, bin_max = 100, 0.5, 100.5
    hist = {}
    for template in normalized_array:
        hist[template] = Hist.new.Regular(nbins, bin_min, bin_max, name="x").Double()
        for i in range(len(normalized_array[template])):
            hist[template][i] = normalized_array[template][i]

    hist["original"] = Hist.new.Regular(nbins, bin_min, bin_max, name="x").Double()
    for i in range(len(acceptance_templates["nominal"])):
        hist["original"][i] = acceptance_templates["nominal"][i]


    binned_template = {}
    for template in hist:
        binned_template[template] = zfit.data.BinnedData.from_hist(hist[template])



    values_alpha, values_beta = {}, {}
    for template_type in ["HLT2", "HLT1", "IDfail", "Stapass", "Stafail"]:
        templates_scale = (binned_template["{0}_minusalpha".format(template_type)], binned_template["nominal"], binned_template["{0}_plusalpha".format(template_type)])
        templates_spread = (binned_template["{0}_minusbeta".format(template_type)], binned_template["nominal"], binned_template["{0}_plusbeta".format(template_type)])

        # define parameters
        alpha = zfit.Parameter('alphafit_{0}{1}'.format(template_type, indexing), 0, -10, 10, step_size=0.01)
        beta = zfit.Parameter('betafit_{0}{1}'.format(template_type, indexing), 0, -10, 10, step_size=0.01)
        # sig_yield = zfit.Parameter('sig_yield_{0}'.format(category), sum(data_hist.view()), 0, sum(data_hist.view())*1.5, step_size=1)

        # define PDFs
        morphing_pdf_scale = zfit.pdf.SplineMorphingPDF(alpha, templates_scale)
        morphing_pdf_spread = zfit.pdf.SplineMorphingPDF(beta, templates_spread)

        alpha.set_value(morphing_params["{0}_alpha".format(template_type)])
        beta.set_value(morphing_params["{0}_beta".format(template_type)])

        values_alpha[template_type] = morphing_pdf_scale.values()
        values_beta[template_type] = morphing_pdf_spread.values()
        
    nominal = acceptance_templates["nominal"]
    # pdb.set_trace()

    values = nominal * values_alpha["HLT2"] * values_alpha["HLT1"] * values_alpha["IDfail"] * values_alpha["Stapass"] * values_alpha["Stafail"] \
     * values_beta["HLT2"] * values_beta["HLT1"] * values_beta["IDfail"] * values_beta["Stapass"] * values_beta["Stafail"]

    acceptance = values.numpy()

    # plot
    if plot == True:
        fig, ax = plt.subplots()
        bins = np.linspace(0.5, 100.5, num=101)
        plt.stairs(acceptance, bins, label="data", color = "red", linestyle="-")
        plt.stairs(acceptance_templates["nominal"], bins, label="reco_nominal", color = "blue", linestyle="-")
        plt.plot([0.5,100.5], [1,1], color = "black")
        ax.set_xlabel("PV", loc="center")
        ax.set_ylabel("Acceptance", loc="center")
        plt.xlim([0.5, 60.5])
        plt.ylim([0.99, 1.01])
        plt.legend()
        plt.savefig("acceptance.png", bbox_inches='tight')

        plt.close()

    return acceptance

def acceptance_uncertainties(morphing_params_nominal, morphing_params_unc, acceptance_templates, acceptance_nominal):
    # load morphing_params:
    morphing_params = {}
    for key in morphing_params_nominal.keys():
        morphing_params[key] = morphing_params_nominal[key]

    # vary one morphing param at a time upwards:
    uncertainty_plus = np.zeros(100)
    i=0 # just a running index to ensure zfit parameters get different names
    for key in morphing_params_nominal.keys():
        morphing_params[key] = morphing_params_nominal[key] + morphing_params_unc[key]
        acceptance = get_acceptance(morphing_params, acceptance_templates, i)
        uncertainty_term = (acceptance - acceptance_nominal) **2
        uncertainty_plus = uncertainty_plus + uncertainty_term

        i = i+1
        morphing_params[key] = morphing_params_nominal[key]
    
    
    uncertainty_minus = np.zeros(100)
    for key in morphing_params_nominal.keys():
        morphing_params[key] = morphing_params_nominal[key] - morphing_params_unc[key]
        acceptance = get_acceptance(morphing_params, acceptance_templates, i)
        uncertainty_term = (acceptance - acceptance_nominal) **2
        uncertainty_minus = uncertainty_minus + uncertainty_term
        
        i = i+1
        morphing_params[key] = morphing_params_nominal[key]
    
    return np.sqrt(uncertainty_plus), np.sqrt(uncertainty_minus)
    
