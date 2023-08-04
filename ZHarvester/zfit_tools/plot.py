import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS") 
import numpy as np

def plot_model(model_total, model_bkg, data, nbins, bin_min, bin_max, name, plot_data=True):  

    x = np.linspace(bin_min, bin_max, nbins+1) 
    y_total = model_total
    y_bkg = model_bkg

    fig, ax = plt.subplots()
    plt.stairs(y_total, x, label="sig + bkg", color="red", linestyle="-")
    plt.stairs(y_bkg, x, label="bkg", color="blue", linestyle="-")
    
    if plot_data:
        x = np.linspace(bin_min+0.5, bin_max-0.5, nbins)

        y = data
        xerr = np.ones(nbins) * 0.5


        plt.errorbar(x, y, xerr=xerr, yerr=np.sqrt(y), label="Data", 
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=4)
    
    ax.set_xlabel("Z boson candidate mass [GeV]", loc="center")
    ax.set_ylabel("Z boson candidates / 1 GeV", loc="center")
    plt.xlim([bin_min, bin_max])
    plt.legend()
    plt.savefig("fit_{0}.png".format(name), bbox_inches='tight')

    plt.close()


# old code below:

# def plot_model(model, data, scale=1, plot_data=True, label="model", color="red", linestyle="-"):  

#     lower, upper = data.data_range.limit1d
#     x = np.linspace(lower, upper, num=1000) 
#     y = model.ext_pdf(x)
#     y *= scale
#     plt.plot(x, y, label=label, color=color, linestyle=linestyle)
    
#     if plot_data:
#         y, x = data.to_hist().to_numpy()
#         xerr = (x[1:] - x[:-1])/2.
#         x = x[:-1] + (x[1:] - x[:-1])/2.

#         plt.errorbar(x, y, xerr=xerr, yerr=np.sqrt(y), label="Data", 
#             fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=4)


# def plot_comp_model(model, data, name="", bkg=None, sig=None):
#     plt.figure()

#     if bkg:
#         plot_model2(bkg, data, plot_data=False, label="bkg", color="green", linestyle="--")
#     if sig:
#         plot_model2(sig, data, plot_data=False, label="sig", color="blue", linestyle="--")

#     plot_model2(model, data)
#     plt.legend()

#     plt.savefig("fit_{0}.png".format(name))