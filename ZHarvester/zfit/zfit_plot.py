import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use("CMS") 
import numpy as np

def plot_model(model, data, scale=1, plot_data=True, label="model", color="red", linestyle="-"):  

    lower, upper = data.data_range.limit1d
    x = np.linspace(lower, upper, num=1000) 
    y = model.ext_pdf(x)
    y *= scale
    plt.plot(x, y, label=label, color=color, linestyle=linestyle)
    
    if plot_data:
        y, x = data.to_hist().to_numpy()
        xerr = (x[1:] - x[:-1])/2.
        x = x[:-1] + (x[1:] - x[:-1])/2.

        plt.errorbar(x, y, xerr=xerr, yerr=np.sqrt(y), label="Data", 
            fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=4)

def plot_comp_model(model, data, name="", bkg=None, sig=None):
    plt.figure()

    if bkg:
        plot_model(bkg, data, plot_data=False, label="bkg", color="green", linestyle="--")
    if sig:
        plot_model(sig, data, plot_data=False, label="sig", color="blue", linestyle="--")

    plot_model(model, data)
    plt.legend()

    plt.savefig("fit_{0}.png".format(name))