from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from scipy.optimize import curve_fit
import uncertainties as unc
import pdb
import numpy as np
import awkward as ak
import vector

# ------------------------------------------------------------------------------
def get_masses(pairs, mass_lo=None, mass_hi=None):
    """
    compute invariant masses
    
    Parameters
    ----------
    pairs : awkward array with pairs
        has to be touples that contain the keys "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"
    mass_lo : float
        lower mass cut
    mass_hi : float
        upper mass cut
    """
    
    # calculate mass of each pair
    l_pt, r_pt = ak.unzip(pairs["Muon_pt"])
    l_eta, r_eta = ak.unzip(pairs["Muon_eta"])
    l_phi, r_phi = ak.unzip(pairs["Muon_phi"])
    l_mass, r_mass = ak.unzip(pairs["Muon_mass"])

    mu1 = vector.obj(pt=l_pt, phi=l_phi, eta=l_eta, mass=l_mass)
    mu2 = vector.obj(pt=r_pt, phi=r_phi, eta=r_eta, mass=r_mass)
    
    masses = (mu1 + mu2).mass
    
    if mass_lo:
        masses = masses[masses > mass_lo]
    if mass_hi:
        masses = masses[masses < mass_hi]
    
    return masses #ak.flatten(masses).to_numpy()

def linear(x, a):
    return a * x 

# # second order polynomial
# def pol2(x, a, b, c):
#     return a * x**2 + b * x + c

# second order polynomial
def pol2(x, a, b):
    return a * x**2 + b * x

def plot_scatter(x, y, xlabel, ylabel, range, zlabel="Density", rangey=None,
                 cutsPtEta='p_\mathrm{T}(\mu) > 30\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4',
                 cutsAdditional=None,
                 title=None,
                 scatter=False,
                 eventWeights=None,
                 saveas='./scat.png'):
    import matplotlib.colors as colors
    plt.clf()
    #cmap = plt.get_cmap("Greys")
    cmap = plt.get_cmap("RdBu")

    if rangey is None:
        rangey = range

    if eventWeights is not None:
        w = eventWeights
    else:
        w = np.ones(len(x))

    if scatter:
        plt.scatter(x, y, weights=w, marker='.', color='k')
    else:
        plt.hist2d(x, y, weights=w, bins=50, range=(range, rangey), cmap=cmap, normed=True, norm=colors.LogNorm())
        plt.colorbar(label=zlabel)
    corr, _ = pearsonr(x, y)
    
    print("make fit")    
    func = linear
    popt, pcov = curve_fit(func, x, y)
    perr = np.sqrt(np.diag(pcov))
    # correlated values
    params = unc.correlated_values(popt, pcov)
    f = lambda ix: func(ix, *params)
    print("Fit params: {0}".format(params))
        
    xx = np.arange(range[0], range[1], 0.1)
    plt.plot(xx, np.array([f(ix).n for ix in xx]), "k--")

    # plot line with slope 1 for comparison
    plt.plot(range, range,"k-")
            
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    xtext= (range[1] - range[0])*0.05 + range[0]
    ytext = (range[1] - range[0])
    box = dict(boxstyle='round', facecolor='white', alpha=0.9)
    textstr = '\n'.join((
        # r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \mu\mu$',
        r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \ell\ell$',
        # r'$56\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
        # r'${0}$'.format(cutsPtEta),
        # r'$\Delta R(\mu, \mu) > 0.4$',
        # r'$\rho_\mathrm{pearson} = $' + '${0}$'.format(round(corr, 3))
    ))
    if cutsAdditional:
        textstr = '\n'.join((textstr, r'${0}$'.format(cutsAdditional)))
    plt.text(xtext, range[0] + 0.95*ytext, textstr, verticalalignment='top', bbox=box)
    if title:
        plt.title(title, loc='right')
    if range:
        plt.xlim(range)
        plt.ylim(rangey)
    print("Plot {0}".format(saveas))
    plt.savefig(saveas)
    plt.close()