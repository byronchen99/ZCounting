from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def plot_scatter(x, y, xlabel, ylabel, range,
                 cutsPtEta='p_\mathrm{t}(\mu) > 30\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4',
                 cutsAdditional=None,
                 title=None,
                 saveas='./scat.png'):
    import matplotlib.colors as colors
    plt.clf()
    #cmap = plt.get_cmap("Greys")
    cmap = plt.get_cmap("RdBu")

    plt.scatter(x, y, marker='.', color='k')
    #plt.hist2d(x, y, bins=50, range=(range, range), cmap=cmap, normed=True, norm=colors.LogNorm())
    #plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    corr, _ = pearsonr(x, y)
    xtext= (range[1] - range[0])*0.05 + range[0]
    ytext = (range[1] - range[0])
    box = dict(boxstyle='round', facecolor='white', alpha=0.9)
    textstr = '\n'.join((
        r'$\mathrm{\mathbf{CMS}}$ Simulation $Z \rightarrow\ \mu\mu$',
        r'$56\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
        r'${0}$'.format(cutsPtEta),
    #    r'$\Delta R(\mu, \mu) > 0.4$',
    #    r'$\rho_{pearson} = $' + '${0}$'.format(round(corr, 3))
    ))
    if cutsAdditional:
        textstr = '\n'.join((textstr, r'${0}$'.format(cutsAdditional)))
    plt.text(xtext, range[0] + 0.95*ytext, textstr, verticalalignment='top', bbox=box)
    if title:
        plt.title(title, loc='right')
    if range:
        plt.xlim(range)
        plt.ylim(range)
    plt.savefig(saveas)
    plt.close()
