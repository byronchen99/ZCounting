from pandas import DataFrame
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

#write minitree into dataframe
def tree_to_df(tree, arrSize=5):
    #if tree has arrays with variable length, split into 'arrSize' new columns and fill empty values with 'NaN'
    df = DataFrame()
    for key in tree.dtype.names:
        if len(tree[key].shape) == 2:
            for i in range(tree[key].shape[1]):
                df[key+"_"+str(i)] = tree[key][:,i]
        else:
            df[key] = tree[key]
            # zero padding because variable length arrays in dataframe columns can not be stored properly
            for key in df.keys():
                if df[key].dtype =='O':
                    for i in range(arrSize):
                        df[key+"_"+str(i)] = df[key].apply(lambda x: x[i] if len(x) > i else float('NaN'))
                    df = df.drop([key],axis=1)
    return df


def plot_scatter(x, y, xlabel, ylabel, range,
                 cutsPtEta='p_\mathrm{t}(\mu) > 27\ \mathrm{GeV} \qquad |\eta(\mu)| < 2.4',
                 cutsAdditional=None,
                 title=None,
                 saveas='./scat.png'):
    plt.clf()
    plt.scatter(x, y, marker='.', color='k')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    corr, _ = pearsonr(x, y)
    xtext= (range[1] - range[0])*0.05 + range[0]
    ytext = (range[1] - range[0])
    box = dict(boxstyle='round', facecolor='white', alpha=0.9)
    textstr = '\n'.join((
        r'$66\ \mathrm{GeV} < \mathrm{M}_{\mu\mu} < 116\ \mathrm{GeV}$',
        r'${0}$'.format(cutsPtEta),
        r'$\Delta R(\mu, \mu) > 0.4$',
        r'$\rho_{pearson} = $' + '${0}$'.format(round(corr, 3))
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
