import pickle
import argparse
import ROOT
import pandas as pd
import pdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(prog='./Efficiencies')
parser.add_argument(
    '-i', '--input', type=str, default=None,
    help='specify input result file'
)
parser.add_argument(
    '-o', '--output', default="./",
    help='specify output'
)
args = parser.parse_args()

results = pickle.load(open(args.input,"r"))


def plot(_x, _y, _y_err, _name, _key, _xaxis, _yaxis):

    def linear(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(linear, _x, _y, sigma=_y_err)
    perr = np.sqrt(np.diag(pcov))

    plt.clf()
    fig, ax = plt.subplots()

    ax.plot(_x, linear(_x, *popt), 'm-', linewidth=2,  label="fit")
    ax.fill_between(_x, linear(_x, *(popt-perr)), linear(_x, *(popt+perr)), color='magenta', alpha=0.3)
    ax.text(0.05, 0.95, r'fit: $(%5.3f\ \pm\ %5.3f) \cdot x + %5.3f\ \pm\ %5.3f $' % (popt[0],perr[0], popt[1],perr[1]),
        color='magenta', transform=ax.transAxes)
    ax.text(0.05, 0.9, r'mean: $ %5.3f\ \pm\ %5.3f $' % (_y.mean(), _y.std()), color='magenta', transform=ax.transAxes)

    ax.errorbar(_x, _y, xerr=np.zeros(len(_x)), yerr=_y_err, fmt='bo', label=_key)

    xmin = min(_x) - 0.01 *(max(_x) - min(_x))
    xmax = max(_x) + 0.01 *(max(_x) - min(_x))

    ax.set(xlim=(xmin, xmax))

    ax.legend()

    ax.set_ylabel(_xaxis)
    ax.set_xlabel(_yaxis)

    plt.savefig(args.output+'/{0}_{1}.png'.format(_name,_key))
    plt.close()

df = pd.concat([df if key != "H" else None for key, df in results.iteritems()])

for key, val in results.iteritems():
    if key == "H":
        print(val)
        continue
    print("plot "+key)
    plot(val['avgPV'], val['fr']*100, val['fr_err']*100, "fakerate", key, 'fakerate in %', "avgPV")
    plot(val['avgPV'], val['tf'], val['tf_err'], "transferfactor", key, 'transferfactor', "avgPV")

plot(df['avgPV'], df['fr']*100, df['fr_err']*100, "fakerate", "2017", 'fakerate in %', "avgPV")
plot(df['avgPV'], df['tf'], df['tf_err'], "transferfactor", "2017", 'transferfactor', "avgPV")
