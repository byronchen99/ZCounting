from pandas import DataFrame
import ROOT
import numpy as np
import pickle
import pdb
import uncertainties as unc
from scipy.stats import poisson

# first order polynomial
def linear(x, a, b):
    return a * x + b

# two first order polynomial with step
def linear_step(x, a1, b1, a2, b2, step=35):
    return (a1 * x + b1) * (x <= step) + (a2 * x + b2) * (x > step)

# second order polynomial
def quad(x, a, b, c):
    return a * x**2 + b * x + c
# exponential
def exp(x, a, b, c, d):
    return a + b * 2**( c * x + d)

# functions folded with poisson

def pexp(l, a, b, c, d):
    p = lambda x: poisson.pmf(x, l)
    f = lambda x: a + b * 2**( c * x + d)

    return sum([p(i) * f(i) for i in range(200)])

def pquad(l, a, b, c):
    p = lambda x: poisson.pmf(x, l)
    f = lambda x: a * x**2 + b * x + c

    return sum([p(i) * f(i) for i in range(200)])

def plinear(l, a, b):
    p = lambda x: poisson.pmf(x, l)
    f = lambda x: a * x + b

    return sum([p(i) * f(i) for i in range(200)])

def plinear_step(l, a, b, c, step):
    p = lambda x: poisson.pmf(x, l)
    f = lambda x: a * x + b + c * (x-step) * (x > step)

    return sum([p(i) * f(i) for i in range(200)])

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
            if arrSize is not None:
                for key in df.keys():
                    if df[key].dtype =='O':
                        for i in range(arrSize):
                            df[key+"_"+str(i)] = df[key].apply(lambda x: x[i] if len(x) > i else float('NaN'))
                        df = df.drop([key],axis=1)
    return df

def to_RootTime(time, currentYear):
    # converts brilcalc time to root TDatime
    time =  time.split(" ")
    return ROOT.TDatime(currentYear, int(time[0].split("/")[0]),
                         int(time[0].split("/")[1]), int(time[1].split(":")[0]),
                         int(time[1].split(":")[1]), int(time[1].split(":")[2])).Convert() + 7200

def getMCCorrection(fIn):
    # input file has to be a picked dictionary with
    #   eta region, function name and parameters for the function
    #   supported functions are linear and linear_step

    functions = pickle.load(open(fIn,"r"))
    corrections = {}
    for key, val in functions.items():
        if val["name"] == "linear":
            corrections[key] = linear(*val["params"])
        if val["name"] == "linear_step":
            corrections[key] = linear_step(*val["params"])

    return corrections

def unorm(x):
    # for counting experiments: define ufloat with poisson uncertainty
    if x > 0.:
        return unc.ufloat(x, np.sqrt(abs(x)))
    else:
        return unc.ufloat(0.0, 0.0)

latex = ROOT.TLatex()
latex.SetNDC()
def workinprogress(x=0.23, y=0.88, space=0.1, textsize=0.04, offset_x=0.13):
    cms(x, y, textsize)
    latex.SetTextAlign(11)
    latex.SetTextFont(52)
    latex.SetTextSize(textsize)
    latex.DrawLatex(x+space, y, 'Work in progress')
def cms(x=0.13, y=0.88, textsize=0.04):
    latex.SetTextAlign(11)
    latex.SetTextFont(61)
    latex.SetTextSize(textsize*1.2)
    latex.DrawLatex(x, y, 'CMS')
def preliminary(x=0.23, y=0.88, space=0.1, textsize=0.04, align=11):
    cms(x, y, textsize)
    latex.SetTextAlign(11)
    latex.SetTextFont(52)
    latex.SetTextSize(textsize)
    latex.DrawLatex(x+space, y, 'Preliminary')
def simulation(x=0.23, y=0.88, textsize=0.04):
    cms(x-0.08*textsize/0.04, y, textsize)
    latex.SetTextAlign(11)
    latex.SetTextFont(52)
    latex.SetTextSize(textsize)
    latex.DrawLatex(x, y, 'Simulation')
def text(text, x, y, color=1, textsize=0.04, align=11, textfont=42):
    latex.SetTextAlign(align)
    latex.SetTextFont(textfont)
    latex.SetTextColor(color)
    latex.SetTextSize(textsize)
    latex.DrawLatex(x, y, text)

def custom_labels_y(_h, _labels, _rotation=-1, _align=-1, _offset=None,):

    _h.GetYaxis().SetNdivisions(len(_labels))

    if _offset:
        _h.GetYaxis().SetLabelOffset(_offset)

    for i in range(len(_labels)):
        _h.GetYaxis().ChangeLabel(i+1,_rotation,-1,_align,-1,-1,_labels[i])
