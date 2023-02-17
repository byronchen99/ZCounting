import matplotlib as mpl
import matplotlib.pyplot as plt

def set_matplotlib_style():

    textsize = 16
    markersize = 4.0
    labelsize =  12.5

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino",],
        "font.size": textsize,
        'text.latex.preamble': [r"""\usepackage{bm}"""]
    })

    mpl.rcParams.update({
        "legend.fontsize" : "medium",
        "axes.labelsize" : "medium",
        "axes.titlesize" : "medium",
        "xtick.labelsize" : "medium",
        "ytick.labelsize" : "medium",
    })
    colors = ["#e74c3c","#f1c40f","#16a085","#27ae60","#2980b9","#8e44ad"]

    return colors, textsize, labelsize, markersize