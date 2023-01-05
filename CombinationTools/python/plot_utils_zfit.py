textsize = 15
markersize = 4

### --- set plotting style
def set_style():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Palatino",],
        "font.size": textsize
        # 'text.latex.preamble': [r"""\usepackage{bm}"""]
    })

    mpl.rcParams.update({
        "legend.fontsize" : "medium",
        "axes.labelsize" : "medium",
        "axes.titlesize" : "medium",
        "xtick.labelsize" : "medium",
        "ytick.labelsize" : "medium",
    })


### --- plot correlation/covariance matrix
def plot_matrix(matrix, labels=[], matrix_type="correlation", name="", outDir="./"):
    import matplotlib.pyplot as plt
    import numpy as np

    print("Make {0} matrix plot".format(matrix_type))

    set_style()


    fig, ax = plt.subplots()
    im = ax.imshow(matrix)

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    if matrix_type == "correlation":
        nround = lambda x: round(x, 3)
    else:
        nround = lambda x: round(x, 2)

    # Loop over data dimensions and create text annotations.
    for i in range(len(labels)):
        for j in range(len(labels)):

            val = nround(matrix[i, j])
            text = ax.text(j, i, val, ha="center", va="center", color="w")

    # ax.set_title("Harvest of local farmers (in tons/year)")
    fig.tight_layout()
    name = name if name.startswith("_") else "_"+name
    suffix = "_"+matrix_type
    plt.savefig("{0}/matrix{1}{2}.png".format(outDir, name, suffix))


### --- nuisance parameters pulls and constraints
def plot_pulls(result, outDir="./"):
    import matplotlib.pyplot as plt
    import numpy as np

    print("Make pulls plot")

    set_style()

    names = np.array([p.name for p in result.params])
    xx = np.array([result.params[p]["value"] for p in result.params])
    xx_hi = np.array([result.params[p]["errors"]["upper"] for p in result.params])
    xx_lo = np.array([result.params[p]["errors"]["lower"] for p in result.params])
    yy = np.arange(len(names))

    # parameters with a name starting with "r_" are rate parameters
    is_rate = np.array([True if n.startswith("r_") else False for n in names])

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=0.0, left=0.4, right=0.97, top=0.99, bottom=0.125)

    ymin = yy[0]-1
    ymax = yy[-1]+1

    ax.plot([1, 1], [ymin,ymax], linestyle="dashed", color="gray")
    ax.plot([-1, -1], [ymin,ymax], linestyle="dashed", color="gray")

    ax.plot([2, 2], [ymin,ymax], linestyle="dashed", color="gray")
    ax.plot([-2, -2], [ymin,ymax], linestyle="dashed", color="gray")


    nround = lambda x: round(x, 3)

    for i, r in enumerate(is_rate):
        if r:
            ax.text(-0.5, yy[i]-0.4, "$"+str(nround(xx[i]))+"^{+"+str(nround(xx_hi[i]))+"}_{"+str(nround(xx_lo[i]))+"}$")

    # only plot nuisance parameters that are constraint (no rate parameters)
    xx = xx[~is_rate]
    xx_hi = xx_hi[~is_rate]
    xx_lo = xx_lo[~is_rate]
    yy = yy[~is_rate]

    ax.errorbar(xx, yy, xerr=(abs(xx_lo), abs(xx_hi)), fmt="ko", ecolor='black', elinewidth=1.0, capsize=1.0, barsabove=True, markersize=markersize)
    ax.set_xlabel("($\\hat{\\Theta} - \\Theta_0 ) / \\Delta \\Theta$")
    ax.set_ylabel("")

    ax.set_yticks(np.arange(len(names)), labels=names)
    ax.set_xlim(-2.5,2.5)
    ax.set_ylim(ymin, ymax)

    plt.savefig("{0}/pulls.png".format(outDir))


### --- plot likelihood scan
def plot_scan(result, loss, minimizer, param, name="param", profile=True, limits=2., outDir="./"):
    import matplotlib.pyplot as plt
    import numpy as np
    import zfit

    set_style()

    print("Make likelihood scan for {0}".format(name))

    val = result.params[param]["value"]
    err = result.params[param]["hesse"]["error"]

    # scan around +/- 2 sigma (prefit) interval
    xLo = val-limits
    xHi = val+limits

    x = np.linspace(xLo, xHi, num=100)
    y = []
    param.floating = False
    for val in x:
        param.set_value(val)
        if profile:
            minimizer.minimize(loss)
        y.append(loss.value())

    y = (np.array(y) - result.fmin)*2

    param.floating = True
    zfit.param.set_values(loss.get_params(), result)

    # ymin = min(y)

    yargmin = np.argmin(y)
    xmin = x[np.argmin(y)]
    # left and right 68% intervals
    xL1 = x[np.argmin(abs(y[:yargmin]-1))]
    xR1 = x[yargmin+np.argmin(abs(y[yargmin:]-1))]
    # left and right 95% intervals
    xL2 = x[np.argmin(abs(y[:yargmin]-4))]
    xR2 = x[yargmin+np.argmin(abs(y[yargmin:]-4))]

    plt.figure()
    
    if max(y) >= 1:
        plt.plot([xLo, xL1], [1,1], linestyle="dashed", color="gray")
        plt.plot([xR1, xHi], [1,1], linestyle="dashed", color="gray")
        plt.plot([xL1, xL1], [0,1], linestyle="dashed", color="gray")
        plt.plot([xR1, xR1], [0,1], linestyle="dashed", color="gray")

    if max(y) >= 4:
        plt.plot([xLo, xL2], [4,4], linestyle="dashed", color="gray")
        plt.plot([xR2, xHi], [4,4], linestyle="dashed", color="gray")
        plt.plot([xL2, xL2], [0,4], linestyle="dashed", color="gray")
        plt.plot([xR2, xR2], [0,4], linestyle="dashed", color="gray")

    plt.plot(x, y, color="red")
    plt.xlabel(name)
    plt.ylabel("-2 $\\Delta$ ln(L)")

    plt.xlim(xLo,xHi)
    plt.ylim(0,5)

    plt.tight_layout()

    suffix = "_profiling" if profile else "_scan"
    plt.savefig("{0}/loss_{1}{2}.png".format(outDir, name, suffix))
