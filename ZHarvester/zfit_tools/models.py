import zfit
import zfit_physics as zphysics

# signal models -----------------------------------------
ZMASS = zfit.Parameter('mass', 91.1876, floating=False)
ZWIDTH = zfit.Parameter('width', 2.4952, floating=False)

# resolution functions -----------------------------------------

def get_func(obs, name):
    if name == "bw":
        return zphysics.pdf.RelativisticBreitWigner(obs=obs, m=ZMASS, gamma=ZWIDTH)

def get_resolution(name, category=""):
    if name == "gauss":
        return zfit.pdf.Gauss(
            obs=zfit.Space('x', (-10., 10.)), 
            mu=zfit.Parameter('sig_mu_{0}'.format(category), 0, -2.5, 2.5), 
            sigma=zfit.Parameter('sig_sigma_{0}'.format(category), 2, 0.1, 5))
    elif name == "cb":
        return zfit.pdf.CrystalBall(
            obs=zfit.Space('x', (-10., 10.)), 
            mu=zfit.Parameter('sig_mu_{0}'.format(category), 0, -2.5, 2.5), 
            sigma=zfit.Parameter('sig_sigma_{0}'.format(category), 2, 0.1, 5),
            alpha=zfit.Parameter('sig_alpha_{0}'.format(category), 5, 0, 20),
            n=zfit.Parameter('sig_n_{0}'.format(category), 5,0.5,10),
            )        


def get_signal(obs, func="bw", resolution="gauss", category=""):
    
    func = get_func(obs, func)

    if resolution == None:
        return func
    
    return zfit.pdf.FFTConvPDFV1(func, get_resolution(resolution, category), obs=obs, interpolation="spline:3")

# background models -----------------------------------------

def get_background(obs, func="exp", category=""):

    if func == "uniform":
        return zfit.pdf.Uniform(obs=obs)
    elif func == "exp":
        return zfit.pdf.Exponential(zfit.Parameter('bkg_lambda_{0}'.format(category), -0.01, -0.5, 1.0), obs=obs)
    elif func == "chebyshev":
        return zfit.pdf.Chebyshev(coeffs=[zfit.Parameter('bkg_{0}_{1}'.format(i, category), 0.0, -1.0, 1.0) for i in range(3)], obs=obs)
