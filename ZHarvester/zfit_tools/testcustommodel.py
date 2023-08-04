import uproot
import matplotlib.pyplot as plt
import vector
from hist import Hist
import numpy as np
import tensorflow as tf
import math

import zfit
from zfit import z

import tensorflow.keras.backend as K

import pdb


class MyGauss(zfit.pdf.ZPDF):
    # _N_OBS = 1  # dimension, can be omitted
    _PARAMS = ['mean1', 'std1', 'mean2', 'std2']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean1 = self.params['mean1']
        std1 = self.params['std1']
        mean2 = self.params['mean2']
        std2 = self.params['std2']
        return z.exp(- ((x - mean1) / std1) ** 2) + z.exp(- ((x - mean2) / std2) ** 2)

class MyGaussBinned(zfit.pdf.ZPDF):
    _PARAMS = ['mean', 'std']  # the name of the parameters

    def _unnormalized_pdf(self, x):
        x = z.unstack_x(x)  # returns a list with the columns: do x, y, z = z.unstack_x(x) for 3D
        mean = self.params['mean']
        std = self.params['std']

        return binned_Gaussian_template(mean, std, -5, 5, 50)


def binned_Gaussian_template(mu, sigma, start, end, bins):
    x_np = np.linspace(start, end, num=bins)
    x = tf.convert_to_tensor(x_np)
    # x = tf.cast(x, tf.float32)
    pi = math.pi
    y = (1 / (sigma*math.sqrt(2 * pi))) * tf.math.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    # hist = Hist.new.Regular(bins, start, end, name = 'x').Double()
    # for i in range(bins):
    #     hist[i] = tf.gather(y,i)
        
    # return zfit.data.BinnedData.from_hist(hist)

    binning = zfit.binned.RegularBinning(bins, start, end, name="x")
    space = zfit.Space("x", binning=binning)
    return zfit.data.BinnedData.from_tensor(space, y)
    # return y


obs = zfit.Space('x', limits=(-10, 10))
binning = zfit.binned.RegularBinning(50, -5, 5, name="x")
obs_binned = zfit.Space("x", binning=binning)

# data
mu_true = 2
sigma_true = 1

data_np = np.random.normal(mu_true, sigma_true, size=1000)
data = zfit.data.Data.from_numpy(obs=obs, array=data_np)
data_binned = zfit.data.BinnedData.from_unbinned(obs_binned, data)
# pdb.set_trace()


# define params
mu = zfit.Parameter("mu", 2.4, -1., 5., step_size=0.001)
sigma = zfit.Parameter("sigma", 1.3, 0, 5., step_size=0.001)

my_Gauss_2 = MyGaussBinned(obs=obs_binned, mean=mu, std=sigma)

pdb.set_trace()
my_Gauss_binned = zfit.pdf.BinnedFromUnbinnedPDF(my_Gauss_2, obs_binned)

loss = zfit.loss.BinnedNLL(my_Gauss_binned, data_binned)
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(loss)



pdb.set_trace()
