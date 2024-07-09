#!/usr/bin/env python
# coding=utf-8
import numpy as np
from scipy.stats import beta
import matplotlib.pyplot as plt
from BayesEvolNumba import log_scaling_param as numpy_log_scaling_param
from BayesEvol import log_scaling_param as scipy_log_scaling_param

x_list = np.logspace(-4, 4, 100)
y_list = [scipy_log_scaling_param(x) for x in x_list]
y_list2 = [numpy_log_scaling_param(x) for x in x_list]
plt.xscale("log")
plt.axvline(1, color="black", linestyle="--")
plt.plot(x_list, y_list, label="Scipy")
plt.plot(x_list, y_list2, label="Numba", linestyle="--")
plt.show()

from BayesEvolNumba import beta_distribution as numpy_beta_distribution
x_list = np.linspace(0, 1, 100)
a, b = 0.5, 2
y_list = [beta.pdf(x, a, b) for x in x_list]
y_list2 = [numpy_beta_distribution(x, a, b) for x in x_list]
plt.plot(x_list, y_list, label="Scipy")
plt.plot(x_list, y_list2, label="Numba", linestyle="--")
plt.axvline(0.5, color="black", linestyle="--")
plt.legend()
plt.show()
