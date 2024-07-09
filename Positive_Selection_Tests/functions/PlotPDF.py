#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def f(x):
    scale = 1
    return stats.norm.logpdf(np.log10(x), scale=scale, loc=0.0)


def norm_distribution_log(x):
    # normal distribution in log scale
    return -0.5 * np.log(2 * np.pi) - 0.5 * np.power(x, 2)


def log_likelihood(x):
    return norm_distribution_log(np.log10(x))


x_list = np.logspace(-4, 4, 100)
y_list = [f(x) for x in x_list]
y_list2 = [log_likelihood(x) for x in x_list]
plt.xscale("log")
plt.axvline(1)
plt.plot(x_list, y_list)
plt.plot(x_list, y_list2)
plt.show()
