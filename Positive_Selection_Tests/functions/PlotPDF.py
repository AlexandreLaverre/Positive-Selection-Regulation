#!/usr/bin/env python
# coding=utf-8
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def f(x):
    scale = 1
    return stats.norm.logpdf(np.log10(x), scale=scale, loc=0.0)


x_list = np.logspace(-4, 4, 100)
y_list = [f(x) for x in x_list]
plt.xscale("log")
plt.axvline(1)
plt.plot(x_list, y_list)
plt.show()
