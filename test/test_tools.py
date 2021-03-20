# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/'))
import extras as profess

class TestTools(unittest.TestCase):

    def test_interpolate_cubic_lagrange(self):

        # define function and generate samples
        def f(x):
            return np.cos(x) * np.exp(-0.01*x*x)
        x, d = np.linspace(0.0, 17.5, 50, retstep=True)
        y = f(x)
        # define fine grid and compute y_interp and y_exact
        xx = np.linspace(0.0, 17.5, 100, endpoint=False)
        y_interp = np.zeros(len(xx))
        for i in range(len(xx)):
            y_interp[i] = profess.interpolate_cubic_lagrange(y,xx[i],0.0,d)
        y_exact = f(xx)
        # make plot
        if False:
            plt.plot(x, y, 'o')
            plt.plot(xx, y_interp, '.')
            plt.plot(xx, y_exact)
            plt.show()
        # test y_interp against y_exact
        self.assertTrue(np.allclose(y_interp, y_exact, atol=0.01))

    def test_interpolate_derivative_cubic_lagrange(self):

        # define function and its derivative, and generate samples
        def f(x):
            return np.cos(x)*np.exp(-0.01*x*x)
        def fp(x):
            return (- np.sin(x)*np.exp(-0.01*x*x)
                        - np.cos(x)*np.exp(-0.01*x*x)*0.01*2.0*x)
        x, d = np.linspace(0.0, 17.5, 50, retstep=True)
        y = f(x)
        yp = fp(x)
        # define fine grid and compute yp_interp and yp_exact
        xx = np.linspace(0.0, 17.5, 100, endpoint=False)
        yp_interp = np.zeros(len(xx))
        for i in range(len(xx)):
            yp_interp[i] = \
                profess.interpolate_derivative_cubic_lagrange(y,xx[i],0.0,d)
        yp_exact = fp(xx)
        # make plot
        if False:
            plt.plot(x, yp, 'o')
            plt.plot(xx, yp_interp, '.')
            plt.plot(xx, yp_exact)
            plt.show()
        # test yp_interp against yp_exact
        self.assertTrue(np.allclose(yp_interp, yp_exact, atol=0.02))

if __name__ == '__main__':

    unittest.main()
