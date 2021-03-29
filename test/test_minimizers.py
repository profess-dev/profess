# Copyright (c) 2019-2020 William C. Witt
# SPDX-License-Identifier: MIT

import numpy as np
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), '../build/'))
import extras as profess

class TestMinimizers(unittest.TestCase):

    def test_two_point_steepest_descent(self):

        A = np.array([[20,  0,  0,  0],
                      [ 0, 10,  0,  0],
                      [ 0,  0,  2,  0],
                      [ 0,  0,  0,  1]])
        b = np.array([1, 1, 1, 1])
        def func(x):
            xa = np.asarray(x)
            f = 0.5*xa.dot(A.dot(xa)) - b.dot(xa)
            g = A.dot(xa) - b
            return f, g
        x0 = np.array([0, 0, 0, 0])
        # test with max_iter=10
        conv, f, x, iters = \
            profess.two_point_steepest_descent(func, x0, 1e-10, 4, 1e-7, 10)
        self.assertTrue(conv == False)
        # test with max_iter=default
        conv, f, x, iters = \
            profess.two_point_steepest_descent(func, x0, 1e-10, 4, 1e-7)
        self.assertTrue(conv == True)
        self.assertAlmostEqual(f, -0.825)
        self.assertTrue(np.allclose(x, np.array([0.05, 0.1, 0.5, 1.0])))
        self.assertEqual(iters, 24)

if __name__ == '__main__':

    unittest.main()
