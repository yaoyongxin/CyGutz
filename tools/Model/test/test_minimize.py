import numpy as np


def fun(x, a, b):
    print x
    res = (x[0]-a)**2 + (x[1]-b)**4
    jac = np.asarray([2*(x[0]-a)**1, 4*(x[1]-b)**3])
    return res, jac

from scipy.optimize import minimize
x0 = [0.1, -0.2]
res = minimize(fun, x0, args=(0., 0.), jac=True, method='CG')
print res
