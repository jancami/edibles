

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model


def custom(b,c):
    amp, cen, wid = 1, 1, 0.3
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(b*c-cen)**2 / (2*wid**2))


mu, sigma = 0, 0.1 # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
b = np.linspace(-1,3,1000)
c = 5
y = custom(b, c)+s
plt.plot(b, y)

gmodel = Model(custom)
result = gmodel.fit(y, b=b, c=c, independent_vars=[ 'b','c']) #, amp=2, cen=1.2, wid=1)

print(result.fit_report())

plt.plot(b, y)
plt.plot(b, result.init_fit, '--', label='initial fit')
plt.plot(b, result.best_fit, '-', label='best fit')
plt.legend()
plt.show()

# df = pd.DataFrame({
#   'A'      : pd.Series([1, 1, 1, 2, 2, 2, 2]),
#   'B'      : pd.Series([5, 4, 6, 6, 5, 6, 5]),
#   'target' : pd.Series([87.79, 40.89, 215.30, 238.65, 111.15, 238.65, 111.15])
# })

# def fun(A, B, p1 = 1, p2 = 1):
#   return p1 * np.exp(A) + p2 * np.exp(B)

# model = Model(fun, independent_vars=['A', 'B'])
# fit = model.fit(df['target'], A = df['A'], B = df['B'])


# print(fit.fit_report())

# plt.plot(x, y)
# plt.plot(b, result.init_fit, '--', label='initial fit')
# plt.plot(b, result.best_fit, '-', label='best fit')
# plt.legend()
# plt.show()















