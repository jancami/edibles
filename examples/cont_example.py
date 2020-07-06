import numpy as np
import matplotlib.pyplot as plt

from edibles.models import ContinuumModel

# #################################################################################
# Example 1

x = np.linspace(0, 3)
y = x**3 - 3 * x**2 + 1

cont_model = ContinuumModel(n_anchors=4)
cont_pars = cont_model.guess(y, x=x)

# ##############################
# Show initial model

out = cont_model.eval(data=y, params=cont_pars, x=x)

init_pars = [cont_pars['y_0'].value,
             cont_pars['y_1'].value,
             cont_pars['y_2'].value,
             cont_pars['y_3'].value]
init_x = np.linspace(0, 3, 4)

plt.scatter(x, y)
plt.plot(x, out, 'C1')
plt.scatter(init_x, init_pars, marker='x', s=80, color='k')
plt.show()

# ##############################
# Fit & show results

result = cont_model.fit(data=y, params=cont_pars, x=x)
out = cont_model.eval(data=y, params=result.params, x=x)
resid = y - out


print(result.fit_report())

result.plot_fit()
plt.scatter(init_x, result.params, marker='x', color='r', s=80, zorder=10, label='Fit params')
plt.scatter(init_x, init_pars, marker='x', s=80, color='k', label='Initial params')
plt.legend()
plt.show()

# #################################################################################
# Example 2

x = np.linspace(0, 8)
y = np.sin(x)

cont_model = ContinuumModel(n_anchors=5)
cont_pars = cont_model.guess(y, x=x)

# ##############################
# Show initial model

out = cont_model.eval(data=y, params=cont_pars, x=x)

init_pars = [cont_pars['y_0'].value,
             cont_pars['y_1'].value,
             cont_pars['y_2'].value,
             cont_pars['y_3'].value,
             cont_pars['y_4'].value]

init_x = np.linspace(0, 8, 5)

plt.scatter(x, y)
plt.plot(x, out, 'C1')
plt.scatter(init_x, init_pars, marker='x', s=80, color='k')
plt.show()

# ##############################
# Fit & show results

result = cont_model.fit(data=y, params=cont_pars, x=x)
out = cont_model.eval(data=y, params=result.params, x=x)
resid = y - out


print(result.fit_report())

result.plot_fit()
plt.scatter(init_x, result.params, marker='x', color='r', s=80, zorder=10, label='Fit params')
plt.scatter(init_x, init_pars, marker='x', s=80, color='k', label='Initial params')
plt.legend()
plt.show()
