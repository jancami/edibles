import numpy as np
import matplotlib.pyplot as plt

from edibles.models import ContinuumModel
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from edibles.utils.edibles_oracle import EdiblesOracle

# #################################################################################
# Example 1

x = np.linspace(0, 3)
y = x**3 - 3 * x**2 + 1

cont_model = ContinuumModel(n_anchors=4)
cont_pars = cont_model.guess(y, x=x)


# ##############################
# Show initial model

out = cont_model.eval(data=y, params=cont_pars, x=x)

y_param_names = []
for i in range(cont_model.n_anchors):
    y_param_names.append('y_' + str(i))

x_param_names = []
for i in range(cont_model.n_anchors):
    x_param_names.append('x_' + str(i))

init_y = []
for par_name in y_param_names:
    init_y.append(cont_pars[par_name].value)

init_x = []
for par_name in x_param_names:
    init_x.append(cont_pars[par_name].value)

plt.scatter(x, y)
plt.plot(x, out, 'C1')
plt.scatter(init_x, init_y, marker='x', s=80, color='k')
plt.show()


# ##############################
# Fit

result = cont_model.fit(data=y, params=cont_pars, x=x)
out = cont_model.eval(data=y, params=result.params, x=x)
# print(result.fit_report())


# ##############################
# Show results

result_y_pars = []
for par_name in y_param_names:
    result_y_pars.append(result.params[par_name].value)

result_x_pars = []
for par_name in x_param_names:
    result_x_pars.append(result.params[par_name].value)

result.plot_fit()
plt.scatter(result_x_pars, result_y_pars, marker='x',
            color='r', s=80, zorder=10, label='Fit params')
plt.scatter(init_x, init_y, marker='x', s=80, color='k', label='Initial params')
plt.legend()
plt.show()


# #################################################################################
# Example 2
# #################################################################################
pythia = EdiblesOracle()
List = pythia.getFilteredObsList(object=["HD 23466"], OrdersOnly=True, Wave=3290)
filename = List.values.tolist()[0]
print(filename)
sp = EdiblesSpectrum(filename)
# sp = EdiblesSpectrum("/HD23466/BLUE_346/HD23466_w346_blue_20180731_O11.fits")
print(sp.wave)
xmin = 3270
xmax = 3305

sp.getSpectrum(xmin=xmin, xmax=xmax)


cont_model = ContinuumModel(n_anchors=4)
cont_pars = cont_model.guess(sp.flux, x=sp.wave)

# ##############################
# Show initial model

out = cont_model.eval(data=sp.flux, params=cont_pars, x=sp.wave)

y_param_names = []
for i in range(cont_model.n_anchors):
    y_param_names.append('y_' + str(i))

x_param_names = []
for i in range(cont_model.n_anchors):
    x_param_names.append('x_' + str(i))

init_y = []
for par_name in y_param_names:
    init_y.append(cont_pars[par_name].value)

init_x = []
for par_name in x_param_names:
    init_x.append(cont_pars[par_name].value)

plt.scatter(sp.wave, sp.flux)
plt.plot(sp.wave, out, 'C1')
plt.scatter(init_x, init_y, marker='x', s=80, color='k')
plt.show()


# ##############################
# Fit

result = cont_model.fit(data=sp.flux, params=cont_pars, x=sp.wave)
out = cont_model.eval(data=sp.flux, params=result.params, x=sp.wave)
# print(result.fit_report())


# ##############################
# Show results

result_y_pars = []
for par_name in y_param_names:
    result_y_pars.append(result.params[par_name].value)

result_x_pars = []
for par_name in x_param_names:
    result_x_pars.append(result.params[par_name].value)

result.plot_fit()
plt.scatter(result_x_pars, result_y_pars, marker='x',
            color='r', s=80, zorder=10, label='Fit params')
plt.scatter(init_x, init_y, marker='x', s=80, color='k', label='Initial params')
plt.legend()
plt.show()
