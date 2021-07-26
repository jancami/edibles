import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import lmfit
from lmfit.model import save_modelresult
from edibles.utils.simulations.SRC.Functions import Signal_Noise_Calculator, WavelengthToWavenumber
from scipy.interpolate import griddata
from edibles.utils.simulations.SimulatedContour import Simulated_Contour as sim
from timeit import default_timer as timer


class Fitted_Contours:

    def __init__(self, DIB, Target, Initial_Params, J_Limit, Q_Branch=True):
        filename = 'SpectraData/'+str(DIB)+'/'+str(Target)+'_avg_spectra.csv'
        df = pd.read_csv(filename)
        self.DIB = DIB
        self.Target = Target
        self.wavelength = df['Wavelength'].to_numpy()
        self.flux = df['Flux'].to_numpy()
        self.parms = Initial_Params
        self.J_Limit = J_Limit
        self.Q_Branch = Q_Branch

    def calc_uncertainty(self):
        cont_x1 = np.array(self.wavelength[-20:])
        cont_y1 = np.array(self.flux[-20:])
        SN1, Fit1 = Signal_Noise_Calculator(cont_x1, cont_y1)

        cont_x2 = np.array(self.wavelength[:15])
        cont_y2 = np.array(self.flux[:15])
        SN2, Fit2 = Signal_Noise_Calculator(cont_x2, cont_y2)

        # Average of the two S/N
        Signal_Noise = 0.5*(SN1+SN2)
        # Create an array of flux uncertainty values based on calculated S/N
        New_Noise = self.flux/Signal_Noise

        self.uncertainty_array = np.full(self.wavelength.shape, np.mean(New_Noise), dtype=float)

    def _residual(self, params, Data_X, Data_Y, Yerr, Jlimit=50, Symmetry_Type=None, Sightline=None, Q_Branch=False):

        # read in the dict of parameters. Apply appropriate constraints depending on symmetry type

        parvals = params.valuesdict()
        A = parvals['B']
        B = parvals['B']
        C = parvals['B']
        Q_scale = parvals['Q_Scale']
        PR_scale = parvals['PR_Scale']
        T = parvals['T']
        delta = parvals['Delta']
        lambda0 = parvals['Lambda_0']
        Y_Scale = parvals['Y_Scale']

        build = sim(A=A, B=B, C=C, Delta_A=delta*A, Delta_B=delta*B, Delta_C=delta*C, Trot=T, Jlimit=Jlimit,
                    Target=Sightline, lambda0=lambda0, Q_Branch=self.Q_Branch, Q_scale=Q_scale, PR_scale=PR_scale)
        old_x = np.asarray(build[0])

        old_y = np.asarray(build[1])

        new_y = griddata(old_x, old_y, Data_X, method='cubic', fill_value=1.0)
        new_y = (new_y/np.max(new_y*Y_Scale))*Y_Scale

        if self.uncertainty_array is None:
            self.residual = new_y-Data_Y
        else:
            self.residual = (new_y-Data_Y)/Yerr
        print(A, B, C, Q_scale, T, delta, lambda0)

#        plt.plot(Data_X,new_y,'c-')
#        plt.plot(Data_X,Data_Y,'b-')
#        plt.show()
        return(self.residual)

    def minimize_contour(self, verbose=True):

        self.results = lmfit.minimize(self._residual, params=self.parms, kws={'Data_X': self.wavelength, 'Data_Y': self.flux, 'Yerr': self.uncertainty_array, 'Jlimit': self.J_Limit,
                                      'Sightline': self.Target, 'Q_Branch': self.Q_Branch}, method='nelder', nan_policy='omit', options={'disp': True, 'fatol': 0.01, 'xatol': 0.01})

        f = open('SavedModels/'+str(self.DIB)+'_'+str(self.Target)+'_Lambda_' +
                 str(self.results.params['Lambda_0'].value)+'.txt', "w")
        f.write(str(self.results.params))
        f.close()

        if verbose == True:
            print(self.results.status)
            print(lmfit.fit_report(self.results))

    def plot_results(self):

        B = self.results.params['B'].value
        delta = self.results.params['Delta'].value
        T_r = self.results.params['T'].value
        Q_s = self.results.params['Q_Scale'].value
        Y_s = self.results.params['Y_Scale'].value
        P_s = self.results.params['PR_Scale'].value
        DIB_l = self.results.params['Lambda_0'].value

        build = sim(A=B, B=B, C=B, Delta_A=delta*B, Delta_B=delta*B, Delta_C=delta*B, Trot=T_r,
                    Jlimit=50, Q_scale=Q_s, PR_scale=P_s, Target=self.Target, lambda0=DIB_l, Q_Branch=True)

        new_y = griddata(build[0], build[1], self.wavelength, method='cubic', fill_value=1.0)
        new_y = (new_y/np.max(new_y*Y_s))*Y_s

        plt.title('Best Fit'+str(self.results.params))
        plt.plot(self.wavelength, self.flux, label='Original Data')
        plt.plot(self.wavelength, new_y, label='Best Fit')
        plt.legend()
        plt.savefig('SavedModels/'+str(self.DIB)+'_'+str(self.Target) +
                    '_Lambda_'+str(self.results.params['Lambda_0'].value)+'.pdf')
        # plt.show()
        plt.close()
################################################################################


#####
def Fit_Contour(DIB, Target, Initial_Params, J_Limit, Q_Branch=True):
    build = Fitted_Contours(DIB, Target, Initial_Params, J_Limit, Q_Branch)
    build.calc_uncertainty()
    build.minimize_contour()
    build.plot_results()
    print('hi')


################################################################################
Sightlines = ['HD144470', 'HD147165', 'HD147683', 'HD149757', 'HD166937',
              'HD170740', 'HD184915', 'HD185418', 'HD185859', 'HD203532', 'HD23180', 'HD24398']

################################################################################
# Create Initial Params##
p = lmfit.Parameters()
p.add('Q_Scale', 0.3, min=0.0)
p.add('PR_Scale', 0.5, min=0.0)
p.add('Y_Scale', 0.005, min=0)
p.add('T', 20, min=0, max=50)
p.add('Delta', 0.002)
p.add('B', 30e-3, min=10e-4, max=10e-1)
for Target in Sightlines:
    for x in range(5):
        start = timer()
        print(start)
        p.add('Lambda_0', 6613.23+(0.01*x), vary=False)

        Fit_Contour(6614, Target, p, 50)
        end = timer()
        print(end)
        print(end-start, ' seconds')
