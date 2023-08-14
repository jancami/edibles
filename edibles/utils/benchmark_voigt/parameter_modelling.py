import numpy as np
import matplotlib.pyplot as plt
from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd
import voigt_profile as vp
import os.path

def file_reader (star_name):
    """
    Will read in any file from the voigt_benchmark folder once given the name of the file.
    Is used in order to read in data for o per, sigma sco, zeta per, and zeta oph stars in example 5
    :paramaters :
        star_name: string
            the file name containing the data that needs to be analysed
    :return:
        data["Wavelength"] : 1d ndarray
            Contains all the measured wavelength values for a particular file
        data["Normfluxval"] : 1d ndarray
            Contains normailised flux values corresponding to each wavelength reading
    """
    # define where the desired file is
    folder = Path(PYTHONDIR + "/data")
    filename = folder / "voigt_benchmarkdata" / star_name

    # state what headers the desired data is under
    Headers = ["Wavelength", "Normfluxval"]

    # read in the data
    data = pd.read_csv(
        filename,
        delim_whitespace=True,
        skiprows=[0],
        header=None,
        names=Headers,
        engine="python",
    )
    return data["Wavelength"], data["Normfluxval"]


def reduced_chi_squared(observed_value, observed_error, expected_value):
    """
    calculates reduced chi squared of data entered

    parameters
    ----------
    observed_value : 1-D numpy array
        value from the real life datd
    observed_error : 1-D numpy array
        error on real data
    expected_value : 1-D numpy array
        value at ths point from fitted function

    returns
    ----------
   (np.sum(((observed_value - expected_value)/observed_error)**2))/(len(observed_value)-2
       reduced chi squared
    """
    return (np.sum(((observed_value - expected_value) / observed_error) ** 2)) / (len(observed_value) - 2)



#define path to file which contains the names of each star and the resolution of intraments used in survey
folder = Path(PYTHONDIR + "/data")
filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / "files_for_parameter_modelling.txt"

# state what headers the desired data is under
Headers = ["star_name", "file_name", "star_file", "resolution"]

# read in the data
file = pd.read_csv(
    filename,
     delim_whitespace=True,
     header=None,
     names=Headers,
     engine="python",
    )

# collect the data from each header in easy to iterate form
files = file['file_name']
resolution = file['resolution']
star_data = file["star_file"]
star_name = file["star_name"]


if __name__ == '__main__' :

    for i in range(len(files)):
    #for i in range(1):

        # use file names from files_for_parameter_modelling.txt to read in flux and wavelength data
        wavelengths, normflux = np.asarray(file_reader(files[i]))

        folder = Path(PYTHONDIR + "/data")
        filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / star_data[i]


        Headers = ["b", "N", "v_rad"]
        star_parameters =pd.read_csv(filename,
         delim_whitespace = True,
         header=None,
         names=Headers,
         engine="python",
        )

        # collect the data from each header in easy to iterate form
        b = np.asarray(star_parameters["b"])
        N = np.asarray(star_parameters["N"])
        v_rad = np.asarray(star_parameters["v_rad"])

        #calculating the normflux of pottasium as predicted using the parameters in the Welty & Hobbs 2001 paper
        WH_flux= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=b,
            N=N,
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=v_rad,
            v_resolution=resolution[i],
        )
        # calculating the flux contribution of each component (mostly used for bug checking)
        WH_flux_1= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=b[0],
            N=N[0],
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=v_rad[0],
            v_resolution=resolution[i],
        )
        WH_flux_2= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=b[1],
            N=N[1],
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=v_rad[1],
            v_resolution=resolution[i],
        )
        WH_flux_3= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=b[2],
            N=N[2],
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=v_rad[2],
            v_resolution=resolution[i],
        )
        WH_flux_4= vp.voigt_absorption_line(
            wavelengths,
            lambda0=7698.974,
            b=b[3],
            N=N[3],
            f=3.393e-1,
            gamma=3.8e7,
            v_rad=v_rad[3],
            v_resolution=resolution[i],
        )
        #WH_flux_5= vp.voigt_absorption_line(
        #    wavelengths,
        #    lambda0=7698.974,
        #    b=b[4],
        #    N=N[4],
        #    f=3.393e-1,
        #    gamma=3.8e7,
        #    v_rad=v_rad[4],
        #    v_resolution=resolution[i],
        #)

        # defining continuum as all points where normflux is less than 0.98
        # this is then used to calculate errors to use for chi squared
        continuum = np.array([])
        for j in range(len(normflux)):
            if normflux[j] > 0.98:
                continuum = np.append(continuum,normflux[j])
        error = np.array([np.std(continuum)] * len(normflux))
        #print('reduced chi squared for {0} from the WH model is {1}'.format(star_name[i],(reduced_chi_squared(normflux, error, WH_flux))))

        # calculating best fit parameters by using the Welty & Hobbs parameters as initial guesses
        fit_flux = vp.fit_multi_voigt_absorptionlines(wavegrid= wavelengths,
                                           ydata= normflux,
                                           restwave= 7698.974,
                                           f= 3.393e-1,
                                           gamma= 3.8e7,
                                           b= star_parameters['b'],
                                           N= star_parameters['N'],
                                           v_rad= star_parameters['v_rad'],
                                           v_resolution= resolution[i],
                                           n_step= 25)

        # plot all of the calculated data from above for analysis
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wavelengths, normflux, color="k", marker="D", fillstyle='none', label= 'Data')
        ax.plot(wavelengths, fit_flux.best_fit, c='r', marker = '*', label= 'best-fit')
        ax.plot(wavelengths, WH_flux_1, c='b', label='WH 1st component')
        ax.plot(wavelengths, WH_flux_2, label='WH 2nd component')
        ax.plot(wavelengths, WH_flux_3, label='WH 3rd component')
        ax.plot(wavelengths, WH_flux_4, label='WH 4th component')
        #ax.plot(wavelengths, WH_flux_5, label='WH 5th component')

        ax.plot(wavelengths, WH_flux,c = 'b', label='Welty & Hobbs')
        fig.suptitle(star_name[i])
        ax.set_xlabel('Wavelength \u00C5')
        ax.set_ylabel('Normalised Flux')
        plt.legend()

        # saving Wlety & Hobbs parameters alongside best fit parameters for easy comparison (saved in parameter_modelling_data)
        save_path = 'c:/Users/user/edibles/edibles/data/voigt_benchmarkdata/parameter_modelling_data/'
        file_name = star_name[i] + '.fit_data.txt'
        file_path = os.path.join(save_path, file_name)
        file = open(file_path, 'w+')

        txt = 'Welty & Hobbs parameteres: \n b: {0}, \n N: {1}, \n v_rad: {2} \n' \
              'Fit parameters: {3}'.format(b, N, v_rad, fit_flux.fit_report())
        file.write(txt)
        file.close()


        # attempting to build code which will find out the number of components in a sightline
        # as well as find the best-fitting parameters for each component
        # currently assigns initial guesses for 1 parameter using parameters from Welty & Hobbs
        # and will then build of this for the rest of the components

        b_components = np.array([b[0]])
        N_components = np.array([N[0]])
        v_rad_components = np.array([v_rad[0]])
        components = 1
        new_fit_flux = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                                          ydata=normflux,
                                                          restwave=7698.974,
                                                          f=3.393e-1,
                                                          gamma=3.8e7,
                                                          b=b_components,
                                                          N=N_components,
                                                          v_rad=v_rad_components,
                                                         v_resolution=resolution[i],
                                                         n_step=25)
        b_components[0] = new_fit_flux.values['b0']
        N_components[0] = new_fit_flux.values['N0']
        v_rad_components[0] = new_fit_flux.values['v_rad0']
        fig_2 = plt.figure()
        ax = fig_2.add_subplot(111)
        ax.plot(wavelengths, normflux, color="k", marker="D", fillstyle='none', label='Data')
        ax.plot(wavelengths, new_fit_flux.best_fit, c='orange', marker='*', label='best-fit')
        fig_2.suptitle(star_name[i])
        ax.set_xlabel('Wavelength \u00C5')
        ax.set_ylabel('Normalised Flux')
        plt.legend()
        plt.show()

        print('b is ', b_components)
        print('N is ', N_components)
        print('v_rad is ', v_rad_components)

        # entering a while loop which will be adapted to use a fitting parameter instead of len(b)
        # which keeps adding components through every run in an attempt to find the best fit parameters
        answer = 'YES'
        #while answer.upper() == 'YES':
        count = 1
        while count < len(b) :
            print('entered while loop')
            for k in range(components):
                b_components[k] = new_fit_flux.values[f'b{k}']
                N_components[k] = new_fit_flux.values[f'N{k}']
                v_rad_components[k] = new_fit_flux.values[f'v_rad{k}']

            components += 1
            print('components',components)
            # use values which are 2 std's higher than the paramter from the first component to create a new component
            new_b = new_fit_flux.params['b0'].stderr * 2 + b_components[0]
            multiplier = 2
            print('calculating new_b')
            # issue with std being too small leading to values which are too close together breaking code so this loop
            # will increase the number of std's away until the new component is suitably different from the first one
            while (new_b - b_components[0])< 0.1 :
                multiplier += 1
                new_b = new_fit_flux.params['b0'].stderr * multiplier + b_components[0]
            print('b multiplier is:', multiplier)
            print('calculated new_b')
            b_components = np.append(b_components, new_b)

            # use values which are 2 std's higher than the paramter from the first component to create a new component
            new_N = new_fit_flux.params['N0'].stderr * 2 + N_components[0]
            N_components = np.append(N_components, new_N)

            multiplier = 2
            # use values which are 2 std's higher than the paramter from the first component to create a new component
            new_v_rad = new_fit_flux.params['v_rad0'].stderr * 2 + v_rad_components[0]
            # issue with std being too small leading to values which are too close together breaking code so this loop
            # will increase the number of std's away until the new component is suitably different from the first one
            while (new_v_rad- v_rad_components[0])< 0.1 :
                multiplier += 1
                new_v_rad = new_fit_flux.params['v_rad0'].stderr * multiplier + v_rad_components[0]
            print('v_rad multiplier is:', multiplier)
            v_rad_components = np.append(v_rad_components, new_v_rad)

            # noticed that some of the b and N values where becoming negative this is not possible so this is my tempary
            # solution to try and push these values back to being postive
            for l in range(len(b_components)):
                if b_components[l] < 0 :
                    b_components[l] = b[0]
                    print(f'reset b{l}')

            for l in range(len(N_components)):
                if N_components[l] < 0 :
                    N_components[l] = N[0]
                    print(f'reset N{l}')

            print('fitting...')
            # use the guess parameters to find best fits
            new_fit_flux = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                                              ydata=normflux,
                                                              restwave=7698.974,
                                                              f=3.393e-1,
                                                              gamma=3.8e7,
                                                              b=b_components,
                                                              N=N_components,
                                                              v_rad=v_rad_components,
                                                              v_resolution=resolution[i],
                                                              n_step=25)
            print('calculated new parameters')
            print('b is now', b_components)
            print('N is now', N_components)
            print('v_rad is now', v_rad_components)
            print('fitted')
            #fig_2 = plt.figure()
            #ax = fig_2.add_subplot(111)
            #ax.plot(wavelengths, normflux, color="k", marker="D", fillstyle='none', label='Data')
            #ax.plot(wavelengths, new_fit_flux.best_fit, c='r', marker='*', label='best-fit')
            #fig_2.suptitle(star_name[i])
            #ax.set_xlabel('Wavelength \u00C5')
            #ax.set_ylabel('Normalised Flux')
            #plt.legend()
            #plt.show()
            #answer = input('Add another component?')
            count += 1
        #plot final result in order to inspect by eye how good of a fit the parametrs are
        fig_2 = plt.figure()
        ax = fig_2.add_subplot(111)
        ax.plot(wavelengths, normflux, color="k", marker="D", fillstyle='none', label='Data')
        ax.plot(wavelengths, new_fit_flux.best_fit, c='orange', marker='*', label='best-fit')
        fig_2.suptitle(star_name[i])
        ax.set_xlabel('Wavelength \u00C5')
        ax.set_ylabel('Normalised Flux')
        plt.legend()
        plt.show()







        #table_data = np.empty((15,3), dtype = np.string_)
        #print(table_data[0,:])
        # store all parameters in table for easy analysis
        table_data = np.asarray([['                                           ' for x in range(3)] for y in range(len(b)*3)])
        row_labels = np.array(['           ' for x in range(len(b)*3)])
        col_labels = ['Welty & Hobbs', 'Fitting with W&H', 'Fitting components']
        fit_flux_array = [fit_flux, new_fit_flux]
        for x in range(3):
            if x == 0:
                for y in range(len(b)):
                    table_data[y][x] = f'{b[y]:#.3f}'
                    table_data[y+len(b)][x] = f'{N[y]:#.3g}'
                    table_data[y+len(b)*2][x] = f'{v_rad[y]:#.3f}'
                    row_labels[y] = f'b{y}'
                    row_labels[y+len(b)] = f'N{y+5}'
                    row_labels[y+len(b)*2]= f'v_rad{y+10}'
            else:
                for y in range(len(b)):
                    table_data[y][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(fit_flux_array[x-1].values[f'b{y}'], fit_flux_array[x-1].params[f'b{y}'].stderr)
                    table_data[y+len(b)][x] = '{0:#.3g} \u00B1 {1:#.3g}'.format(fit_flux_array[x-1].values[f'N{y}'],fit_flux_array[x-1].params[f'N{y}'].stderr)
                    table_data[y+len(b)*2][x] = '{0:#.3f} \u00B1 {1:#.3f}'.format(fit_flux_array[x-1].values[f'v_rad{y}'], fit_flux_array[x-1].params[f'v_rad{y}'].stderr)

        #print(fit_b)
        #print('b=',b)
        #print(table_data)
        # 'potting' the table so that the data can be viewed
        fig_3 = plt.figure()
        ax = fig_3.add_subplot(111)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        table = ax.table(table_data, rowLabels = row_labels, colLabels = col_labels, loc = 'center')
        plt.box(on=False)
        table.scale(1,1)
        fig_3.suptitle(star_name[i])
        plt.show()


#print(new_fit_flux.fit_report())
print('done')

#save all data in tables to make for easy comparison
# attempt to start parameter guessing with variable components, use uncertainties to help when adding multiple components
# one outside of while loop then always split last point