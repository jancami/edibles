import astropy.constants as cst


def convertToVelocity(data):


# this could use the edibles_spectrum object



    wave, flux = data

    delta_lambda = 

    vel_wave = delta_lambda / wave * cst.c.to('km/s').value