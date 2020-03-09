import numpy as np
import astropy.constants as cst
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter
import sherpa.ui as UI
import edibles.src.math as eMath
import edibles.src.datahandling as DataHandling
import math

class VoigtLine(ArithmeticModel):

    def __init__(self, name='VoigtLine'):

        self.lam_0 = Parameter(name, 'lam_0', 5000., frozen=False, min=0.0)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=0)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(self, name,
                                 (self.lam_0, self.b, self.d, self.N, self.f, self.tau_0))


    def calc(self, pars, x, *args, **kwargs):
        lam_0, b, d, N, f, tau_0 = pars

        if N != 999:
            transmission = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=lam_0, b=b, d=d, tau_0=tau_0)

        return transmission

class VoigtLine_KnownWav(ArithmeticModel):
    # use lambda_0 and v_offset, for lines with known wavelength
    # by default, the line is restrict to pm 100 km/s

    def __init__(self, name='VoigtLine', lam_0=0, v_max = 100):

        self.v_offset = Parameter(name, 'v_offset', 0., frozen=False, min=-v_max, max=v_max)
        self.lam_0 = Parameter(name, 'lam_0', lam_0, frozen=True)
        self.b = Parameter(name, 'b', 3.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=0)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(self, name,
                                 (self.v_offset, self.lam_0, self.b, self.d, self.N, self.f, self.tau_0))


    def calc(self, pars, x, *args, **kwargs):
        v_offset, lam_0, b, d, N, f, tau_0 = pars
        lam_0 = (1 + v_offset / cst.c.to('km/s').value) * lam_0
        if N != 999:
            transmission = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=lam_0, b=b, d=d, tau_0=tau_0)

        return transmission

class Sightline(ArithmeticModel):

    def __init__(self, star_name, cont):
        self.star_name = star_name
        self.model = cont
        self.resolution = 100000
        self.component=[cont]

    def addLine(self, lam_0, *kwords, line_name="Line", b=3.5, d=0.0005, tau_0=0.1):
        if "known_wave" in kwords:
            new_line = VoigtLine_KnownWav(name=line_name)
        else:
            new_line = VoigtLine(name=line_name)

        new_line.lam_0 = lam_0
        new_line.b = b
        new_line.d = d
        new_line.tau_0 = tau_0
        self.model *= new_line
        return new_line

    def addLineSeries(self, lam_0, *kwords, line_name="Line", b=3.5, d=0.0005, tau_0=0.1):
        current_component = len(self.component)
        lam_0 = DataHandling.parseInput(1, lam_0, checklen=False)
        (line_name_list, b, d, tau_0) = DataHandling.parseInput(len(lam_0), line_name, b, d, tau_0)
        if type(line_name) is str:
            for i, name in enumerate(line_name_list):
                line_name_list[i] = name + "_{i}".format(i=i)

        for i, name in enumerate(line_name_list):
            new_line = self.addLine(lam_0[i], *kwords, line_name=name, b=b[i], d=d[i], tau_0=tau_0[i])
            self.component.append(new_line)

        if "link_velocity" in kwords and "known_wave" in kwords:
            for i in range(len(line_name_list) - 1):
                UI.link(self.component[current_component+i+1].v_offset, self.component[current_component].v_offset)

        if "link_b" in kwords:
            for i in range(len(line_name_list) - 1):
                UI.link(self.component[current_component + i + 1].b, self.component[current_component].b)


class DualTelluric(ArithmeticModel):

    def __init__(self, name='DuleTelluric', same_b=False, same_d=False):
        self.Cst_Cont = Parameter(name, 'Cst_Cont', 1., frozen=False, min=0.95, max=1.05)

        self.t1_lam_0 = Parameter(name, 't1_lam_0', 5000., frozen=False, min=0.0)
        self.t1_b = Parameter(name, 't1_b', 3.5, frozen=False, min=1e-12)
        self.t1_d = Parameter(name, 't1_d', 0.0005, frozen=False, min=0)
        self.t1_tau_0 = Parameter(name, 't1_tau_0', 0.1, frozen=False, min=0.0)

        self.t2_lam_0 = Parameter(name, 't2_lam_0', 5000., frozen=False, min=0.0)
        self.t2_b = Parameter(name, 't2_b', 3.5, frozen=False, min=1e-12)
        self.t2_d = Parameter(name, 't2_d', 0.0005, frozen=False, min=0)
        self.t2_tau_0 = Parameter(name, 't2_tau_0', 0.1, frozen=False, min=0.0)

        self.resolution = 80000

        if same_b: self.t1_b = self.t2_b
        if same_d: self.t1_d = self.t2_d

        ArithmeticModel.__init__(self, name,
                                 (self.Cst_Cont, self.t1_lam_0, self.t1_b, self.t1_d, self.t1_tau_0,
                                  self.t2_lam_0, self.t2_b, self.t2_d, self.t2_tau_0))

    def calc(self, pars, x, *args, **kwargs):
        Cst_Cont, t1_lam_0, t1_b, t1_d, t1_tau_0, t2_lam_0, t2_b, t2_d, t2_tau_0  = pars

        transmission1 = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=t1_lam_0, b=t1_b, d=t1_d, tau_0=t1_tau_0)
        transmission2 = eMath.voigtOpticalDepthAbsorption(lam=x, lam_0=t2_lam_0, b=t2_b, d=t2_d, tau_0=t2_tau_0)

        transmission = np.ones_like(x) * Cst_Cont
        transmission = transmission * transmission1 * transmission2



        ## instrumental brodening
        dx = (np.max(x) - np.min(x)) / len(x)
        gsigma = np.mean(x) / self.resolution / 2.35482
        gx = np.arange(-3*gsigma, 3*gsigma, dx)
        gaussian_weight = np.exp(-(gx/gsigma)**2 / 2)
        gaussian_weight = gaussian_weight / np.sum(gaussian_weight)
        m = len(gx)
        n = len(x)

        transmission_out = np.append(np.ones(m-1)*Cst_Cont, transmission)
        transmission_out = np.append(transmission_out, np.ones(m-1)*Cst_Cont)
        transmission_out = np.convolve(transmission_out,gaussian_weight,mode="same")
        transmission_out = transmission_out[m-1:n+m-1]

        ## broadening to be added

        return transmission_out
