import numpy as np
import astropy.constants as cst
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter
import edibles_Working.edibles_function as eF
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
            transmission = eF.voigtAbsorptionLine(lam=x, lam_0=lam_0, b=b, d=d, N=N, f=f)
        else:
            transmission = eF.voigtAbsorptionLine(lam=x, lam_0=lam_0, b=b, d=d, tau_0=tau_0)

        return transmission


class Sightline(ArithmeticModel):

    def __init__(self, star_name, cont):
        self.star_name = star_name
        self.model = cont
        self.resolution = 100000



    def addLine(self, line_name, lam_0, b=3.5, d=0.0005, tau_0=0.5):
        new_line = VoigtLine(name=line_name)
        new_line.lam_0 = lam_0
        new_line.b = b
        new_line.d = d
        new_line.tau_0 = tau_0
        self.model *= new_line


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

        transmission1 = eF.voigtAbsorptionLine(lam=x, lam_0=t1_lam_0, b=t1_b, d=t1_d, tau_0=t1_tau_0)
        transmission2 = eF.voigtAbsorptionLine(lam=x, lam_0=t2_lam_0, b=t2_b, d=t2_d, tau_0=t2_tau_0)

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
