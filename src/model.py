import numpy as np
import astropy.constants as cst
from sherpa import models
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter
from sherpa.data import Data1D
from sherpa.instrument import Kernel, ConvolutionKernel, ConvolutionModel
import sherpa.ui as UI
import edibles.src.math as eMath
import edibles.src.datahandling as DataHandling
import math

class VoigtLine(ArithmeticModel):

    def __init__(self, name='VoigtLine'):

        self.lam_0 = Parameter(name, 'lam_0', 5000., frozen=False, min=0.0)
        self.b = Parameter(name, 'b', 1.5, frozen=False, min=1e-12)
        self.d = Parameter(name, 'd', 0.0005, frozen=False, min=0)
        self.N = Parameter(name, 'N', 999, frozen=True, hidden=True, min=0.0)
        self.f = Parameter(name, 'f', 999, frozen=True, hidden=True, min=0.0)
        self.tau_0 = Parameter(name, 'tau_0', 0.1, frozen=False, min=0.0)

        ArithmeticModel.__init__(self, name,
                                 (self.lam_0, self.b, self.d,self.N, self.f, self.tau_0))


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
        self.b = Parameter(name, 'b', 1.5, frozen=False, min=1e-12)
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

class Cloud(ArithmeticModel):

    def __init__(self, name="cloud", velocity = 0.0):
        self.name = name
        self.lines = []
        self.velocity = velocity
        self.instrumental = None
        #self.model = None
        #self.compiled = False

    def set_velocity(self, velocity):
        self.velocity = velocity
        if len(self.lines) > 0:
            self.lines[0].v_offset = velocity

    def addLines(self, lam_0, line_name=None, b=1.5, d = 0.0005, tau_0=0.1):
        lam_0 = DataHandling.parseInput(1,lam_0,checklen=False)
        line_name, b, d, tau_0 = DataHandling.parseInput(len(lam_0), line_name, b, d, tau_0)
        for i, name in enumerate(line_name):
            if name is None:
                name = self.name + "_known_wav_line"+str(len(self.lines))
            new_line = VoigtLine_KnownWav(name=name)
            new_line.lam_0 = lam_0[i]
            new_line.b = b[i]
            new_line.d = d[i]
            new_line.tau_0 = tau_0[i]

            if  len(self.lines)>0:
                UI.link(new_line.v_offset, self.lines[0].v_offset)

            self.lines.append(new_line)

    def importInstrumental(self, kernel):
        x_g = np.ones_like(kernel)
        kernel = Data1D("kernel", x_g, kernel)
        self.instrumental = ConvolutionKernel(kernel, name="Conv")

    def compileModel(self, *kwords, add_instrumental = True, sightline=False, conv_correction=None):

        if conv_correction is not None and (sightline or self.instrumental is not None):
            self.set_velocity(self.velocity - conv_correction)

        self.lines[0].v_offset = self.velocity

        link_b = False
        freeze_d = False
        if "link_b" in kwords:
            link_b = True
        if "freeze_d" in kwords:
            freeze_d = True

        model = self.lines[0]
        for i, line in enumerate(self.lines[1:]):
            if link_b:
                UI.link(line.b, self.lines[0].b)
            else:
                line.b.link = None

            if freeze_d:
                line.d.frozen = True
            else:
                line.d.frozen = False
            model = model * line

        if not sightline:
            model *= models.Const1D()

        if add_instrumental and self.instrumental is not None:
            model = self.instrumental(model)

        return model

    def __str__(self):
        str_out = "="*15 + " Cloud: " + self.name + " " + "="*15 + "\n"
        str_out += "Velocity: {v:.2f} km/s\n".format(v = self.velocity)
        str_out += "{n} lines in this cloud\n".format(n = len(self.lines))
        str_out += "="*20 + "\n"
        for i,line in enumerate(self.lines):
            str_out += "Line {idx}: {name} at {wav}\n".format(idx=i, name=line.name, wav=line.lam_0.val)

        return str_out + "\n\n"

class Sightline(ArithmeticModel):

    def __init__(self, name="Model"):
        self.name = name
        self.clouds = []
        self.cloud_names = []
        self.clouds_velocities = []
        self.lines = []
        self.instrumental = None

    def importCloud(self, new_cloud):
        self.clouds.append(new_cloud)
        self.clouds_velocities.append(new_cloud.velocity)
        if new_cloud.name.lower() == "cloud":
            self.cloud_names.append("Cloud_{i}".format(i=len(self.clouds)))
        else:
            self.cloud_names.append(new_cloud.name)

    def addLine_tolines(self, lam_0, line_name=None, b=1.5, d=0.0005, tau_0=0.1):
        # create a Voigtline instance and append to self.lines, i.e. lines defined by wavelength rather than velocity
        lam_0 = DataHandling.parseInput(1, lam_0, checklen=False)
        line_name, b, d, tau_0 = DataHandling.parseInput(len(lam_0), line_name, b, d, tau_0)
        for i, name in enumerate(line_name):
            if name is None:
                name = "Line_{idx}".format(idx = len(self.lines))
            new_line = VoigtLine(name=name)
            new_line.lam_0 = lam_0[i]
            new_line.b = b[i]
            new_line.d = d[i]
            new_line.tau_0 = tau_0[i]
            self.lines.append(new_line)

    def addLine_toclouds(self, lam_0, clouds=None, cloud_names=None, velocity=None, line_name=None, b=1.5, d=0.0005, tau_0=0.1):
        assert clouds is not None or velocity is not None or cloud_names is not None,\
            "Please specific which cloud"
        # velocity override cloud_name override cloud(idx)
        if clouds is not None:
            clouds = DataHandling.parseInput(1,clouds,checklen=False)
            clouds_existing = []
            for cloud in clouds:
                if cloud < len(self.clouds):
                    clouds_existing.append(cloud)
                else:
                    print("Cloud{idx} not found".format(idx=cloud))

        if cloud_names is not None:
            cloud_names = DataHandling.parseInput(1, cloud_names, checklen=False)
            clouds_existing = []
            for cloud_name in cloud_names:
                if cloud_name in self.cloud_names:
                    clouds_existing.append(self.cloud_names.index(cloud_name))
                else:
                    print("Cloud " + cloud_name + " not found")

        if velocity is not None:
            velocity = DataHandling.parseInput(1, velocity, checklen=False)
            clouds_existing, clouds_new = [], []
            for v in velocity:
                if v in self.clouds_velocities:
                    clouds_existing.append(self.clouds_velocities.index(v))
                else:
                    clouds_new.append(v)

        if len(clouds_existing) > 0:
            for cloud in clouds_existing:
                self.__addLine_existingcloud__(lam_0, cloud=cloud, line_name=line_name, b=b, d=d, tau_0=tau_0)

        if len(clouds_new) > 0:
            for v in clouds_new:
                self.__addLine_newcloud__(lam_0, velocity=v, line_name=line_name, b=b, d=d, tau_0=tau_0)
                print("New cloud at {v:.2f} km/s created".format(v=v))

    def __addLine_existingcloud__(self, lam_0, cloud=None, line_name=None, b=1.5, d=0.0005, tau_0=0.1):
        self.clouds[cloud].addLines(lam_0, line_name=line_name, b=b, d=d, tau_0=tau_0)

    def __addLine_newcloud__(self, lam_0, velocity=None, line_name=None, b=1.5, d=0.0005, tau_0=0.1):
        new_cloud_name = "Cloud{i}".format(i = len(self.clouds))
        new_cloud = Cloud(name=new_cloud_name, velocity=velocity)
        new_cloud.addLines(lam_0,line_name=line_name, b=b, d=d, tau_0=tau_0)
        self.importCloud(new_cloud)

    def importInstrumental(self, kernel):
        x_g = np.ones_like(kernel)
        kernel = Data1D("kernel", x_g, kernel)
        self.instrumental = ConvolutionKernel(kernel, name="conv")

    def compileModel(self, *kwords, add_instrumental = True, conv_correction=None):

        if conv_correction is not None and self.instrumental is None:
            conv_correction = None

        CstCont = models.Const1D()
        model = CstCont
        for line in self.lines:
            model *= line

        for cloud in self.clouds:
            model *= cloud.compileModel(*kwords, sightline=True, add_instrumental=False,
                                        conv_correction=conv_correction)

        if add_instrumental and self.instrumental is not None:
            model = self.instrumental(model)

        return model

    def __str__(self):
        str_out = "="*15 + " Cloud: " + self.name + " " + "="*15 + "\n"
        str_out += "{n_line} lines and {n_cloud} clouds\n\n".format(n_line=len(self.lines), n_cloud=len(self.clouds))
        str_out += "="*20 + "\n\n"

        str_out += "=" * 15 + " Lines " + "="*15 + "\n"
        for i, line in enumerate(self.lines):
            str_out += "Line {idx}: {name} at {wav}\n".format(idx=i, name=line.name, wav=line.lam_0.val)

        str_out += "\n"

        for i, cloud in enumerate(self.clouds):
            str_out += cloud.__str__()

        return str_out + "\n\n"


class Sightline_old(ArithmeticModel):

    def __init__(self, star_name, cont):
        self.star_name = star_name
        self.model = cont
        self.resolution = 100000
        self.component=[cont]

    def addLine(self, lam_0, *kwords, line_name="Line", b=1.5, d=0.0005, tau_0=0.1):
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

    def addLineSeries(self, lam_0, *kwords, line_name="Line", b=1.5, d=0.0005, tau_0=0.1):
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

    def addInstrumental(self, kernel, name="Kernel"):
        kernel = np.array(kernel)
        if np.sum(kernel) != 1:
            kernel = kernel / np.sum(kernel)
        kernel = Data1D("kernel", kernel, kernel)
        conv = ConvolutionKernel(kernel, name="conv")
        self.model = conv(self.model)






class DualTelluric(ArithmeticModel):

    def __init__(self, name='DuleTelluric', same_b=False, same_d=False):
        self.Cst_Cont = Parameter(name, 'Cst_Cont', 1., frozen=False, min=0.95, max=1.05)

        self.t1_lam_0 = Parameter(name, 't1_lam_0', 5000., frozen=False, min=0.0)
        self.t1_b = Parameter(name, 't1_b', 1.5, frozen=False, min=1e-12)
        self.t1_d = Parameter(name, 't1_d', 0.0005, frozen=False, min=0)
        self.t1_tau_0 = Parameter(name, 't1_tau_0', 0.1, frozen=False, min=0.0)

        self.t2_lam_0 = Parameter(name, 't2_lam_0', 5000., frozen=False, min=0.0)
        self.t2_b = Parameter(name, 't2_b', 1.5, frozen=False, min=1e-12)
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
