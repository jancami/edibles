##import statements
from edibles.utils.simulations.SimulatedContour import Simulated_Contour as sim
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
def WavelengthToWavenumber(values):
    wavenumbers=1/(values*(10**-8))
    return(wavenumbers)
colormap=plt.cm.cool
plt.rcParams.update({'font.size': 10})
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

#####################
#declare any constants#
Jlimit=100
lambda0=6614
Q_Branch=True
Sightline='ParameterSurvey'

########################
##Declare Innitial Values.
SymmetryType='Spherical'
A_init=20*(10**-3)
B_init=20*(10**-3)
C_init=20*(10**-3)


T_init=15
delta_init=0
Q_scale_init=1
##############################
params_to_vary=['B','T','delta']
vary_1=np.linspace(1.0,2.1,11)

vary_2=np.linspace(0.0,0.011,11)


for param in params_to_vary:
    if param=='delta':
        vary=vary_2
    else:
        vary=vary_1
    colours=[colormap(i) for i in np.linspace(0, 0.9, len(vary))]
    for i in range(len(vary)):
        Q_scale=Q_scale_init

        if param=='B':
            B=B_init*vary[i]
            A=C=B
#            A=A_init
#            C=C_init
            T=T_init
            delta=delta_init
        elif param=='T':
            B=B_init
            A=A_init
            C=C_init
            T=T_init*vary[i]
            delta=delta_init
        elif param=='delta':
            B=B_init
            A=A_init
            C=C_init
            T=T_init
            
            delta=delta_init+vary[i]
        
        build=sim(A=A,B=B,C=C,Delta_A=delta*A,Delta_B=delta*B,Delta_C=delta*C,Trot=T,Jlimit=Jlimit, Target=Sightline,lambda0=lambda0,Q_Branch=Q_Branch,Q_scale=Q_scale)
        x_vals=WavelengthToWavenumber(np.asarray(build[0]))
        y_vals=build[1]
        
        plt.plot(x_vals,y_vals,ls='-',color=colours[i],label=r'T: '+(str(round(T,3)))+' B: '+(str(round(B,3)))+' Delta: '+(str(round(delta,5))))
        
    plt.legend(fontsize='x-small')
    plt.title('Parameter Survey - Changing '+param)

    plt.grid(axis='x', color='0.95')
    plt.savefig("KnobTurns/"+str(SymmetryType)+"/Changing"+param+".pdf")
    plt.close()
