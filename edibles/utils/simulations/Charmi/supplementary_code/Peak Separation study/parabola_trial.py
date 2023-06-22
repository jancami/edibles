# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 19:05:50 2022

@author: Charmi Bhatt


"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import argrelextrema


smooth_data = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\smooth_1.txt", delim_whitespace = True)

x = smooth_data['Position']
y = np.array(1-0.1*(smooth_data['Intensity']/max(smooth_data['Intensity'])))
plt.plot(x,y)

minima = [argrelextrema(y, np.less)]
minima_ind = minima[0][0]
minima_ind_plus_one = np.array(minima_ind) + 1
minima_ind_minus_one = np.array(minima_ind) - 1

P_parabola_x = (x[minima_ind_minus_one[0]], x[minima_ind[0]], x[minima_ind_plus_one[0]])
Q_parabola_x = (x[minima_ind_minus_one[1]], x[minima_ind[1]], x[minima_ind_plus_one[1]])
R_parabola_x = (x[minima_ind_minus_one[2]], x[minima_ind[2]], x[minima_ind_plus_one[2]])


P_parabola_y = (y[minima_ind_minus_one[0]], y[minima_ind[0]], y[minima_ind_plus_one[0]])
Q_parabola_y = (y[minima_ind_minus_one[1]], y[minima_ind[1]], y[minima_ind_plus_one[1]])
R_parabola_y = (y[minima_ind_minus_one[2]], y[minima_ind[2]], y[minima_ind_plus_one[2]])


   


vertex = []

def vertex_of_parabola(x,y):
    x_1 = x[0]
    x_2 = x[1]
    x_3 = x[2]
    y_1 = y[0]
    y_2 = y[1]
    y_3 = y[2]

    a = y_1/((x_1-x_2)*(x_1-x_3)) + y_2/((x_2-x_1)*(x_2-x_3)) + y_3/((x_3-x_1)*(x_3-x_2))

    b = (-y_1*(x_2+x_3)/((x_1-x_2)*(x_1-x_3))
         -y_2*(x_1+x_3)/((x_2-x_1)*(x_2-x_3))
         -y_3*(x_1+x_2)/((x_3-x_1)*(x_3-x_2)))

    c = (y_1*x_2*x_3/((x_1-x_2)*(x_1-x_3))
        +y_2*x_1*x_3/((x_2-x_1)*(x_2-x_3))
        +y_3*x_1*x_2/((x_3-x_1)*(x_3-x_2)))
    
    vertex.append(-b/(2*a))



vertex_of_parabola(P_parabola_x, P_parabola_y)
vertex_of_parabola(Q_parabola_x, Q_parabola_y)
vertex_of_parabola(R_parabola_x, R_parabola_y)


peak_sep_pq = (vertex[1] - vertex[0]) 
peak_sep_qr = (vertex[2] - vertex[1])
peak_sep_pr = (vertex[2] - vertex[0])

print(peak_sep_pq)
print(peak_sep_qr)
print(peak_sep_pr)



