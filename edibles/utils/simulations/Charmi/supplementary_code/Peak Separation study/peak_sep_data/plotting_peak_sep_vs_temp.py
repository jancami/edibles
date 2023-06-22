# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 17:47:09 2022

@author: Charmi Bhatt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


B005z4 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.4_B=0.005_delta_B =0.txt", delim_whitespace=(True), header=None)
B005z4_pr =B005z4.iloc[:,3]
B005z4_pq =B005z4.iloc[:,1]
B005z4_qr =B005z4.iloc[:,2]
B005z4_temp = B005z4.iloc[:,0]

B005z45 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.45_B=0.005_delta_B =0.txt", delim_whitespace=(True), header=None)
B005z45_pr =B005z45.iloc[:,3]
B005z45_pq =B005z45.iloc[:,1]
B005z45_qr =B005z45.iloc[:,2]
B005z45_temp = B005z45.iloc[:,0]

B005z5 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.5_B=0.005_delta_B =0.txt", delim_whitespace=(True), header=None)
B005z5_pr =B005z5.iloc[:,3]
B005z5_pq =B005z5.iloc[:,1]
B005z5_qr =B005z5.iloc[:,2]
B005z5_temp = B005z5.iloc[:,0]


B005z55 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.55_B=0.005_delta_B =0.txt", delim_whitespace=(True), header=None)
B005z55_pr =B005z55.iloc[:,3]
B005z55_pq =B005z55.iloc[:,1]
B005z55_qr =B005z55.iloc[:,2]
B005z55_temp = B005z55.iloc[:,0]

'===================================================='

B01z4 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.4_B=0.01_delta_B =0.txt", delim_whitespace=(True), header=None)
B01z4_pr = B01z4.iloc[:,3]
B01z4_pq = B01z4.iloc[:,1]
B01z4_qr = B01z4.iloc[:,2]
B01z4_temp = B01z4 .iloc[:,0]

B01z45 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.45_B=0.01_delta_B =0.txt", delim_whitespace=(True), header=None)
B01z45_pr = B01z45.iloc[:,3]
B01z45_pq = B01z45.iloc[:,1]
B01z45_qr = B01z45.iloc[:,2]
B01z45_temp = B01z45 .iloc[:,0]

B01z5 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.5_B=0.01_delta_B =0.txt", delim_whitespace=(True), header=None)
B01z5_pr = B01z5.iloc[:,3]
B01z5_pq = B01z5.iloc[:,1]
B01z5_qr = B01z5.iloc[:,2]
B01z5_temp = B01z5 .iloc[:,0]

B01z55 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.55_B=0.01_delta_B =0.txt", delim_whitespace=(True), header=None)
B01z55_pr = B01z55.iloc[:,3]
B01z55_pq = B01z55.iloc[:,1]
B01z55_qr = B01z55.iloc[:,2]
B01z55_temp = B01z55 .iloc[:,0]

"========================================="

B05z4 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.4_B=0.05_delta_B =0.txt", delim_whitespace=(True), header= None)
B05z4_pr = B05z4.iloc[:,3]
B05z4_pq = B05z4.iloc[:,1]
B05z4_qr = B05z4.iloc[:,2]
B05z4_temp = B05z4.iloc[:,0]

B05z45 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.45_B=0.05_delta_B =0.txt", delim_whitespace=(True), header= None)
B05z45_pr = B05z45.iloc[:,3]
B05z45_pq = B05z45.iloc[:,1]
B05z45_qr = B05z45.iloc[:,2]
B05z45_temp = B05z45.iloc[:,0]


B05z5 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.5_B=0.05_delta_B =0.txt", delim_whitespace=(True), header= None)
B05z5_pr = B05z5.iloc[:,3]
B05z5_pq = B05z5.iloc[:,1]
B05z5_qr = B05z5.iloc[:,2]
B05z5_temp = B05z5.iloc[:,0]


B05z55 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\delta_B=0\peak_seps_zeta=-0.55_B=0.05_delta_B =0.txt", delim_whitespace=(True), header=None)
B05z55_pr = B05z55.iloc[:,3]
B05z55_pq = B05z55.iloc[:,1]
B05z55_qr = B05z55.iloc[:,2]
B05z55_temp = B05z55.iloc[:,0]


# plt.plot(B005z4_temp, B005z4_pq, color = 'red')
# plt.plot(B005z4_temp, B005z4_qr)
# plt.show()
# plt.plot(B005z45_temp, B005z45_pq, color = 'red')
# plt.plot(B005z45_temp, B005z45_qr)
# plt.show()
# plt.plot(B005z5_temp, B005z5_pq, color = 'red')
# plt.plot(B005z5_temp, B005z5_qr)
# plt.show()
# plt.plot(B005z55_temp, B005z55_pq, color = 'red')
# plt.plot(B005z55_temp, B005z55_qr)
# plt.show()

plt.plot(B005z4_temp, B005z4_qr - B005z4_pq)
plt.plot(B005z45_temp, B005z45_qr - B005z45_pq, color = 'black')
plt.plot(B005z5_temp, B005z5_qr - B005z5_pq)
plt.plot(B005z55_temp, B005z55_qr - B005z55_pq)


# plt.plot(B01z4_temp, B01z4_pq, color = 'red')
# plt.plot(B01z4_temp, B01z4_qr)
# plt.show()
# plt.plot(B01z45_temp, B01z45_pq, color = 'red')
# plt.plot(B01z45_temp, B01z45_qr)
# plt.show()
# plt.plot(B01z5_temp, B01z5_pq, color = 'red')
# plt.plot(B01z5_temp, B01z5_qr)
# plt.show()
# plt.plot(B01z55_temp, B01z55_pq, color = 'red')
# plt.plot(B01z55_temp, B01z55_qr)
# plt.show()

plt.plot(B01z4_temp, B01z4_qr - B01z4_pq)
plt.plot(B01z45_temp, B01z45_qr - B01z45_pq)
plt.plot(B01z5_temp, B01z5_qr - B01z5_pq, color = 'black')
plt.plot(B01z55_temp, B01z55_qr - B01z55_pq)


# plt.plot(B05z4_temp, B05z4_pq, color = 'red')
# plt.plot(B05z4_temp, B05z4_qr)
# plt.show()
# plt.plot(B05z45_temp, B05z45_pq, color = 'red')
# plt.plot(B05z45_temp, B05z45_qr)
# plt.show()
# plt.plot(B05z5_temp, B05z5_pq, color = 'red')
# plt.plot(B05z5_temp, B05z5_qr)
# plt.show()
# plt.plot(B05z55_temp, B05z55_pq, color = 'red')
# plt.plot(B05z55_temp, B05z55_qr)
# plt.show()

plt.plot(B05z4_temp, B05z4_qr - B05z4_pq, color ='black')
plt.plot(B05z45_temp, B05z45_qr - B05z45_pq, color = 'red')
plt.plot(B05z5_temp, B05z5_qr - B05z5_pq, color ='blue')
plt.plot(B05z55_temp, B05z55_qr - B05z55_pq, color ='green')


plt.xlim(0,50)




'PR vs. temp'
plt.figure(figsize=(15,6))

with sns.color_palette("flare", n_colors=4):
    plt.plot(B005z4_temp, B005z4_pr, marker = 'o', label = 'B = 0.005 $cm^{-1}$, zeta = -0.4')
    plt.plot(B005z45_temp, B005z45_pr, marker = 'o', label = 'B = 0.005 $cm^{-1}$, zeta = -0.45')
    plt.plot(B005z5_temp, B005z5_pr, marker = 'o', label = 'B = 0.005 $cm^{-1}$, zeta = -0.5')
    plt.plot(B005z55_temp, B005z55_pr, marker = 'o', label = 'B = 0.005 $cm^{-1}$, zeta = -0.55')


with sns.color_palette("dark:salmon_r", n_colors=4):
    plt.plot(B01z4_temp, B01z4_pr, marker = 'o', label = 'B = 0.01 $cm^{-1}$, zeta = -0.4')
    plt.plot(B01z45_temp, B01z45_pr, marker = 'o', label = 'B = 0.01 $cm^{-1}$, zeta = -0.45')
    plt.plot(B01z5_temp, B01z5_pr, marker = 'o', label = 'B = 0.01 $cm^{-1}$, zeta = -0.5')
    plt.plot(B01z55_temp, B01z55_pr, marker = 'o', label = 'B = 0.01 $cm^{-1}$, zeta = -0.55')

plt.plot(B05z4_temp, B05z4_pr, color = '#abcae4' , marker = 'o', label = 'B = 0.05 $cm^{-1}$, zeta = -0.4')
plt.plot(B05z45_temp, B05z45_pr, color = '#619bcc' , marker = 'o', label = 'B = 0.05 $cm^{-1}$, zeta = -0.45')
plt.plot(B05z5_temp, B05z5_pr, color = '#316a9a' , marker = 'o', label = 'B = 0.05 $cm^{-1}$, zeta = -0.5')
plt.plot(B05z55_temp, B05z55_pr, color = '#193750', marker = 'o', label = 'B = 0.005 $cm^{-1}$, zeta = -0.55')




plt.axhline(y = 1.46, color = 'black', linestyle = '-')
plt.axhline(y = 1.27, color = 'black', linestyle = '-')
plt.ylim(1,3)

plt.legend()


