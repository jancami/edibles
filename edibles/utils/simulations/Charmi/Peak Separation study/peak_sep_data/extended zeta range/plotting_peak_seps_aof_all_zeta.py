# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 00:20:43 2022

@author: Charmi Bhatt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle

plt.figure(figsize=(15,6))

plt.axhline(y = 1.27, color = 'black', linestyle = '-')
plt.axhline(y = 1.46, color = 'black', linestyle = '-')


B005z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z02.iloc[:,0], B005z02.iloc[:,3], color = '#216a00', linewidth = 2, marker = 'o')
plt.errorbar(B005z02.iloc[:,0], B005z02.iloc[:,3], color = '#216a00', yerr = 0.1)
# B005z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z04.iloc[:,0], B005z04.iloc[:,3], color = '#6fe415', linewidth = 2, marker = 'o')

# B005z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z06.iloc[:,0], B005z06.iloc[:,3], color = '#5fc312', linewidth = 2, marker = 'o')

# B005z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z08.iloc[:,0], B005z08.iloc[:,3], color = '#4fa30f', linewidth = 2, marker = 'o', label = 'B = 0.005 $cm^{-1}$')

# B005z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z1.iloc[:,0], B005z1.iloc[:,3], color = '#4fa30f', linewidth = 2, marker = 'o')

# B005z12 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.2_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z12.iloc[:,0], B005z12.iloc[:,3], color = '#3f820c', linewidth = 2, marker = 'o')

# B005z14 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.4_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z14.iloc[:,0], B005z14.iloc[:,3], color = '#2f6109', linewidth = 2, marker = 'o')

B005z16 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.6_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z16.iloc[:,0], B005z16.iloc[:,3], color = '#216a00', linewidth = 2, marker = 'o', label = 'B = 0.005 $cm^{-1}$')
plt.errorbar(B005z16.iloc[:,0], B005z16.iloc[:,3], color = '#216a00', yerr =0.1)
# B005z18 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.8_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B005z18.iloc[:,0], B005z18.iloc[:,3], color = '#216a00', linewidth = 2, marker = 'o', label = 'B = 0.005 $cm^{-1}$')

# good_B005z02 = B005z02[(B005z02.iloc[:,0] <= 65)]
# good_B005z04 = B005z04[(B005z04.iloc[:,0] <= 65)]
# good_B005z06 = B005z06[(B005z06.iloc[:,0] <= 65)]
# good_B005z08 = B005z08[(B005z08.iloc[:,0] <= 65)]
# good_B005z1 = B005z1[(B005z1.iloc[:,0] <= 65)]
# good_B005z12 = B005z12[(B005z1.iloc[:,0] <= 65)]
# good_B005z14 = B005z14[(B005z12.iloc[:,0] <= 65)]
# good_B005z16 = B005z16[(B005z14.iloc[:,0] <= 65)]
# good_B005z18 = B005z18[(B005z16.iloc[:,0] <= 65)]


# list_a = []
# list_a.append(good_B005z02.iloc[:,3])
# list_a.append(good_B005z04.iloc[:,3])
# list_a.append(good_B005z06.iloc[:,3])
# list_a.append(good_B005z08.iloc[:,3])
# list_a.append(good_B005z1.iloc[:,3])
# list_a.append(good_B005z12.iloc[:,3])
# list_a.append(good_B005z14.iloc[:,3])
# list_a.append(good_B005z16.iloc[:,3])
# list_a.append(good_B005z18.iloc[:,3])



# #print(list_a)

# i_sum_a = []
# for i in range(len(list_a)):
#     i_sum_a.append(np.sum(list_a[i]/0.005))
# sum_of_list_a = np.sum(i_sum_a)
# #print(sum_of_list_a)


# i_length_a = []
# for i in range(len(list_a)):
#     i_length_a.append(len(list_a[i]))
# length_of_list_a = np.sum(i_length_a)
# #print(length_of_list_a)


# mean_of_list_a = sum_of_list_a/length_of_list_a
# print(mean_of_list_a)


# '================================'
    

B01z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z02.iloc[:,0], B01z02.iloc[:,3], color = '#fca300', linewidth = 2, marker = 'o')
plt.errorbar(B01z02.iloc[:,0], B01z02.iloc[:,3], color = '#fca300', yerr = 0.1)
# B01z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z04.iloc[:,0], B01z04.iloc[:,3], color = '#ffcc70', linewidth = 2, marker = 'o')

# B01z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z06.iloc[:,0], B01z06.iloc[:,3], color = '#ffbb3f', linewidth = 2, marker = 'o')

# B01z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z08.iloc[:,0], B01z08.iloc[:,3], color = '#ffaa0e', linewidth = 2, marker = 'o')

# B01z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z1.iloc[:,0], B01z1.iloc[:,3], color = '#dc8e00', linewidth = 2, marker = 'o', label = 'B = 0.01 $cm^{-1}$')

# B01z12 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.2_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z12.iloc[:,0], B01z12.iloc[:,3], color = '#dc8e00', linewidth = 2, marker = 'o')

# B01z14 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.4_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z14.iloc[:,0], B01z14.iloc[:,3], color = '#ab6e00', linewidth = 2, marker = 'o')

B01z16 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.6_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z16.iloc[:,0], B01z16.iloc[:,3], color = '#fca300', linewidth = 2, marker = 'o', label = 'B = 0.01 $cm^{-1}$')
plt.errorbar(B01z16.iloc[:,0], B01z16.iloc[:,3], color = '#fca300', yerr = 0.1)
# B01z18 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.8_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B01z18.iloc[:,0], B01z18.iloc[:,3], color = '#dc8e00', linewidth = 2, marker = 'o', label = 'B = 0.005 $cm^{-1}$')
    


# good_B01z02 = B01z02[(B01z02.iloc[:,0] <= 32.5)]
# good_B01z04 = B01z04[(B01z04.iloc[:,0] <= 32.5)]
# good_B01z06 = B01z06[(B01z06.iloc[:,0] <= 32.5)]
# good_B01z08 = B01z08[(B01z08.iloc[:,0] <= 32.5)]
# good_B01z1 = B01z1[(B01z1.iloc[:,0] <= 32.5)]
# good_B01z12 = B01z12[(B01z12.iloc[:,0] <= 32.5)]
# good_B01z14 = B01z14[(B01z14.iloc[:,0] <= 32.5)]
# good_B01z16 = B01z16[(B01z16.iloc[:,0] <= 32.5)]
# good_B01z18 = B01z18[(B01z18.iloc[:,0] <= 32.5)]



# list_b = []
# list_b.append(good_B01z02.iloc[:,3])
# list_b.append(good_B01z04.iloc[:,3])
# list_b.append(good_B01z06.iloc[:,3])
# list_b.append(good_B01z08.iloc[:,3])
# list_b.append(good_B01z1.iloc[:,3])
# list_b.append(good_B01z12.iloc[:,3])
# list_b.append(good_B01z14.iloc[:,3])
# list_b.append(good_B01z16.iloc[:,3])
# list_b.append(good_B01z18.iloc[:,3])

# i_sum_b = []
# for i in range(len(list_b)):
#     i_sum_b.append(np.sum(list_b[i]/0.01))
# sum_of_list_b = np.sum(i_sum_b)
# #print(sum_of_list_b)


# i_length_b = []
# for i in range(len(list_b)):
#     i_length_b.append(len(list_b[i]))
# length_of_list_b = np.sum(i_length_b)
# #print(length_of_list_b)


# mean_of_list_b = sum_of_list_b/length_of_list_b 
# print(mean_of_list_b)


# '===================================='


B05z0 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=0_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z0.iloc[:,0], B05z0.iloc[:,3], color = '#255075', linewidth = 2, marker = 'o')
plt.errorbar(B05z0.iloc[:,0], B05z0.iloc[:,3], color = '#255075', yerr = 0.1)
# B05z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B05z02.iloc[:,0], B05z02.iloc[:,3], color = '#3d84bf', linewidth = 2, marker = 'o')

# B05z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B05z04.iloc[:,0], B05z04.iloc[:,3], color = '#316a9a', linewidth = 2, marker = 'o', label = 'B = 0.05 $cm^{-1}$')

# B05z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B05z06.iloc[:,0], B05z06.iloc[:,3], color = '#255075', linewidth = 2, marker = 'o')

# B05z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
# plt.plot(B05z08.iloc[:,0], B05z08.iloc[:,3], color = '#193750', linewidth = 2, marker = 'o')

B05z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z1.iloc[:,0][0:2], B05z1.iloc[:,3][0:2], color = '#255075', linewidth = 2, marker = 'o', label = 'B = 0.05 $cm^{-1}$')
plt.errorbar(B05z1.iloc[:,0][0:2], B05z1.iloc[:,3][0:2], color = '#255075', yerr = 0.1)

plt.ylim(1,2)

plt.xlabel('Rotational temperature $T_{rot}$ ( in K)', labelpad=(10))
plt.ylabel('Peak Separation $\Delta \\nu_{PR}$ (in $cm^{-1}$)', labelpad=(10))
plt.legend()
plt.title('                                             Peak Separation $\Delta \\nu_{PR}$ as a function of Rotational temperature $T_{rot}$                                   ', pad = 8, backgroundcolor = 'white', size = 'x-large')


plt.ylim(1,2)
plt.gca().add_patch(Rectangle((2.7,1.27), 6.9, 0.19, linewidth=3, color = '#96d2f1', alpha = 0.7))
plt.gca().add_patch(Rectangle((9.6,1.27), 22.9, 0.19, linewidth=1, color = '#ffd890', alpha = 0.9))
plt.gca().add_patch(Rectangle((19.5,1.27), 45.5, 0.19, linewidth=1, color = '#91ff5f', alpha = 0.6))


#rect2 = matplotlib.patches.Rectangle((2.7,1.27), 9.2, 1.46, color='red')
#rect = plt.Rectangle((2.7,1.27), 9.2, 1.46, color='red')

#ax.add_patch(rect2)
# good_B05z0 = B05z0[(B05z0.iloc[:,0] <= 10)]
# good_B05z02 = B05z02[(B05z02.iloc[:,0] <= 10)]
# good_B05z04 = B05z04[(B05z04.iloc[:,0] <= 10)]
# good_B05z06 = B05z06[(B05z06.iloc[:,0] <= 10)]
# good_B05z08 = B05z08[(B05z08.iloc[:,0] <= 10)]
# good_B05z1 = B05z1[(B05z1.iloc[:,0] <= 10)]

# list_c = []
# list_c.append(good_B05z0.iloc[:,3])
# list_c.append((good_B05z02.iloc[:,3]))
# list_c.append(good_B05z04.iloc[:,3])
# list_c.append(good_B05z06.iloc[:,3])
# list_c.append(good_B05z08.iloc[:,3])
# list_c.append(good_B05z1.iloc[:,3])
# list_c = (np.array(list_c))/0.05
# mean_list_c = np.sum(list_c)/((list_c).shape[0]*(list_c).shape[1])
# # print((list_c).shape)
# # print((list_c))
# print('mean_list_c = ' , mean_list_c)

# print('-------------')
# print(sum_of_list_a +sum_of_list_b + np.sum(list_c))
# print()
# Total_mean = (sum_of_list_a +sum_of_list_b + np.sum(list_c))/(length_of_list_a + length_of_list_b + ((list_c).shape[0]*(list_c).shape[1]))
# print(Total_mean)

