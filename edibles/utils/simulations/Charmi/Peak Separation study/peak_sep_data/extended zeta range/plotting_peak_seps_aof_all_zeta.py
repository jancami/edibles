# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 00:20:43 2022

@author: Charmi Bhatt
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(15,6))

# plt.axhline(y = 1.27, color = 'black', linestyle = '-')
# plt.axhline(y = 1.46, color = 'black', linestyle = '-')




B005z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z02.iloc[:,0], B005z02.iloc[:,5])

B005z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z04.iloc[:,0], B005z04.iloc[:,5])

B005z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z06.iloc[:,0], B005z06.iloc[:,5])

B005z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z08.iloc[:,0], B005z08.iloc[:,5])

B005z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z1.iloc[:,0], B005z1.iloc[:,5])

B005z12 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.2_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z12.iloc[:,0], B005z12.iloc[:,5])

B005z14 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.4_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z14.iloc[:,0], B005z14.iloc[:,5])

B005z16 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.6_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z16.iloc[:,0], B005z16.iloc[:,5])

B005z18 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.8_B=0.005_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B005z18.iloc[:,0], B005z18.iloc[:,5], color='black')

good_B005z02 = B005z02[(B005z02.iloc[:,0] <= 65)]
good_B005z04 = B005z04[(B005z04.iloc[:,0] <= 65)]
good_B005z06 = B005z06[(B005z06.iloc[:,0] <= 65)]
good_B005z08 = B005z08[(B005z08.iloc[:,0] <= 65)]
good_B005z1 = B005z1[(B005z1.iloc[:,0] <= 65)]
good_B005z12 = B005z12[(B005z1.iloc[:,0] <= 65)]
good_B005z14 = B005z14[(B005z12.iloc[:,0] <= 65)]
good_B005z16 = B005z16[(B005z14.iloc[:,0] <= 65)]
good_B005z18 = B005z18[(B005z16.iloc[:,0] <= 65)]


list_a = []
list_a.append(good_B005z02.iloc[:,5])
list_a.append(good_B005z04.iloc[:,5])
list_a.append(good_B005z06.iloc[:,5])
list_a.append(good_B005z08.iloc[:,5])
list_a.append(good_B005z1.iloc[:,5])
list_a.append(good_B005z12.iloc[:,5])
list_a.append(good_B005z14.iloc[:,5])
list_a.append(good_B005z16.iloc[:,5])
list_a.append(good_B005z18.iloc[:,5])



#print(list_a)

i_sum_a = []
for i in range(len(list_a)):
    i_sum_a.append(np.sum(list_a[i]/0.005))
sum_of_list_a = np.sum(i_sum_a)
#print(sum_of_list_a)


i_length_a = []
for i in range(len(list_a)):
    i_length_a.append(len(list_a[i]))
length_of_list_a = np.sum(i_length_a)
#print(length_of_list_a)


mean_of_list_a = sum_of_list_a/length_of_list_a
print(mean_of_list_a)


# '================================'

B01z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z02.iloc[:,0], B01z02.iloc[:,5])

B01z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z04.iloc[:,0], B01z04.iloc[:,5])

B01z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z06.iloc[:,0], B01z06.iloc[:,5])

B01z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z08.iloc[:,0], B01z08.iloc[:,5])

B01z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z1.iloc[:,0], B01z1.iloc[:,5])

B01z12 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.2_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z12.iloc[:,0], B01z12.iloc[:,5])

B01z14 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.4_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z14.iloc[:,0], B01z14.iloc[:,5])

B01z16 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.6_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z16.iloc[:,0], B01z16.iloc[:,5])

B01z18 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1.8_B=0.01_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B01z18.iloc[:,0], B01z18.iloc[:,5],color ='black')



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
# list_b.append(good_B01z02.iloc[:,5])
# list_b.append(good_B01z04.iloc[:,5])
# list_b.append(good_B01z06.iloc[:,5])
# list_b.append(good_B01z08.iloc[:,5])
# list_b.append(good_B01z1.iloc[:,5])
# list_b.append(good_B01z12.iloc[:,5])
# list_b.append(good_B01z14.iloc[:,5])
# list_b.append(good_B01z16.iloc[:,5])
# list_b.append(good_B01z18.iloc[:,5])

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
plt.plot(B05z0.iloc[:,0], B05z0.iloc[:,5]/0.05)

B05z02 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.2_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z02.iloc[:,0], B05z02.iloc[:,5]/0.05, color = 'black')

B05z04 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.4_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z04.iloc[:,0], B05z04.iloc[:,5]/0.05)

B05z06 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.6_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z06.iloc[:,0], B05z06.iloc[:,5]/0.05)

B05z08 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-0.8_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z08.iloc[:,0], B05z08.iloc[:,5]/0.05)

B05z1 = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Peak Separation study\peak_sep_data\extended zeta range\peak_seps_zeta=-1_B=0.05_delta_B =0.txt", delim_whitespace=(True), header = None)
plt.plot(B05z1.iloc[:,0], B05z1.iloc[:,5]/0.05)


# good_B05z0 = B05z0[(B05z0.iloc[:,0] <= 10)]
# good_B05z02 = B05z02[(B05z02.iloc[:,0] <= 10)]
# good_B05z04 = B05z04[(B05z04.iloc[:,0] <= 10)]
# good_B05z06 = B05z06[(B05z06.iloc[:,0] <= 10)]
# good_B05z08 = B05z08[(B05z08.iloc[:,0] <= 10)]
# good_B05z1 = B05z1[(B05z1.iloc[:,0] <= 10)]

# list_c = []
# list_c.append(good_B05z0.iloc[:,5])
# list_c.append((good_B05z02.iloc[:,5]))
# list_c.append(good_B05z04.iloc[:,5])
# list_c.append(good_B05z06.iloc[:,5])
# list_c.append(good_B05z08.iloc[:,5])
# list_c.append(good_B05z1.iloc[:,5])
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

print(list_a)