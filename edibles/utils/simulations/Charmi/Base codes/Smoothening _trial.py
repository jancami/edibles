import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

linelist = pd.read_csv(r"C:\Users\Charmi Bhatt\OneDrive\Desktop\my_local_github\edibles\edibles\utils\simulations\Charmi\Kerr's conditions\condition_c\pgopher_kerr_condition_c_Jmax_300_A1g_E1u.txt", delim_whitespace=(True))

print(linelist['strength'])

waveno_stepsize = 0.075
grid_size = int(((np.max(linelist['position']) + 0.5) - (np.min(linelist['position']) - 0.5))/waveno_stepsize)  

smooth_wavenos = np.linspace(np.min(linelist['position']) - 0.5 , np.max(linelist['position']) + 0.5, grid_size)
smooth_intensities = np.zeros(smooth_wavenos.shape)

resolution = 1e5

# for idx,wavepoint in np.ndenumerate(smooth_wavenos):
#     w_int = ss.norm.pdf(linelist['position'], wavepoint, wavepoint/(2.355*resolution)) * (linelist['strength']) 
    
#     smooth_intensities[idx] = w_int.sum()

   

plt.figure(figsize=(30,6))
plt.stem(linelist['position'], 1-0.1*(linelist['strength']/max(linelist['strength'])),  label='calculated', linefmt='y', markerfmt='yo', bottom=1)
#plt.plot(smooth_wavenos, 1-0.1*(smooth_intensities/max(smooth_intensities)), color = 'black', linewidth = 3)
#plt.title('Calculated: T = ' + str(T) + 'K ,  ground_B =  ' + str(ground_B))
# plt.xlim(15119.0, 15119.7)
plt.show()
