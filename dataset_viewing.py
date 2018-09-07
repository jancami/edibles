import matplotlib.pyplot as plt
import pickle

with open('datasets.dict', 'rb') as handle:
    dictionary = pickle.load(handle)

plt.plot(dictionary['dataset1'][0], dictionary['dataset1'][1])
plt.suptitle('Dataset1')
plt.show()


plt.plot(dictionary['dataset2'][0], dictionary['dataset2'][1])
plt.suptitle('Dataset2')
plt.show()


plt.plot(dictionary['dataset3'][0], dictionary['dataset3'][1])
plt.suptitle('Dataset3')
plt.show()


plt.plot(dictionary['dataset4'][0], dictionary['dataset4'][1])
plt.suptitle('Dataset4')
plt.show()


plt.plot(dictionary['dataset5'][0], dictionary['dataset5'][1])
plt.suptitle('Dataset5')
plt.show()


plt.plot(dictionary['dataset6'][0], dictionary['dataset6'][1])
plt.suptitle('Dataset6')
plt.show()


plt.plot(dictionary['dataset7'][0], dictionary['dataset7'][1])
plt.suptitle('Dataset7')
plt.show()


plt.plot(dictionary['dataset8'][0], dictionary['dataset8'][1])
plt.suptitle('Dataset8')
plt.show()
