import matplotlib.pyplot as plt
import numpy as np

with open("test.txt", "r") as f:
    data = f.readlines()

data = [float(i) for i in data]

count, bins, ignored = plt.hist(data, 50, normed=True)

mu = 0
rho = 3
sigma = 1/np.sqrt(rho * 4/3*np.pi*2.5**3)

# mu = 2.5
# sigma = 0.5
plt.plot(bins, 
       1 / (sigma * np.sqrt(2 * np.pi)) * 
       np.exp(-(bins - mu)**2 / (2 * sigma**2)), linewidth=2, c='r')

plt.show()