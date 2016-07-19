import numpy as np
import matplotlib.pyplot as plt

##with open("shear.out","r") as s:
##    shear_data = s.readlines()
##
##    # Split by comma and grab all columns except the last, which is a linebreak
##    shears = np.array([shear_entry.split(',')[:-1] for shear_entry in shear_data])
##    shears = shears.astype(float)

with open("water.out", "r") as f:
    water_dat = f.readlines()

water_dat = water_dat[1:]
water_dat = [line.split(",") for line in water_dat]
water_dat = [[float(i) for i in line] for line in water_dat]
water_dat = np.array(water_dat)

y_vx = water_dat[:, [3, 5]]
##y_vx.sort(axis=0)

##x_min = -5
##x_max = 5
##bins = 10
##x_bounds = np.linspace(x_min, x_max, bins + 1)
##grouped = [[] for _ in range(bins)]
##
##list_idx = 0 
##for row in y_vx:
##    while not x_bounds[list_idx] < row[0] <= x_bounds[list_idx + 1]:
##        list_idx += 1
##    assert(list_idx < len(x_bounds))
##    grouped[list_idx].append(row[1])
##
##mean_vel = [sum(row) / len(row) for row in grouped]
##num_water = [len(row) for row in grouped]
##
##plt.subplot(121)
##plt.plot(mean_vel, 'bo')
##plt.subplot(122)
##plt.plot(num_water, 'ro')
##plt.show()

        
