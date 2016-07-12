import numpy as np

with open("shear.out","r") as s:
    shear_data = s.readlines()

    # Split by comma and grab all columns except the first, which is the time
    shears = np.array([shear_entry.split(',')[1:] for shear_entry in shear_data])
    shears = shears.astype(float)

    shear_xy = shears[:,0]
    shear_xz = shears[:,1]
    shear_yz = shears[:,2]
    
