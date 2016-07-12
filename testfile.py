import numpy as np

with open("shear.out","r") as s:
    shear_data = s.readlines()

    # Split by comma and grab all columns except the last, which is a linebreak
    shears = np.array([shear_entry.split(',')[:-1] for shear_entry in shear_data])
    shears = shears.astype(float)

    
