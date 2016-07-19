import os
import fileinput
import numpy as np

def calc_pressure():
    # Read pressure file to get the pressure results
    with open("pressure.out","r") as p:
        pressure_data = p.readlines()

    # Split by space and grab the last part, which is the pressure, and then 
    # convert to float.
    pressure = [float(pressure_entry.split()[-1]) for pressure_entry in pressure_data]  

    quart_press = pressure[3*len(pressure)//4:]
    avg_quart_press = np.mean(quart_press)

    return avg_quart_press

def get_viscosity():
    # File contains only 1 number
    with open("viscosity.out", "r") as v:
        viscosity = v.readline()

    return float(viscosity)

def modify_file(line_num, new_var):
    # Modify the input file
    with fileinput.input("inputs.in", inplace=True) as f:
        for linenum, line in enumerate(f):
            
            # For the line that we want to modify
            if linenum == line_num:
                print(new_var)
            # Every other line is same
            else:
                print(line, end="")    

# os.chdir("C:\\Users\\Qi\\Desktop\\C Codes\\Clone")
os.system("make clean && make")

with open("results.txt", "w") as f:
    f.write("Particle Number, Pressure, Viscosity\n")

# 38 ~ 47
for part_num in range(38, 50, 3):
    for trial in range(3):

        # Modify the 2nd line, replace with part_num
        modify_file(line_num=1, new_var=part_num)

        # Keep looping if we have abnormal execution from assertion fail
        ret = 1 
        while ret != 0:
            # Run the simulation            
            ret = os.system("dpd_sim.exe < inputs.in")

        pressure = calc_pressure()
        viscosity = get_viscosity()

        with open("results.txt", "a") as f: 
            f.write("{0}, {1:.5}, {2}\n".format(part_num, pressure, viscosity))

        print("Particle number = {0}, Trial = {1}".format(part_num, trial))
        print("{0}, {1:.5}, {2}".format(part_num, pressure, viscosity))
