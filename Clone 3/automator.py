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

    # Round to 4 dp
    return round(avg_quart_press, 4)

def calc_temp():
    # Read temp file to get the temp results
    with open("temp.out","r") as t:
        temp_data = t.readlines()

    # Split by space and grab the 2nd and 3rd parts - translation and 
    # rotational temps respectively
    temp_tuples = [temp_entry.split(",") for temp_entry in temp_data]  
    temp_tuples = [(float(temp_entry[1]), float(temp_entry[2])) for temp_entry in temp_tuples]

    trans_temp, rot_temp = zip(*temp_tuples)

    quart_trans = trans_temp[3*len(trans_temp)//4:]
    avg_quart_trans = np.mean(quart_trans)

    quart_rot = rot_temp[3*len(rot_temp)//4:]
    avg_quart_rot = np.mean(quart_rot)

    # Round to 4 dp
    return round(avg_quart_trans, 4), round(avg_quart_rot, 4)


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

#####################
# MODIFY THESE
#####################
# Variable name - Used when we are printing into the output files
var_name = "Particle Number"
# Number of repeats we want to do
num_trial = 5
# Variable that we are modifying
var_range = [66, 68, 70, 72, 74]

# Amount of water present - can either be constant or varying.
wat_range = [10000 - 64*i for i in var_range]
#wat_range = [7000] * len(var_range)

# Line in inputs that the var resides 
# Particle number = 1; Shear rate = 5
var_line_dict = {"Particle Number": 1, "Shear Rate": 5}
var_line = var_line_dict[var_name]
#####################

with open("results.txt", "w") as f:
    f.write("{}, Viscosity, Pressure, Trans Temp, Rot Temp\n".format(var_name))

for trial in range(num_trial):
    for variable, water in zip(var_range, wat_range):

        # Modify the var_line-th line and replace with the new variable
        modify_file(line_num=0, new_var=water)
        modify_file(line_num=var_line, new_var=variable)

        # Keep looping if we have abnormal execution from assertion fail
        ret = 1 
        while ret != 0:
            # Run the simulation            
            ret = os.system("dpd_sim.exe < inputs.in")

        pressure = calc_pressure()
        viscosity = get_viscosity()
        trans_temp, rot_temp = calc_temp()

        with open("results.txt", "a") as f: 
            f.write("{0}, {1}, {2}, {3}, {4}\n".format(variable, viscosity, pressure, trans_temp, rot_temp))

        print("\n\n/*******************************************************************/")
        print("{} = {}, Trial = {}".format(var_name, variable, trial))
        print("{} = {}, Viscosity = {}, Presure = {}, Trans Temp = {}, Rot Temp = {}".format(var_name, variable, viscosity, pressure, trans_temp, rot_temp))
        print("/*******************************************************************/\n")
