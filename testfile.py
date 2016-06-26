##import numpy as np
##
##with open("outputs.out","r") as f:
##    data = f.readlines()
##
### Chop off first line which is the column header
##data = data[1:]
##
### Split each line by commas, and then convert each entry from str->float
##for index, line in enumerate(data):
##    data[index] = line.split(",")
##    for inner_index, entry in enumerate(data[index]):
##        data[index][inner_index] = float(entry)
##
##data = np.array(data)
##print(data[0:3,6])
##print(data[0:3,7])

##with open("inputs.in", "r") as g:
##    inputs = g.readlines()
##
### Only have 11 inputs
##input_data = inputs[:11]
##input_data = list(map(float, input_data))
##
### Skip 1 line because it's a space
##labels = inputs[12:]
### Chop out the part before the first space
##for index, string in enumerate(labels):
##    string = string.strip()
##    bef,space,aft = string.partition(" ")
##    labels[index] = aft
##
##input_dict = dict(zip(labels,input_data))

with open("temp.out","r") as t:
    temp_data = t.readlines()

temp = [float(temp_entry.split()[-1]) for temp_entry in temp_data]
print(temp)


