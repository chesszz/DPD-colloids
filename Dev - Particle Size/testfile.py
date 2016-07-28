with open("temp.out","r") as t:
    temp_data = t.readlines()

# Split by space and grab the 2nd and 3rd parts - translation and 
# rotational temps respectively
temp_tuples = [temp_entry.split(",") for temp_entry in temp_data]  
temp_tuples = [(float(temp_entry[1]), float(temp_entry[2])) for temp_entry in temp_tuples]

trans_temp, rot_temp = zip(*temp_tuples)
