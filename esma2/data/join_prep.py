import numpy as np
prep1_file = open("preparedness_us_16a_16b.txt", "r")

prep2_file = open("preparedness_us_16a_18b.txt", "r")
prep_file = open("preparedness_us_16a_34b.txt", "w")

num_bases, num_times, n_types_amb, n_supply1, n_supply2, n_sub_regions = [int(x) for x in prep1_file.readline().split()]
prep_file.write(f"34 {num_times} {n_types_amb} {n_supply1} {n_supply2}\n")
prep_array = np.zeros((34, num_times, n_supply1, n_supply2))
for line in prep1_file.readlines():
    tokens = line.split()
    b = int(tokens[1])
    t = int(tokens[2])
    supply0 = int(tokens[3])
    supply1 = int(tokens[4])
    val = float(tokens[5])
    prep_array[b,t,supply0,supply1] = val
    # prep_file.write(f"{b} {t} {supply0} {supply1} {val}\n")
prep1_file.close()
bases_offset = 16
prep2_file.readline()
for line in prep2_file.readlines():
    tokens = line.split()
    b = int(tokens[1])
    t = int(tokens[2])
    supply0 = int(tokens[3])
    supply1 = int(tokens[4])
    val = float(tokens[5])
    
    prep_array[bases_offset+ b,t,supply0,supply1] = val
    # prep_file.write(f"{bases_offset + b} {t} {supply0} {supply1} {val}\n")
prep2_file.close()

for b in range(34):
    for t in range(num_times):
        for s0 in range(n_supply1):
            for s1 in range(n_supply2):
                prep_file.write(f"{b} {t} {s0} {s1} {prep_array[b,t,s0,s1]}\n")

prep_file.close()