#/usr/bin/env python

import numpy as np


# in_path = "calibration/Bases_voronoi/arrivals.dat"
in_path = "calibration/Rect10x10_km/sorted_arrivals.dat"


arq = open(in_path, "r")

G,T,R,P = [int(x) for x in arq.readline().split()]

nb_obs_total = 104
nb_observations = nb_obs_total*np.ones((T,G,R,P))
nb_arrivals = np.zeros((T,G,R,P))
sample = np.frompyfunc(list, 0, 1)(np.empty((T,G,R,P), dtype=object))

for line in arq.readlines():
    g,t,r,p,j,val = [int(x) for x in line.split()]
    sample[t,g,r,p].append(val)
    nb_arrivals[t,g,r,p] += val
arq.close()

lam = np.zeros((T,G,R,P))
duration = 0.5

for t in range(T):
    for g in range(G):
        for r in range(R):
            for p in range(P):
                if nb_observations[t,g,r,p] > 0:
                    lam[t,g,r,p] = (nb_arrivals[t,g,r,p] / nb_observations[t,g,r,p])
                    print(t,g,r,p,"arrivals =",nb_arrivals[t,g,r,p],"obs =",nb_observations[t,g,r,p], "lam =",lam[t,g,r,p])





# out_path = "calibration/Bases_voronoi/empirical_estimation.dat"
out_path = "calibration/Rect10x10_km/sorted_empirical_estimation.dat"
out_arq = open(out_path, "w")
out_arq.write(f"{T} {G} {R} {P} {nb_obs_total}\n")
for t in range(T):
    for g in range(G):
        for r in range(R):
            for p in range(P):
                out_arq.write(f"{t} {g} {r} {p} {lam[t,g,r,p]}\n")

out_arq.close()
