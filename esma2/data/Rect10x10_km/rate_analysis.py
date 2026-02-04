import numpy as np
import matplotlib.pyplot as plt

arq = open("sorted_empirical_estimation.dat", "r")
T,G,R,P,nb_obs = [int(x) for x in arq.readline().strip().split()]
lam =  np.zeros((T,G,R,P))
for line in arq.readlines():
    t,g,r,p,rate = line.strip().split()
    lam[int(t),int(g),int(r),int(p)] = float(rate)
arq.close()

lam_tg = np.sum(lam, axis=(2,3))
lam_friday = lam_tg[36:40,4]

lam_crt = np.zeros((P,R,G*T))
for t in range(T):
    for g in range(G):
        for r in range(R):
            for p in range(P):
                lam_crt[p,r,g*T + t] = lam[t,g,r,p]

lam_t = np.sum(lam_crt, axis=(0,1))[228:236]
#plot rates for lam_t and lam_friday together
plt.figure(figsize=(10, 6))
plt.plot(lam_t, marker='o', label='lam_t', linestyle='--')
plt.plot(lam_friday, marker='s', label='lam_friday',linestyle="-.")
plt.title("Comparison of Rates")
plt.xlabel("Time Slot")
plt.ylabel("Rate")
plt.legend()
plt.grid()
plt.show()