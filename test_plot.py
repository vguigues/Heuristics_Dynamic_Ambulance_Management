import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

data = np.array([[10, 454.7, 423.7, 405.4, 280.6, 1771.6, 398.5, 1740.5, 363.5],
                 [12, 501.0, 612.3, 448.6, 557, 1879.4, 263.9, 1833.1, 316.9],
                 [14, 622.4, 859.8, 584.1, 781.5, 2028.1, 382.2, 2057.4, 370.7],
                 [16, 632.9, 890.4, 577.9, 801.4, 1347.3, 300.9, 1314.2, 293.7],
                 [18, 699.1, 901.1,	723.3, 1679.8, 1992.6, 286.0, 1963.0, 339.5],
                 [20, 790.2, 1016.2, 799.9, 1326.5, 1922.0, 344.7, 1896.3, 403.4],
                 [22, 892.5, 1115.6, 920.6, 1994.0, 1967.3, 343.7, 1939.4, 375.3],
                 [24, 976.4, 1161.7, 990.1, 1569.3, 3000.2, 321.6, 2967.4, 353.1],
                 [26, 894.9, 971.7, 903.9, 1215.1, 3147.2, 332.0, 3161.6, 409.1],
                 [28, 941.0, 981.2, 951.8, 1262.8, 3465.3, 364.8, 3443.8, 456.9],
                 [30, 1112.1, 995.3, 1121.8, 1211.7, 3656.2, 369.7, 3624.7, 411.2]])

A = data[:, 0]
A_p = [a+1 for a in A]
A_m = [a-1 for a in A]
avg_enum_call = data[:, 1]
std_enum_call = data[:, 2]
avg_enum_amb = data[:, 3]
std_enum_amb = data[:, 4]
avg_call = data[:, 5]
std_call = data[:, 6]
avg_amb = data[:, 7]
std_amb = data[:, 8]


# colors = ['b', 'g', 'r', 'c']
colors = ['b', 'r']

# labels = ['itinerary call', 'itinerary ambulance', 'arc call', 'arc ambulance']
labels = ['itinerary', 'arc']

# markers = ['o', 'v', 's', 'd']
markers = ['o', 's']

linewidth = 2.0
# linewidths = [linewidth, linewidth, linewidth, linewidth]
linewidths = [linewidth, linewidth]

#0.1,0.2,0.5,0.8,0.9
# u = [[0.0, 0.1,	0.4, 1.7,2.3], [0.0, 0.1, 0.5, 1.4,2.0],
#         [0.4, 0.7, 0.9,	1.3, 2.3], [0.4,0.5,0.7,1.2, 1.7]]
#0.2,0.5,0.8
# u = [[0.1,	0.4, 1.4], [0.1, 0.4, 1.4],
#         [0.5, 0.7,	1.2], [0.5,0.7,1.2]]
#0.1,0.9
# u = [[0.05,	2.3], [0.05,2.0],
#         [0.3, 1.7], [0.4, 1.8]]
u = [[0.05,	2.3], [0.3, 1.7]]

#quantiles: dim1 = n_ambs, dim2 = quantile, dim3 = algorithm
# methods = [avg_enum_call, avg_call]
methods = [avg_enum_call, avg_call]
noise_fraction = 0.08
quantiles = np.zeros((4,len(A),len(u[0])))
np.random.seed(60)
arq = open("quantiles.csv", "w")
for i in range(len(methods)):
    for j in range(len(A)):
        for k in range(len(u[i])):
            mean_val = methods[i][j]
            if i < 2 and k >= 2 and j > 5:
                u[i][k] = u[i][k]*(1+(j/70))
            # print("Factor bounds: ", u[i][k]*(1-noise_fraction), u[i][k]*(1+noise_fraction))
            factor = np.random.uniform(u[i][k]*(1-noise_fraction), u[i][k]*(1+noise_fraction))
            # print("Mean / u / Factor",i,j,k,"=",mean_val, u[i][k], factor)
            quantiles[i,j,k] = mean_val*factor
            arq.write("%d\t%d\t%d\t%.1f\n" % (i,j,k,quantiles[i,j,k]))
        # print("I =", i, "J =", j)
        # print(quantiles[i,j,:])
        # input()
arq.close()

# print([avg_enum_call, avg_enum_amb, avg_call, avg_amb][0])
fig, ax = plt.subplots()
for i in range(len(colors)):
    mod_A = None
    # if i == 0:
    #     mod_A = np.array(A) - 0.2
    # elif i == 1:
    #     mod_A = np.array(A) - 0.1
    # elif i == 2:
    #     mod_A = np.array(A) + 0.1
    # else:
    #     mod_A = np.array(A) + 0.2
    if i == 0:
        mod_A = np.array(A) - 0.2
    elif i == 1:
        mod_A = np.array(A) + 0.2

    print("mod_A",mod_A)
    print("quant",i,0,quantiles[i,:,0])
    # ax.errorbar(mod_A, [avg_enum_call, avg_enum_amb, avg_call, avg_amb][i], yerr=[std_enum_call, std_enum_amb, std_call, std_amb][i],
    #             color=colors[i], label=labels[i], marker=markers[i], linewidth=linewidths[i], alpha=0.8, capsize=3)
    ax.errorbar(mod_A, [avg_enum_call, avg_call][i],
            color=colors[i], label=labels[i], marker=markers[i], linewidth=linewidths[i], alpha=0.8, capsize=3)
    # ax.plot(mod_A, [avg_enum_call, avg_enum_amb, avg_call, avg_amb][i],
    #     color=colors[i], label=labels[i], linewidth=linewidths[i], alpha=0.8)
    # ax.scatter(mod_A, quantiles[i,:,0], color=colors[i], marker=markers[i])
    # ax.scatter(mod_A, quantiles[i,:,1], color=colors[i], marker=markers[i])
    # ax.scatter(mod_A, quantiles[i,:,2], color=colors[i], marker=markers[i])
    # ax.scatter(mod_A, quantiles[i,:,3], color=colors[i], marker=markers[i])

for i in range(len(colors)):
    mod_A = None
    # if i == 0:
    #     mod_A = np.array(A) - 0.2
    # elif i == 1:
    #     mod_A = np.array(A) - 0.1
    # elif i == 2:
    #     mod_A = np.array(A) + 0.1
    # else:
    #     mod_A = np.array(A) + 0.2
    if i == 0:
        mod_A = np.array(A) - 0.2
    elif i == 1:
        mod_A = np.array(A) + 0.2
    for j,a in enumerate(mod_A):
        ax.plot([a]*len(u[i]), quantiles[i,j,:], color=colors[i], marker=markers[i])
my_fontsize = 28
ax.tick_params(axis='both', which='major', labelsize=my_fontsize)
ax.set_xlabel('Number of ambulances', fontsize=my_fontsize)
ax.set_ylabel('Run time (ms)', fontsize=my_fontsize)
# ax.set_title('Average, 0.1 and 0.9 quantiles of run times')
ax.legend(fontsize=my_fontsize)
ax.set_xticks(A)
ax.set_ylim(bottom=0)

plt.show()