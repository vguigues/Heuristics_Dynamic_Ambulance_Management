#/usr/bin/env python3
#-*-encoding: utf-8-*-

import sys
import pandas as pd
import numpy as np
from scipy.stats import beta
import math
from random import randint

ambs_list = [10, 12,14,16,18, 20, 22, 24, 26, 28, 30]
heuristic_list = ["queue", "prep2", "preparedness", "ordered", "district", "coverage"]


def get_diffs(data_h1, data_h2):
    actual_diffs = []
    penalized_diffs = []
    if len(data_h1[0]) != len(data_h2[0]):
        raise ValueError("Number of actual response times are different between heuristics") 
    
    if len(data_h1[1]) != len(data_h1[1]):
        raise ValueError("Number of penalized response times are different between heuristics") 
    
    if len(data_h1[0]) != len(data_h1[1]) or len(data_h2[0]) != len(data_h2[1]):
        raise ValueError("Number of actual and penalized response times arre different")
    

    for i in range(len(data_h1[0])):
        actual_diffs.append(data_h1[0][i] - data_h2[0][i])
        penalized_diffs.append(data_h1[1][i] - data_h2[1][i])
    
    return (actual_diffs, penalized_diffs)


def read_one_stage(heuristic, n_ambs, n_scen):
    file_arq = "results/one_stage_%s_a%d_s%d.txt" % (heuristic, n_ambs, n_scen)
    result_file = open(file_arq, "r")
    actual_values = []
    penalized_values = []

    for s in range(n_scen):
        num_calls = int(result_file.readline())
        for c in range(num_calls):
            actual, pen = [float(x) for x in result_file.readline().split()]
            actual_values.append(actual)
            penalized_values.append(pen)
    
    result_file.close()
    file_arq = "results/one_stage/one_stage_%s_a%d_s%d.txt" % (heuristic, n_ambs, n_scen)

    return (actual_values, penalized_values)


def read_values(heuristic, n_ambs, H = 2):
    file_arq = ""
    if H == 2:
        file_arq = "results/results_before_T_test/results_random_times/two_stage_%s_a%d_s5.txt" % ("heuristic", n_ambs)
    else:
        file_arq = "results/two_stage_%s_a%d_s5_H%d.txt" % (heuristic, n_ambs, H)
    result_file = open(file_arq, "r")
    actual_values = []
    penalized_values = []
    N = result_file.readline()
    for line in result_file.readlines():
        actual_value, penalized_value = [float(x) for  x in line.split()]
        actual_values.append(actual_value)
        penalized_values.append(penalized_value)
    result_file.close()
    # result_file2 = open("results/results_random_times/two_stage_%s_a%d_s5.txt" % (heuristic, n_ambs), "r")
    # N = result_file2.readline()
    # for line in result_file2.readlines():
    #     actual_value, penalized_value = [float(x) for  x in line.split()]
    #     actual_values.append(actual_value)
    #     penalized_values.append(penalized_value)
    # result_file2.close()

    return (actual_values, penalized_values)

def compare_diffs(heuristic1):
    print("Comparison %s" % (heuristic1))
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        for heuristic2 in heuristic_list:
        # for heuristic2 in ["cg"]:
            if heuristic2 == heuristic1:
                continue
            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            print("%.1f" % (pen_diff.std() / np.sqrt(len(pen_diff))), end="\t")
        print()
    print("Actual")
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        for heuristic2 in heuristic_list:
        # for heuristic2 in ["cg"]:
            if heuristic2 == heuristic1:
                continue
            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            print("%.1f" % (actual_diff.std() / np.sqrt(len(actual_diff))), end="\t")
        print();

def compare_means(heuristic1):
    print("Comparison %s" % (heuristic1))
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        for heuristic2 in heuristic_list:
        # for heuristic2 in ["cg"]:
            if heuristic2 == heuristic1:
                continue
            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            print("%.1f" % (pen_diff.mean()), end="\t")
        print()
    print("Actual")
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        for heuristic2 in heuristic_list:
            if heuristic2 == heuristic1:
                continue

            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            print("%.1f" % (actual_diff.mean()), end="\t")
        print();


def num_stds(heuristic1):
    print("Comparison %s penalized" % (heuristic1))
    header = [x for x in heuristic_list if x != heuristic1]
    print("RT\t%s" % ("\t".join(header)))
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        print("%d\t%.0f" % (n_ambs,pd.Series(data1[1]).mean()), end="\t")
        for heuristic2 in heuristic_list:
        # for heuristic2 in ["cg"]:
            if heuristic2 == heuristic1:
                continue
            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            pen_mean = pen_diff.mean()
            avg_std = pen_diff.std() / np.sqrt(len(pen_diff))
            n_stds = pen_mean / avg_std
            print("%.1f" % (n_stds), end="\t")
        print()

    print("Comparison %s actual" % (heuristic1))
    for n_ambs in ambs_list:
        data1 = read_values(heuristic1, n_ambs)
        print("%d\t%.0f" % (n_ambs,pd.Series(data1[0]).mean()), end="\t")
        for heuristic2 in heuristic_list:
            if heuristic2 == heuristic1:
                continue
            data2 = read_values(heuristic2, n_ambs)
            actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
            actual_mean = actual_diff.mean()
            avg_std = actual_diff.std() / np.sqrt(len(actual_diff))
            n_stds = actual_mean / avg_std
            print("%.1f" % (n_stds), end="\t")
        print();

def print_response_times():
    print("Response Times Penalized")
    header = [x for x in heuristic_list]
    print("%s" % ("\t".join(header)))
    k = 1
    for n_ambs in ambs_list:
        print("%d" % (n_ambs), end="\t")
        for heuristic in heuristic_list:
            data = read_values(heuristic, n_ambs)
            vals = data[k]
            print("%.0f" % (pd.Series(vals).mean()), end="\t")
        print()
    print("Response Times actual")
    k = 0
    for n_ambs in ambs_list:
        print("%d" % (n_ambs), end="\t")
        for heuristic in heuristic_list:
            data = read_values(heuristic, n_ambs)
            vals = data[k]
            print("%.0f" % (pd.Series(vals).mean()), end="\t")
        print()


def num_stds_h(heuristic1):
    print("penalized")
    n_ambs = 20
    data1 = read_values(heuristic1, n_ambs, 3)
    print("%d\t%.0f" % (n_ambs,pd.Series(data1[1]).mean()), end="\t")
    data2 = read_values(heuristic1, n_ambs, 4)
    actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
    pen_mean = pen_diff.mean()
    avg_std = pen_diff.std() / np.sqrt(len(pen_diff))
    n_stds = pen_mean / avg_std
    print("%.1f" % (n_stds), end="\t")
    print()
    print("actual")
    n_ambs = 20
    data1 = read_values(heuristic1, n_ambs, 3)
    print("%d\t%.0f" % (n_ambs,pd.Series(data1[0]).mean()), end="\t")
    data2 = read_values(heuristic1, n_ambs, 4)
    actual_diff, pen_diff = [ pd.Series(x) for x in get_diffs(data1, data2)]
    actual_mean = actual_diff.mean()
    avg_std = actual_diff.std() / np.sqrt(len(actual_diff))
    n_stds = actual_mean / avg_std
    print("%.1f" % (n_stds), end="\t")
    print()


def main():
    # if len(sys.argv) < 2:
    #     print("Usage: compare [heuristic1] [heuristic2]")
    #     return
    # algo = "enumerate"
    # # actual, pen = read_values(algo, 20, 2)

    # # print("H2 Actual mean:", pd.Series(actual).mean())
    # # print("H2 Penalized mean:", pd.Series(pen).mean())
    # actual, pen = read_values(algo, 20, 3)

    # print("H3 Actual mean:", pd.Series(actual).mean())
    # print("H3 Penalized mean:", pd.Series(pen).mean())

    # actual, pen = read_values(algo, 20, 4)

    # print("H3 Actual mean:", pd.Series(actual).mean())
    # print("H3 Penalized mean:", pd.Series(pen).mean())

    # input()

    # num_stds_h(algo)



    # print_response_times()
    # heuristic1 = sys.argv[1]
    # num_stds(heuristic1)
    # compare_diffs(heuristic1)
    # compare_means(heuristic1)
    print("penalized")
    print("A\t%s\n" % ("\t".join(heuristic_list)))
    for n_ambs in ambs_list:
        print("%d" % (n_ambs), end="\t")
        for heuristic in heuristic_list:
            actual, pen = read_one_stage(heuristic, n_ambs, 15)
            print("%.0f" % (np.mean(actual)), end="\t")
        print()

    print("actual")
    print("A\t%s\n" % ("\t".join(heuristic_list)))
    for n_ambs in ambs_list:
        print("%d" % (n_ambs), end="\t")
        for heuristic in heuristic_list:
            actual, pen = read_one_stage(heuristic, n_ambs, 15)
            print("%.0f" % (np.std(actual)), end="\t")
        print()





if __name__ == "__main__":
    main()