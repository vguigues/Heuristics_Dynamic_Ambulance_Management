#/usr/bin/env python3
#-*-encoding: utf-8 -*-


import sys
import pandas as pd
import numpy as np
from scipy.stats import beta
import math


def print_quantiles(n_ambs):
    file_name = "results/run_times_%d.txt" % (n_ambs)
    try:
        run_times = pd.read_csv(file_name, sep=" ", skiprows=1, index_col=False)
    except FileNotFoundError:
        return "%d" % (n_ambs) + "\t-\t-"*4, [math.inf, math.inf, math.inf, math.inf]
    
    run_times.columns = ["num_calls", "index_solver", "run_time"]
        
    enum_call = run_times[run_times["index_solver"] == 0]
    enum_amb = run_times[run_times["index_solver"] == 4]
    run_times_call = run_times[run_times["index_solver"] == 1]
    run_times_amb = run_times[run_times["index_solver"] == 2]

    # print("Enum call:")
    # print(enum_call.describe())
    # input()
    # print("Enum amb:")
    # print(enum_amb.describe())
    # input()
    # print("Arc call:")
    # print(run_times_call.describe())
    # input()
    # print("Arc amb:")
    # print(run_times_amb.describe())
    # input()
    
    Q = [0,0.2,0.5,0.8,0.9,1]
    print("n_ambs = %d\tmean\t0\t   0.2\t0.5\t0.8\t0.9\t1" % n_ambs)
    print("enum_calls\t%.1f" % (enum_call["run_time"].mean()), end="\t")
    mean_val = enum_call["run_time"].mean()
    for q in Q:
        print("%.1f" % (enum_call["run_time"].quantile(q) / mean_val),end="\t")
    print("\nenum_amb\t%.1f" % (enum_amb["run_time"].mean()), end="\t")
    mean_val = enum_amb["run_time"].mean()
    for q in Q:
        print("%.1f" % (enum_amb["run_time"].quantile(q) / mean_val),end="\t")
    print("\narc_call\t%.1f" % (run_times_call["run_time"].mean()), end="\t")
    mean_val = run_times_call["run_time"].mean()
    for q in Q:
        print("%.1f" % (run_times_call["run_time"].quantile(q) / mean_val),end="\t")
    print("\narc_amb\t\t%.1f" % (run_times_amb["run_time"].mean()), end="\t")
    mean_val = run_times_amb["run_time"].mean()
    for q in Q:
        print("%.1f" % (run_times_amb["run_time"].quantile(q) / mean_val),end="\t")
    print("\n=====================================")

def read_run_times_new(n_ambs):
    file_name = "results/run_times_%d.txt" % (n_ambs)
    try:
        run_times = pd.read_csv(file_name, sep=" ", skiprows=1, index_col=False)
    except FileNotFoundError:
        return "%d" % (n_ambs) + "\t-\t-"*4, [math.inf, math.inf, math.inf, math.inf]
    
    run_times.columns = ["num_calls", "index_solver", "run_time"]
        
    enum_call = run_times[run_times["index_solver"] == 0]
    enum_amb = run_times[run_times["index_solver"] == 4]
    run_times_call = run_times[run_times["index_solver"] == 1]
    run_times_amb = run_times[run_times["index_solver"] == 2]

    sample_stds = [enum_call["run_time"].std() / math.sqrt(enum_call["run_time"].count()),
                enum_amb["run_time"].std() / math.sqrt(enum_amb["run_time"].count()),
                run_times_call["run_time"].std() / math.sqrt(run_times_call["run_time"].count()),
                run_times_amb["run_time"].std() / math.sqrt(run_times_amb["run_time"].count())]


    result = [str(n_ambs)]
    result.append("%.1f" % (enum_call["run_time"].mean()))
    result.append("%.1f" % (enum_call["run_time"].std()))
    result.append("%.1f" % (enum_amb["run_time"].mean()))
    result.append("%.1f" % (enum_amb["run_time"].std()))
    result.append("%.1f" % (run_times_call["run_time"].mean()))
    result.append("%.1f" % (run_times_call["run_time"].std()))
    result.append("%.1f" % (run_times_amb["run_time"].mean() - 700))
    result.append("%.1f" % (run_times_amb["run_time"].std() - 700))

    return "\t".join(result), sample_stds

def read_run_times_old(n_ambs):
    file_name = "results/results_runtimes/run_times_%d.txt" % (n_ambs)
    try:
        run_times = pd.read_csv(file_name, sep=" ", skiprows=1, index_col=False)
    except FileNotFoundError:
        return "%d" % (n_ambs) + "\t-\t-"*4, [math.inf,math.inf,math.inf]
    
    run_times.columns = ["num_calls", "index_solver", "run_time"]
        
    run_times_enum = run_times[run_times["index_solver"] == 0]
    run_times_call = run_times[run_times["index_solver"] == 1]
    run_times_amb = run_times[run_times["index_solver"] == 2]


    indexes = run_times[run_times["index_solver"].between(1,2,inclusive="both")]["index_solver"].copy(deep=True)
    # enum_call = pd.DataFrame(columns=["num_calls", "index_solver", "run_time"])
    # enum_amb = pd.DataFrame(columns=["num_calls", "index_solver", "run_time"])
    call_rows = []
    amb_rows = []
    for i,row in run_times_enum.iterrows():
        newrow = row
        if i % 2 == 0:
            # pd.concat([enum_call, row], axis=1)
            # enum_call.append(row, ignore_index=True)
            call_rows.append(newrow.values)
        else:
            # enum_amb.append(row, ignore_index=True)
            # pd.concat([enum_amb, row], axis=1)
            amb_rows.append(newrow.values)

    # enum_call = run_times_enum.append(pd.DataFrame(call_rows, columns=run_times_enum.columns)).reset_index()
    enum_call = pd.concat([run_times_enum, pd.DataFrame(call_rows, columns=run_times_enum.columns)]).reset_index(drop=True)
    # enum_amb = run_times_enum.append(pd.DataFrame(amb_rows, columns=run_times_enum.columns)).reset_index()
    enum_amb = pd.concat([run_times_enum, pd.DataFrame(amb_rows, columns=run_times_enum.columns)]).reset_index(drop=True)

    sample_stds = [enum_call["run_time"].std() / math.sqrt(enum_call["run_time"].count()),
                enum_amb["run_time"].std() / math.sqrt(enum_amb["run_time"].count()),
                run_times_call["run_time"].std() / math.sqrt(run_times_call["run_time"].count()),
                run_times_amb["run_time"].std() / math.sqrt(run_times_amb["run_time"].count())]


    result = [str(n_ambs)]
    result.append("%.1f" % (run_times_enum["run_time"].mean()))
    result.append("%.1f" % (run_times_enum["run_time"].std()))
    result.append("%.1f" % (enum_amb["run_time"].mean()))
    result.append("%.1f" % (enum_amb["run_time"].std()))
    result.append("%.1f" % (run_times_call["run_time"].mean()))
    result.append("%.1f" % (run_times_call["run_time"].std()))
    result.append("%.1f" % (run_times_amb["run_time"].mean() - 700))
    result.append("%.1f" % (run_times_amb["run_time"].std() - 700))

    return "\t".join(result), sample_stds
  

def main():
    if len(sys.argv) < 2:
        print("Usage: read_run_times [n_amb]")
        return
    
    all_n_ambs = [10]
    lines = []
    enum_call = []
    enum_amb = []
    model_call = []
    model_amb = []
    for n_ambs in all_n_ambs:
        result = read_run_times_old(n_ambs) if n_ambs < 30 and n_ambs >= 10 else read_run_times_new(n_ambs)
        tokens = result[0].split("\t")
        lines.append(result[0]+"\t"+"\t".join(["%.1f" % (x) for x in result[1]]))
        try:
            enum_call.append((float(tokens[1]), float(tokens[2])))
            enum_amb.append((float(tokens[3]), float(tokens[4])))
            model_call.append((float(tokens[5]), float(tokens[6])))
            model_amb.append((float(tokens[7]), float(tokens[8])))
        except:
            continue    
        

    for line in lines:
        print(line)

    for n_ambs in all_n_ambs:
        print_quantiles(n_ambs)
    
if __name__ == "__main__":
    main()