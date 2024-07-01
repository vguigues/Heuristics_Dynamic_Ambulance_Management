#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as df
import math
from itertools import cycle


def get_data(prefix, nb_ambs, heuristic, setup, num_scenarios, best_base):
    scen, which_amb, w_scene, w_pen, end = [], [], [], [], []
    valid = []
    bb = "_bb" if best_base else ""
    filename = f"{prefix}/{setup}_{heuristic}_{nb_ambs}_{num_scenarios}{bb}.txt"
    try:
        arq = open(filename, "r")
    except FileNotFoundError:
        print(f"Missing experiment {filename}")
        return filename
    # print(f"Reading {filename}")
    for s in range(num_scenarios):
        num_calls = int(arq.readline().strip())
        w_amb, w_s, w_p, this_end = [], [], [], []
        this_valid = True
        for i in range(num_calls):
            tokens = arq.readline().split()
            wa = int(tokens[0])
            ws = float(tokens[1])
            wp = float(tokens[2])
            e = float(tokens[3])

            w_amb.append(wa)
            w_s.append(ws)
            w_p.append(wp)
            this_end.append(e)
            if wa == -1:
                this_valid = False
        this_scen = [s] * num_calls
        scen += this_scen
        which_amb += w_amb
        w_scene += w_s
        w_pen += w_p
        end += this_end
        valid += [this_valid] * num_calls

    arq.close()

    return df.DataFrame(
        {
            "scen": scen,
            "w_scene": w_scene,
            "w_pen": w_pen,
            "which_amb": which_amb,
            "end": end,
            "s_complete": valid,
        }
    )


# Gráfico para número de ambulâncias
# Média, máximo, q(0.9) para penalizado e não penalizado


def write_table(algo, param, nb_ambs, heuristics, setup, num_scenarios):
    best_base = True
    has_markov = "markov_prep" in heuristics
    bb = "bb" if best_base else ""
    arq = open(f"tables/{algo}/{setup}_{param}_{num_scenarios}_{bb}.txt", "w")
    arq.write("\t".join(heuristics) + "\n")
    arq.write("A\t" + ("mean\tQ(0.9)\tmax\t" * len(heuristics)) + "\n")
    # arq.write("Heuristic\t")
    # for nb_amb in nb_ambs:
    #     arq.write(f"A = {nb_amb}\t" * 3)
    # arq.write("\n")
    for nb_amb in nb_ambs:
        line = []
        line.append("%d" % nb_amb)
        for heuristic in heuristics:
            data = get_data(algo, nb_amb, heuristic, setup, num_scenarios, best_base)
            line.append("%.0f" % data[data["which_amb"] != -1][param].mean())
            line.append("%.0f" % data[data["which_amb"] != -1][param].quantile(0.9))
            line.append("%.0f" % data[data["which_amb"] != -1][param].max())
        full_line = "\t".join(line)
        arq.write(full_line + "\n")
    arq.close()


heuristics = [
    "queue",
    "forward",
    "priorities",
    "minmax",
    "non_miopyc",
    "preparedness",
    "prep2",
    "district",
    "ordered",
    "markov_prep",
]
LEGENDS = {
    "queue": "CA",
    "forward": "BM",
    "priorities": "GHP1",
    "minmax": "GHP2",
    "non_miopyc": "NM",
    "preparedness": "Lee",
    "prep2": "Andersson",
    "district": "Mayorga",
    "ordered": "Bandara",
    "markov_prep": "Markov",
}


def plot_results(algo, param, nb_ambs, heuristics, setup, num_scenarios, best_base):
    ambs = []
    means = {}
    quantiles = {}
    maximums = {}

    best_base_compatible = ["forward", "priorities", "minmax", "non_myopic"]

    for heuristic in heuristics:
        means[heuristic] = []
        quantiles[heuristic] = []
        maximums[heuristic] = []

    for nb_amb in nb_ambs:
        line = []
        ambs.append(nb_amb)
        for heuristic in heuristics:
            data = get_data(
                algo,
                nb_amb,
                heuristic,
                setup,
                num_scenarios,
                best_base and heuristic in best_base_compatible,
            )
            if not isinstance(data, str):
                means[heuristic].append(data[data["which_amb"] != -1][param].mean())
                quantiles[heuristic].append(
                    data[data["which_amb"] != -1][param].quantile(0.9)
                )
                maximums[heuristic].append(data[data["which_amb"] != -1][param].max())
            else:
                print(f"Warning: Experiment {data} is missing!")
                means[heuristic].append("-")
                quantiles[heuristic].append("-")
                maximums[heuristic].append("-")
    for heuristic in heuristics:
        print(f"algo {algo}, heuristic {heuristic}, means = {means[heuristic]}")
    fig, ax = plt.subplots(figsize=(12, 10))
    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)
    for i, heuristic in enumerate(heuristics):
        plt.plot(
            ambs, means[heuristic], linestyle=next(linecycler), label=LEGENDS[heuristic]
        )
    m_font_size = 28
    tick_size = 20
    ax.set_xlabel("# Ambulances", fontsize=m_font_size)
    ylabel = "PRT" if param == "w_pen" else "RT"
    ax.set_ylabel(f"Mean {ylabel}", fontsize=m_font_size)
    ax.legend(fontsize=m_font_size)
    new_list = range(math.floor(min(ambs)), math.ceil(max(ambs)) + 1)
    plt.xticks(new_list, fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.savefig(
        f"tables/{algo}/{setup}_{param}_mean_{num_scenarios}{'_bb' if best_base else ''}_{algo}.pdf",
        bbox_inches="tight",
    )
    fig, ax = plt.subplots(figsize=(12, 10))
    for i, heuristic in enumerate(heuristics):
        plt.plot(
            ambs,
            maximums[heuristic],
            linestyle=next(linecycler),
            label=LEGENDS[heuristic],
        )
    ax.set_xlabel("# Ambulances", fontsize=m_font_size)
    ylabel = "PRT" if param == "w_pen" else "RT"
    ax.set_ylabel(f"Max {ylabel}", fontsize=m_font_size)
    ax.legend(fontsize=m_font_size)
    new_list = range(math.floor(min(ambs)), math.ceil(max(ambs)) + 1)
    plt.xticks(new_list, fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.savefig(
        f"tables/{algo}/{setup}_{param}_max_{num_scenarios}{'_bb' if best_base else ''}_{algo}.pdf",
        bbox_inches="tight",
    )
    fig, ax = plt.subplots(figsize=(12, 10))
    for i, heuristic in enumerate(heuristics):
        plt.plot(
            ambs,
            quantiles[heuristic],
            linestyle=next(linecycler),
            label=LEGENDS[heuristic],
        )
    ax.set_xlabel("# Ambulances", fontsize=m_font_size)
    ylabel = "PRT" if param == "w_pen" else "RT"
    ax.set_ylabel(f"Q0.9 {ylabel}", fontsize=m_font_size)
    ax.legend(fontsize=m_font_size)
    new_list = range(math.floor(min(ambs)), math.ceil(max(ambs)) + 1)
    plt.xticks(new_list, fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.savefig(
        f"tables/{algo}/{setup}_{param}_q09_{num_scenarios}{'_bb' if best_base else ''}_{algo}.pdf",
        bbox_inches="tight",
    )


def display_results(args):
    setup = args[-2]
    num_scenarios = int(args[-1])
    nb_ambs = []

    heuristics = []
    if args[0] == "markov":
        nb_ambs = [6, 8, 10]
        heuristics = [
            "queue",
            "markov_prep",
            "preparedness",
            "district",
            "ordered",
        ]
        heuristics.append("markov_prep")
    else:
        heuristics = [
            "forward",
            "priorities",
            "minmax",
            "non_miopyc",
            # "preparedness",
            # "district",
            # "ordered",
        ]
        nb_ambs = [10, 12, 14, 16, 18, 20]
        # nb_ambs = [10]

    write_table("one_stage", "w_scene", nb_ambs, heuristics, setup, num_scenarios)
    write_table("one_stage", "w_pen", nb_ambs, heuristics, setup, num_scenarios)
    write_table("two_stage", "w_scene", nb_ambs, heuristics, setup, num_scenarios)
    write_table("two_stage", "w_pen", nb_ambs, heuristics, setup, num_scenarios)

    # best_base = False
    # plot_results(
    #     "one_stage", "w_pen", nb_ambs, heuristics, setup, num_scenarios, best_base
    # )
    # plot_results(
    #     "two_stage", "w_pen", nb_ambs, heuristics, setup, num_scenarios, best_base
    # )
    # plot_results(
    #     "one_stage", "w_scene", nb_ambs, heuristics, setup, num_scenarios, best_base
    # )
    # plot_results(
    #     "two_stage", "w_scene", nb_ambs, heuristics, setup, num_scenarios, best_base
    # )


def validate_args(args):
    if len(args) > 4:
        print(
            "ERROR: at most 3 arguments expected: python plot_vals.py [markov] amb_setup num_scenarios"
        )
        exit(1)
    elif len(args) == 4:
        if args[1] != "markov":
            print(
                "ERROR: if 3 arguments provided, first argument must be markov: python plot_vals.py markov amb_setup num_scenarios"
            )
            exit(2)
    elif len(args) < 4:
        if not args[1] in ["us", "rj"]:
            print(
                "ERROR: supported setups are us and rj: python plot_vals.py [markov] (us|rj) num_scenarios"
            )
            exit(3)


def main():
    validate_args(sys.argv)
    display_results(sys.argv[1:])


if __name__ == "__main__":
    main()
