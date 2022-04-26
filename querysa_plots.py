import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_k_data(data, qtype, save_path=""):
    for l in data["length"].unique():
        print(l)
        plt.figure()
        curr_data = data[data["length"] == l]
        # print(curr_data.head())
        ks = [0, 5, 10]
        averaged_naive_data = [
            curr_data[(curr_data["k"] == 0) & (curr_data["mode"] == "naive")]["time"].mean(),
            curr_data[(curr_data["k"] == 5) & (curr_data["mode"] == "naive")]["time"].mean(),
            curr_data[(curr_data["k"] == 10) & (curr_data["mode"] == "naive")]["time"].mean(),
        ]
        averaged_accel_data = [
            curr_data[(curr_data["k"] == 0) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
            curr_data[(curr_data["k"] == 5) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
            curr_data[(curr_data["k"] == 10) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
        ]
        plt.plot(ks, averaged_naive_data, label="naive")
        plt.plot(ks, averaged_accel_data, label="simple accel")
        plt.xlabel("k")
        plt.ylabel("Query time (ms)")
        if (qtype == "rand"):
            title = "Query time for reference length {}, random queries".format(l)
        else:
            title = "Query time for reference length {}, reference queries".format(l)
        plt.title(title)
        plt.legend()
        plt.savefig(save_path+"query_time_rlen_"+str(l)+"_k_"+qtype+".png")


def plot_qlen_data(data, qtype, save_path=""):
    for l in data["length"].unique():
        print(l)
        plt.figure()
        curr_data = data[data["length"] == l]
        # print(curr_data.head())
        qlens = [10, 100]
        averaged_naive_data = [
            curr_data[(curr_data["q_len"] == 10) & (curr_data["mode"] == "naive")]["time"].mean(),
            curr_data[(curr_data["q_len"] == 100) & (curr_data["mode"] == "naive")]["time"].mean(),
        ]
        averaged_accel_data = [
            curr_data[(curr_data["q_len"] == 10) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
            curr_data[(curr_data["q_len"] == 100) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
        ]
        plt.plot(qlens, averaged_naive_data, label="naive")
        plt.plot(qlens, averaged_accel_data, label="simple accel")
        plt.xlabel("Query Length")
        plt.ylabel("Query time (ms)")
        if (qtype == "rand"):
            title = "Query time for reference length {}, random queries".format(l)
        else:
            title = "Query time for reference length {}, reference queries".format(l)
        plt.title(title)        
        plt.legend()
        plt.savefig(save_path+"query_time_rlen_"+str(l)+"_qlen_"+qtype+".png")

def plot_qn_data(data, qtype, save_path=""):
    for l in data["length"].unique():
        print(l)
        plt.figure()
        curr_data = data[data["length"] == l]
        # print(curr_data.head())
        qn = [100, 1000, 10000]
        averaged_naive_data = [
            curr_data[(curr_data["q_n"] == 100) & (curr_data["mode"] == "naive")]["time"].mean(),
            curr_data[(curr_data["q_n"] == 1000) & (curr_data["mode"] == "naive")]["time"].mean(),
            curr_data[(curr_data["q_n"] == 10000) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
        ]
        averaged_accel_data = [
            curr_data[(curr_data["q_n"] == 100) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
            curr_data[(curr_data["q_n"] == 1000) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
            curr_data[(curr_data["q_n"] == 10000) & (curr_data["mode"] == "simpaccel")]["time"].mean(),
        ]
        plt.plot(qn, averaged_naive_data, label="naive")
        plt.plot(qn, averaged_accel_data, label="simple accel")
        plt.xlabel("Query Length")
        plt.ylabel("Query time (ms)")
        plt.xscale("log")
        if (qtype == "rand"):
            title = "Query time for reference length {}, random queries".format(l)
        else:
            title = "Query time for reference length {}, reference queries".format(l)
        plt.title(title)        
        plt.legend()
        plt.savefig(save_path+"query_time_rlen_"+str(l)+"_qn_"+qtype+".png")

for qtype in ["ref", "rand"]:
    time_file = "results/query_times_{}.csv".format(qtype)

    data = []

    with open(time_file, "r") as f:
        for line in f.readlines():
            row = line.strip().split(",")
            row = [int(row[0]), int(row[1]), row[2], int(row[3]), int(row[4]), float(row[5])]
            data.append(row)
        
    df = pd.DataFrame(data, columns=["length", "k", "mode", "q_len", "q_n", "time"])

    plot_k_data(df, qtype, "plots/")
    plot_qlen_data(df, qtype, "plots/")
    plot_qn_data(df, qtype, "plots/")