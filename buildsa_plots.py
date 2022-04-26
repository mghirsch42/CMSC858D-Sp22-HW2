import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_data(data, save_path=""):
    for l in data["length"].unique():
        print(l)
        plt.figure()
        curr_data = data[data["length"] == l]
        plt.plot(curr_data[curr_data["mode"]=="naive"]["k"], curr_data[curr_data["mode"]=="naive"]["time"], label="naive")
        plt.plot(curr_data[curr_data["mode"]=="simpaccel"]["k"], curr_data[curr_data["mode"]=="simpaccel"]["time"], label="simple accel")
        plt.xlabel("k")
        plt.ylabel("Build time (ms)")
        plt.title("Build time for reference length {}".format(l, "naive"))
        plt.legend()
        plt.savefig(save_path+"build_time_"+str(l)+".png")

time_file = "results/build_times.csv"

data = []

with open(time_file, "r") as f:
    for line in f.readlines():
        row = line.strip().split(",")
        row = [int(row[0]), int(row[1]), row[2], float(row[3])]
        data.append(row)
     
df = pd.DataFrame(data, columns=["length", "k", "mode", "time"])

plot_data(df, "plots/")