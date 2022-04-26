import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def plot_data(data, save_path=""):
    for l in data["length"].unique():
        print(l)
        plt.figure()
        curr_data = data[data["length"] == l]
        plt.plot(curr_data[curr_data["mode"]=="naive"]["k"], curr_data[curr_data["mode"]=="naive"]["size"], label="naive")
        plt.plot(curr_data[curr_data["mode"]=="simpaccel"]["k"], curr_data[curr_data["mode"]=="simpaccel"]["size"], label="simple accel")
        plt.xlabel("k")
        plt.ylabel("Index size (KB)")
        plt.title("Index size for reference length {}".format(l, "naive"))
        plt.legend()
        plt.savefig(save_path+"index_size_"+str(l)+".png")

size_file = "results/index_size.csv"

data = []

with open(size_file, "r") as f:
    for line in f.readlines():
        row = line.strip().split(",")
        if (row[0] == "rlen"): continue  # skip header row
        row = [int(row[0]), int(row[1]), row[2], int(row[3])]
        data.append(row)
     
df = pd.DataFrame(data, columns=["length", "k", "mode", "size"])

plot_data(df, "plots/")