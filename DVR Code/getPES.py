import os
import csv
import numpy as np

def getPES(filename):
    grid = np.zeros(0)
    PE = np.zeros(0)
    if os.path.exists(filename):
        with open(filename, "r", newline="") as csv_file:
            reader = csv.reader(csv_file, delimiter=",")
            for row in reader:
                grid = np.append(grid, float(row[0]))
                PE = np.append(PE, float(row[1]))
    return grid, PE
