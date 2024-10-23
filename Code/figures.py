import numpy as np
import asymrotor
import matplotlib.pyplot as plt
from rotor import Rotor
from matplotlib.pyplot import cm
import os
import csv
import itertools

def getRAP():
    jmax=3
    jlist = list(range(jmax+1))
    a = 20
    c = 0.5
    brange = np.linspace(c, a, 100)
    energylist = []
    k = []
    for b in brange:
        r2 = Rotor(jmax, [a,b,c])
        r2max = np.max(r2.fullenergies)
        energies = []
        k.append((2*b-a-c)/(a-c))
        for j in jlist:
            r = Rotor(j, [a, b, c])
            energies.append(np.sort(r.fullenergies)/r2max)
        if not energylist:
            energylist = energies
        else:
            energylist = [np.vstack((matrix, new_row)) for matrix, new_row in zip(energylist, energies)]


    # colors = cm.get_cmap('tab10', jmax+1)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#b03060']
    color = itertools.cycle(colors)
    # color = iter([colors(i) for i in range(colors.N)])
    plt.figure(figsize=(10,6.5))
    for (matrix, j) in zip(energylist,jlist):
        c = next(color)
        for row in matrix.T:
            plt.plot(k,row, color=c, label='J = ' +str(j))

    plt.title('Energy levels vs Asymmetry Parameter', fontsize=18)
    plt.xlabel('Asymmetry Parameter [unitless]', fontsize=16)
    plt.ylabel('Relative Energy', fontsize=16)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=14)
    plt.xlim([-1,0] )
    plt.show()


def getspec():
    filename = 'IHODBinnedOutput.csv'
    spec2 = np.array([])
    # filename = 'I-HOD_10K_lowres.csv'
    freq2 = np.array([])
    if os.path.exists(filename):
        with open(filename, "r", newline="", encoding='utf-8-sig') as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                freq2 = np.append(freq2, float(row[0]))
                spec2 = np.append(spec2, float(row[1]))
    spec2 = spec2 / np.max(spec2)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 5))

    # plt.plot( freq + shift, -spec, color='red', linewidth=0.5)
    ax = axs[1]
    ax.plot(freq2, spec2, color='blue', linewidth=0.5)
    ax.set_title('Rotational spectrum of I-HOD', fontsize=18)
    ax.set_xlabel('Energy (cm-1)', fontsize=14)
    ax.set_ylabel('Intensity', fontsize=14)

    spec2 = np.array([])
    filename = 'I-HOD_10K_lowres.csv'
    freq2 = np.array([])
    if os.path.exists(filename):
        with open(filename, "r", newline="", encoding='utf-8-sig') as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                freq2 = np.append(freq2, float(row[0]))
                spec2 = np.append(spec2, float(row[1]))
    spec2 = spec2 / np.max(spec2)

    ax = axs[0]
    ax.plot(freq2, spec2, color='blue', linewidth=2)
    ax.set_title('Vibrational spectrum of I-HOD', fontsize=18)
    ax.set_xlabel('Energy (cm-1)', fontsize=14)
    ax.set_ylabel('Intensity', fontsize=14)
    # Highlight a vertical region between x=2 and x=4
    ax.axvspan(3700, 3708, color='green', alpha=0.3)

    plt.subplots_adjust(top=0.88, hspace=0.75)
    plt.show()

getRAP()