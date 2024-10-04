import symrotor
import numpy as np
import asymrotor
import matplotlib.pyplot as plt
import os
import csv

filename = 'IHODBinnedOutput.csv'
spec2 = np.array([])
# filename = 'I-HOD_10K_lowres.csv'
freq2 = np.array([])
if os.path.exists(filename):
    with open(filename, "r", newline="",encoding='utf-8-sig') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            freq2 = np.append(freq2, float(row[0]))
            spec2 = np.append(spec2, float(row[1]))
spec2 = spec2/np.max(spec2)

shift = 3687.037+0.0483
ihod = [17.177, 0.085704, 0.0827841]
uihod = [16.517, 0.0881541, 0.0841155]
print(ihod)
print(uihod)
lims = [np.min(freq2)-shift,np.max(freq2)-shift]
# freq, spec = asymrotor.spectra(ihod, uihod, [0.3,1,0],30,10,'ihoh', lims, 0.01, False, stats = [1,1,1,1])
# spec = spec/np.max(spec)


fig, axs = plt.subplots(nrows=2, ncols=1, figsize = (8,5))

# plt.plot( freq + shift, -spec, color='red', linewidth=0.5)
ax = axs[1]
ax.plot(freq2, spec2, color='blue', linewidth=0.5)
ax.set_title('Rotational spectrum of I-HOD')
ax.set_xlabel('Energy (cm-1)')
ax.set_ylabel('Intensity')

spec2 = np.array([])
filename = 'I-HOD_10K_lowres.csv'
freq2 = np.array([])
if os.path.exists(filename):
    with open(filename, "r", newline="",encoding='utf-8-sig') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            freq2 = np.append(freq2, float(row[0]))
            spec2 = np.append(spec2, float(row[1]))
spec2 = spec2/np.max(spec2)

ax = axs[0]
ax.plot(freq2, spec2, color='blue', linewidth=2)
ax.set_title('Vibrational spectrum of I-HOD')
ax.set_xlabel('Energy (cm-1)')
ax.set_ylabel('Intensity')
# Highlight a vertical region between x=2 and x=4
ax.axvspan(3700, 3708, color='green', alpha=0.3)

plt.subplots_adjust(top=0.88, hspace=0.45)
plt.show()

