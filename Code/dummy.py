import symrotor
import numpy as np
import asymrotor
import matplotlib.pyplot as plt
import os
import csv

ihod = [15.5, 0.0871, 0.0866]
uihod = [14.5, 0.0872, 0.0867]
lims = [-100,100]
freq, spec = asymrotor.spectra(ihod, uihod, [0.5,0,1],30,50,'ihoh', lims, 2, False, stats = [1,1,1,1])


filename = 'I-HOD_10K_lowres.csv'
spec2 = np.array([])
freq2 = np.array([])
if os.path.exists(filename):
    with open(filename, "r", newline="",encoding='utf-8-sig') as csv_file:
        reader = csv.reader(csv_file)
        for row in reader:
            freq2 = np.append(freq2, float(row[0]))
            spec2 = np.append(spec2, float(row[1]))

plt.figure()
plt.plot(freq + 3691, 10*spec, color='red', linewidth=0.5)
plt.plot(freq2, spec2, color='blue', linewidth=0.5)
plt.title('Rotational spectrum of ihod')
plt.xlabel('Energy (cm-1)')
plt.ylabel('Intensity')
plt.xlim([3600,3800])
plt.ylim([-5,30])
plt.show()

# ihoh = [16.785, 0.087166, 0.0842808]
# uihoh= [16.734, 0.08316, 0.08613]
# lims = [45,55]
# freq, spec = asymrotor.spectra(ihoh, uihoh, [0,1,0],20,10,'ihoh', lims, 0.01, False, stats = [1,1,1,1])
# plt.figure()
# plt.plot(freq, spec, color='red', linewidth=0.5)
# plt.title('Rotational spectrum of ihoh')
# plt.xlabel('Energy (cm-1)')
# plt.ylabel('Intensity')
# plt.show()