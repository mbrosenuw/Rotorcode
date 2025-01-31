import symrotor
import numpy as np
import asymrotor
import matplotlib.pyplot as plt
import os
import csv

# filename = 'IHODBinnedOutput.csv'
# spec2 = np.array([])
# # filename = 'I-HOD_10K_lowres.csv'
# freq2 = np.array([])
# if os.path.exists(filename):
#     with open(filename, "r", newline="",encoding='utf-8-sig') as csv_file:
#         reader = csv.reader(csv_file)
#         for row in reader:
#             freq2 = np.append(freq2, float(row[0]))
#             spec2 = np.append(spec2, float(row[1]))
#
# c1 = (3703.497 < freq2) & (freq2 < 3703.883)
# freq2 = freq2[c1]
# spec2 = spec2[c1]
# spec2 = spec2/np.max(spec2)
#
# shift = 3687.037+0.0483
# ihod = [17.177, 0.085704, 0.0827841]
# uihod = [16.517, 0.0881541, 0.0841155]
#
#
# lims = [np.min(freq2)-shift,np.max(freq2)-shift]
# freq, spec = asymrotor.spectra(ihod, uihod, [0.3,1,0],20,10,'ihoh', lims, 0.003, False, stats = [1,1,1,1])
# # freq, spec = asymrotor.spectranoLS(ihod, uihod, [0.3,1,0],30,10,'ihoh', lims, 0.005, False, stats = [1,1,1,1])
# spec = spec/np.max(spec)
# fig = plt.figure(figsize=(9,5))
# plt.plot(freq+shift, -spec, color = 'red', label = 'Theory', linewidth = 0.5)
# plt.plot(freq2, spec2, color = 'blue', label = 'Experiment', linewidth = 0.5)
# plt.title('Q Branch of I$^-$HOD', fontsize=18)
# plt.ylabel('Intensity', fontsize = 14)
# plt.xlabel('Energy $[cm^{-1}]$', fontsize = 14)
# plt.legend(loc = 'best')
# plt.xlim([3703.497, 3703.883])
# plt.show()

#
# ihod = [17.177, 0.085704, 0.0827841]
# uihod = [16.517, 0.0881541, 0.0841155]
# freq2, spec2 = asymrotor.spectra([0.08], s1, [0,0,1], jmax=20, T = 100, name="prolate", lims = [-100,100], width = 0.1, showE=False)
# fig, axs = plt.subplots(nrows=2, ncols=1, figsize = (8,5))

# plt.plot( freq + shift, -spec, color='red', linewidth=0.5)
# ax = axs[1]
# ax.plot(freq2, spec2, color='blue', linewidth=0.5)
# ax.set_title('Rotational spectrum of I$^-$HOD', fontsize=18)
# plt.ylabel('Intensity', fontsize = 14)
# plt.xlabel('Energy $[cm^{-1}]$', fontsize = 14)
# plt.legend(loc = 'best')
# plt.show()
#
# spec2 = np.array([])
# filename = 'I-HOD_10K_lowres.csv'
# freq2 = np.array([])
# if os.path.exists(filename):
#     with open(filename, "r", newline="",encoding='utf-8-sig') as csv_file:
#         reader = csv.reader(csv_file)
#         for row in reader:
#             freq2 = np.append(freq2, float(row[0]))
#             spec2 = np.append(spec2, float(row[1]))
# spec2 = spec2/np.max(spec2)
#
# # s2 = [20,20,1]
# # freq, spec = asymrotor.spectra(s2, s2, [0,0,1], jmax=20, T = 100, name="oblate", lims = [-300, 300], width = 0.1, showE=False)
# ax = axs[0]
# ax.plot(freq2, spec2, color='blue', linewidth=2)
# ax.set_title('Vibrational spectrum of I-HOD', fontsize= 18)
# ax.set_xlabel('Energy (cm-1)', fontsize= 14)
# ax.set_ylabel('Intensity', fontsize= 14)
# # Highlight a vertical region between x=2 and x=4
# ax.axvspan(3700, 3708, color='green', alpha=0.3)
#
# plt.subplots_adjust(top=0.88, hspace=0.75)
# plt.show()

dms = [0.59406, 0.25421,0.19073]
mu = [0,1,0]
jmax = 20
T = 20
lims = [-12,12]
# width = 0.0067*1.5
width = 0.03
freq, spec = asymrotor.spectra(dms, dms, mu,jmax,T,'Dimethyl Sulfide', lims, width, False, stats = [1,1,1,1])
spec = spec/np.max(spec)
fig = plt.figure(figsize=(7,5))
shift = 1033
plt.plot(shift+freq, spec, color = 'red', label = 'Theory', linewidth = 0.5)
plt.title('Rotational Spectrum of Dimethyl Sulfide', fontsize=18)
plt.ylabel('Intensity', fontsize = 14)
plt.xlabel('Energy $[cm^{-1}]$', fontsize = 14)
plt.legend(loc = 'best')
plt.xlim(shift+np.array(lims))
plt.show()
np.savez('dms_purerotor2.npz', freq = freq, spec = spec, shift = shift)