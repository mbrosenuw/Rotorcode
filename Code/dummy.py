import symrotor
import numpy as np
import asymrotor
import matplotlib.pyplot as plt

ihod = [15.5, 0.0871, 0.0866]
uihod = [14.5, 0.0872, 0.0867]
freq, spec = asymrotor.spectra(ihod, uihod, [1,1,1],20,10,'ihoh', [-50,50], 2, False, stats = [1,1,1,1])

plt.figure()
plt.plot(freq + 3691, spec, color='red', linewidth=0.5)
plt.title('Rotational spectrum of ihod')
plt.xlabel('Energy (cm-1)')
plt.ylabel('Intensity')
plt.show()
