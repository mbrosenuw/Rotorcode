import symrotor
import numpy as np
import asymrotor
from rotor import Rotor
chl = np.array([0.322, 0.322, 5.246])
# asymrotor.spectra(chl, chl, [0, 0, 1], 30, 10, 'chloromethane', (-40, 40), 0.05,False)

water = [27.28,14.5,9.28]
# water = [14.5,27.28,9.95]
sys = Rotor(4, water)
print(np.round(sys.oldH, decimals=2))
# uwater = water - np.multiply(0.05,water)
# print(water, uwater)
# asymrotor.spectra(water,uwater, [1,1,1],20,100,'water', [-200,200], 0.5, False)


# ihoh = [16.785, 0.08717, 0.08428]
# uihoh = [16.734, 0.08613, 0.08316]
# asymrotor.spectra(ihoh, uihoh, [0,0,10000],20,100,'ihoh', [-4,4], 0.005, False)
ihoh = [0.08717, 0.08428, 16.785]
uihoh = [0.08613, 0.08316, 16.734]
sys = Rotor(2, ihoh)
for wfn in sys.wfns:
    print(np.round(wfn,5))
# symrotor.spectra(ihoh, uihoh, [0,0,1],20,10,'water', [-4,4], 0.005, True)