import symrotor
import numpy as np
import asymrotor
chl = np.array([0.322, 0.322, 5.246])
asymrotor.spectra(chl, chl, [0, 0, 1], 30, 10, 'chloromethane', (-40, 40), 0.05,False)

# water = [27.28,14.5,9.95]
# # water = [14.5,27.28,9.95]
# uwater = water - np.multiply(0.05,water)
# print(water, uwater)
# asymrotor.spectra(water,water, [1,0,0],30,100,'water', [-200,200], 0.5, False)
