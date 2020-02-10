# -*- coding: utf-8 -*-
"""
Adam Brzosko
This is the main file used to run a simulation
of a gas in a container.
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np
from Classes import Simulation as s
from scipy.optimize import curve_fit

Number_of_Frames = 1000
Number_of_Balls = 150

Boltzmann_distribution = False

I = s(NoB = Number_of_Balls)
I.run(Number_of_Frames-1, True)

with open('masses.data', 'rb') as filehandle:
   masses = pickle.load(filehandle)
with open('velocities.data', 'rb') as filehandle:
   velocities = pickle.load(filehandle)
with open('instantenous_pressure.data', 'rb') as filehandle:
   pressure = pickle.load(filehandle)

#%%
avg_press = (sum(pressure))/(Number_of_Frames-1)
if Boltzmann_distribution == False:
    velocities_x = [x[0] for x in velocities]
    velocities_y = [y[1] for y in velocities]
    vx = [x**2 for x in velocities_x]
    vy = [y**2 for y in velocities_y]
    v = [a+b for a,b in zip(vx,vy)]
    bolts_v = []
    bolts_m = []
    for i in range(200*Number_of_Balls,len(v)):
        bolts_v.append(np.sqrt(v[i]))
        bolts_m.append(masses[i])
def max_boltz(v, T, k):
    return k*v*(np.exp(-(0.5*(v*v))/T))
'''   
boltzmann = plt.figure()
ax1 = boltzmann.add_subplot(1, 1, 1)
plt.hist(bolts_v, bins = 30)       
plt.title('Maxwell-Boltzmann distribution for our gas')
ax1.set_xlabel('Energy of a molecule')
ax1.set_ylabel('Number of instances')
plt.show()'''

bin_heights, bin_borders, _ = plt.hist(bolts_v, bins=30, label='histogram')
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
popt, _ = curve_fit(max_boltz, bin_centers, bin_heights)

x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 30)
plt.plot(x_interval_for_fit, max_boltz(x_interval_for_fit, *popt), label='fit')
plt.legend()

avg_press = (sum(pressure))/1000
avg_press