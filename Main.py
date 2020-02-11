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

'''Change these values as you wish'''

Number_of_Frames = 1000
Number_of_Balls = 20

'''Toggle these True/False to obtain graphs or switch the aniation off'''

Boltzmann_distribution = False
Pressure_graph = False
Energy_graph = False
animate = True




I = s(NoB = Number_of_Balls)
I.run(Number_of_Frames-1, animate)

with open('masses.data', 'rb') as filehandle:
   masses = pickle.load(filehandle)
with open('velocities.data', 'rb') as filehandle:
   velocities = pickle.load(filehandle)
with open('instantenous_pressure.data', 'rb') as filehandle:
   pressure = pickle.load(filehandle)
with open('time.data', 'rb') as filehandle:
   time = pickle.load(filehandle)

print('Simulation volume is', 2*np.pi*s.volume)
print('Number of particles is', Number_of_Balls)

eq_pressure = pressure[200:]  
eq_time = time[200:]    #only take values after reaching thermal equilibrium
avg_press = (sum(eq_pressure))/(sum(eq_time))
print('The calculated pressure is', avg_press)

energies = []
for i in range(Number_of_Frames-2):
    inst_energy = []
    for j in range(i*Number_of_Balls, (i+1)*Number_of_Balls):
        inst_energy.append(0.5*masses[j]*((velocities[j])**2))
    energies.append(sum(inst_energy))
    #total energy is same in every frame (conservation of energy), so
    #calculation slightly redundant, but just to perform a graphical check
tot_energy = sum(energies)/(Number_of_Frames-2)
kB = 1
Temperature = tot_energy/(Number_of_Balls*kB)
print('Temperature is', Temperature)

Energy_density = tot_energy/(np.pi*(s.volume)**2)
print('Energy density (actual pressure for an ideal gas) is', Energy_density)

cumulative_time = []
for i in range(len(time)):
    inst_time = []
    j=0
    while j<=i:
        inst_time.append(time[j])
        j+=1
    cumulative_time.append(sum(inst_time))

print('The total time of the simulation is', cumulative_time[-1])

if Energy_graph == True:
    f2 = plt.figure()
    plt.plot(energies)
    plt.title('Total energy of the system')
    plt.xlabel('Number of collisions since start')
    plt.ylabel('Energy')
    plt.show()
    
if Pressure_graph == True:
    f3 = plt.figure()
    plt.plot(cumulative_time,pressure)
    plt.title('Plot of instantenous pressure changes since start')
    plt.xlabel('Time since start')
    plt.ylabel('Instantenous pressure')
    plt.show()

if Boltzmann_distribution == True:
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

    f4 = plt.figure()
    bin_heights, bin_borders, _ = plt.hist(bolts_v, bins=20, label='histogram')
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    popt, _ = curve_fit(max_boltz, bin_centers, bin_heights)

    x_interval_for_fit = np.linspace(bin_borders[0], bin_borders[-1], 20)
    plt.plot(x_interval_for_fit, max_boltz(x_interval_for_fit, *popt),\
             label='fit')
    plt.legend()
    plt.title('Maxwell-Boltzmann distribution for our gas')
    plt.xlabel('Speed of a molecule')
    plt.ylabel('Number of instances')
    plt.show()
    print('Fit parameters are ([temperature, normalisation factor])', popt)