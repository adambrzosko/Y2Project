#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 20:57:14 2020
@author: Adam Brzosko
Y2 project

This is a test file, containing snippets of code that were
used to test the code used in the main part 
"""

from classes import Simulation as s

I = s()

I.run(10, True)

#%%
import random as rand
NoB = 3
x = []
y = []
for i in range(-10,11,2):    #sqrt of half of the square of radius, step two to avoid overlap
    x.append(i)
    y.append(i)
for i in range(NoB):
    pos_x = rand.choice(x)
    pos_y = rand.choice(y)
    print(pos_x)
    print(pos_y)
    x.remove(pos_x)
    y.remove(pos_y)
    print(x)
    print(y)
#%%
from Classes import Ball
import numpy as np
def collide(self, other):
    '''change vel components parallel to the axis
    containing centres of the balls in both balls'''
    pos_self_vector = other._pos - self._pos  #vectors along the axis
    pos_other_vector = - pos_self_vector       
    
    paral_vel_self = ((np.dot(self._vel,pos_self_vector))*(pos_self_vector))/(np.dot(pos_self_vector,pos_self_vector))   #component splitting self
    perp_vel_self = self._vel - paral_vel_self      
        
    if other._mass == float('inf'):
        new_paral_vel_other = 0
        perp_vel_other = 0
        new_paral_vel_self = - paral_vel_self
    else:
        paral_vel_other = ((np.dot(other._vel,pos_other_vector))*(pos_other_vector))/(np.dot(pos_other_vector,pos_other_vector))    #component splitting other
        perp_vel_other = other._vel - paral_vel_other
        new_paral_vel_self = (((self._mass - other._mass)/(self._mass + other._mass))*paral_vel_self) + (((2*other._mass)/(self._mass + other._mass))*paral_vel_other)  #calculate new component self
        new_paral_vel_other = (((2*self._mass)/(self._mass + other._mass))*paral_vel_self) + (((other._mass-self._mass)/(self._mass + other._mass))*paral_vel_other)   #calculate new component other
    self._vel = perp_vel_self + new_paral_vel_self
    other._vel = perp_vel_other + new_paral_vel_other    #change the componenets


b1 = Ball(r=np.array([4,4]), v=np.array([-1,-1]))
b2 = Ball(v=np.array([1,1]))
t = b1.time_to_collision(b2)
print(t)
print(b1)
print(b2)
b1.move(t)
b2.move(t)
b1.collide(b2)
print(b1)
print(b2)
#%%
import numpy as np
import pylab as pl
import random as rand
from Classes import Ball

class Simulation:
    ball_list = []
    volume = 20
    
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 11
        self._container = Ball(container=True)
        x = []
        y = []
        for i in range(-10,11,2):    #sqrt of half of the square of radius, step three (two is touching) to avoid overlap spawn
            x.append(i)
            y.append(i)
        for i in range(self._NoB):
            pos_x = rand.choice(x)
            pos_y = rand.choice(y)
            Simulation.ball_list.append(Ball(r= np.array([pos_x,pos_y]), v=np.array([rand.choice(x),rand.choice(y)])))
            x.remove(pos_x)
            y.remove(pos_y)

    def next_collision(self):
        tb = []
        t = Simulation.ball_list[self._NoB-1].time_to_collision(self._container)
        print('t is', t)
        ball1_id = self._NoB-1    #set last ball as default to collide with a container
        Container_collision = True
        for k in range(self._NoB-1):     #calculate all other balls times to collisions
            tc = Simulation.ball_list[k].time_to_collision(self._container)  #first collision time with container
            print('tc is:', tc)
            for i in range(k+1, self._NoB):  #then with all balls after it on the list
                tb.append(Simulation.ball_list[k].time_to_collision(Simulation.ball_list[i]))
                if tb != [] and tb[i-k-1] < t and tb[i-k-1] < tc and tb[i-k-1] > 1e-10:
                    t = tb[i-k-1]
                    ball1_id = k
                    ball2_id = i
                    Container_collision = False
                    print(k,i)
                elif tc < t and tc > 1e-10:
                    t = tc
                    ball1_id = k
                    Container_collision = True
            print('tb is:', tb)
            tb = []
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t)
        if Container_collision == False:
            Simulation.ball_list[ball1_id].collide(Simulation.ball_list[ball2_id])
            print('ball', ball1_id, 'collides with', ball2_id)
        else:
           Simulation.ball_list[ball1_id].collide(self._container)
           print('ball', ball1_id, 'collides with container')
        
    def run(self, num_frames, animate=False):
        if animate:
            ax = pl.axes(xlim=(-Simulation.volume, Simulation.volume), ylim=(-Simulation.volume, Simulation.volume))
            ax.add_artist(pl.Circle([0., 0.], Simulation.volume, ec='b', fill=False, ls='solid'))
            p = []
            for i in range(self._NoB):
                p.append(ax.add_patch(Simulation.ball_list[i].get_patch()))
        for frame in range(num_frames):
            self.next_collision()
            '''Energy = (self._ball._mass*np.dot(self._ball._vel, self._ball._vel))/2
            print(Energy)
            Pressure = (Energy)/(np.pi*(self._container._rad)**2)
            print(Pressure)'''
            if animate:
                for i in range(self._NoB):
                    print(p[i])
                    p[i].center = (Simulation.ball_list[i]._pos)
                    print(p[i])
                    pl.pause(0.001)

I = Simulation()
I.run(20, True)
#%%
import matplotlib.pyplot as plt

class Simulation:
    ball_list = []
    volume = 20
    pressure_hits = []
    T = 0
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 11
        self._container = Ball(container=True)
        x = []
        y = []
        for i in range(-10,11,2):    #sqrt of half of the square of radius, step three (two is touching) to avoid overlap spawn
            x.append(i)
            y.append(i)
        for i in range(self._NoB):
            pos_x = rand.choice(x)
            pos_y = rand.choice(y)
            Simulation.ball_list.append(Ball(r= np.array([pos_x,pos_y]), v=np.array([rand.choice(x),rand.choice(y)])))
            x.remove(pos_x)
            y.remove(pos_y)

    def next_collision(self):
        tb = []
        t = Simulation.ball_list[self._NoB-1].time_to_collision(self._container)
        print('t is', t)
        ball1_id = self._NoB-1    #set last ball as default to collide with a container
        Container_collision = True
        for k in range(self._NoB-1):     #calculate all other balls times to collisions
            tc = Simulation.ball_list[k].time_to_collision(self._container)  #first collision time with container
            print('tc is:', tc)
            for i in range(k+1, self._NoB):  #then with all balls after it on the list
                tb.append(Simulation.ball_list[k].time_to_collision(Simulation.ball_list[i]))
                if tb != [] and tb[i-k-1] < t and tb[i-k-1] < tc and tb[i-k-1] > 1e-15:
                    t = tb[i-k-1]
                    ball1_id = k
                    ball2_id = i
                    Container_collision = False
                    print(k,i)
                elif tc < t and tc > 1e-15:
                    t = tc
                    ball1_id = k
                    Container_collision = True
            print('tb is:', tb)
            tb = []
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t)
        if Container_collision == False:
            Simulation.ball_list[ball1_id].collide(Simulation.ball_list[ball2_id])
            print('ball', ball1_id, 'collides with', ball2_id)
            Simulation.pressure_hits.append(0)
        else:
           Simulation.ball_list[ball1_id].collide(self._container)
           Simulation.pressure_hits.append((2*Simulation.ball_list[ball1_id]._mass*np.sqrt(abs(np.dot(Simulation.ball_list[ball1_id]._vel,Simulation.ball_list[ball1_id]._vel))))/(2*np.pi*Simulation.volume))
           print('ball', ball1_id, 'collides with container')
        
    def run(self, num_frames, animate=False):
        if animate:
            ax = pl.axes(xlim=(-Simulation.volume, Simulation.volume), ylim=(-Simulation.volume, Simulation.volume))
            ax.add_artist(self._container.get_patch())
            p = []
            for i in range(self._NoB):
                p.append(ax.add_patch(Simulation.ball_list[i].get_patch()))
        dist_list = []
        centr_list = []
        tot_energy_list = []
        tot_momentum_list = []
        velocities = []
        for frame in range(num_frames):
            energy_list = []
            momenta = []
            for i in range(self._NoB):
                centr_list.append(np.sqrt(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[i].pos())))
                energy_list.append((Simulation.ball_list[i]._mass*np.dot(Simulation.ball_list[i]._vel, Simulation.ball_list[i]._vel))/2)
                momenta.append(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel)
                velocities.append(Simulation.ball_list[i]._vel)
                for k in range(self._NoB):
                    if k != i:
                        dist_list.append(np.sqrt(abs(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[k].pos()))))
            tot_energy_list.append(sum(energy_list))
            tot_momentum_list.append(sum(momenta))
            self.next_collision()
            '''
            print(Energy)
            Pressure = (Energy)/(np.pi*(self._container._rad)**2)
            print(Pressure)'''
            if animate:
                    pl.pause(0.001)
        '''kB = 1
        T = (tot_energy_list*2)/(3*kB*NoB)
        bolts = []
        for i in range(NoB):
            for j in range(i, i+num_frames):
                bolts.append(velocities[j]*np.exp(-(0.5*(velocities[j]**2))/(kB*T[i])))'''
        hist1 = plt.figure()
        ax1 = hist1.add_subplot(1, 1, 1)
        ax1.hist(centr_list)       
        plt.title('Histogram of distances of the balls from the center')
        ax1.set_xlabel('Distance from the center')
        ax1.set_ylabel('Number of instances')
        plt.show()
        hist2 = plt.figure()
        ax1 = hist2.add_subplot(1, 1, 1)
        plt.hist(dist_list)
        plt.title('Histogram of distances between the balls')
        plt.xlabel('Distance from other balls')
        plt.ylabel('Number of instances')
        plt.show()
        energy = plt.figure()
        ax1 = energy.add_subplot(1, 1, 1)
        ax1.plot(tot_energy_list)       
        plt.title('Total kinetic energy of the gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Energy')
        plt.show()
        momentum = plt.figure()
        ax1 = momentum.add_subplot(1, 1, 1)
        ax1.plot(tot_momentum_list)       
        plt.title('Total momentum of the gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()
        pressure = plt.figure()
        ax1 = pressure.add_subplot(1, 1, 1)
        ax1.plot(Simulation.pressure_hits)       
        plt.title('Pressure on the container')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Pressure')
        plt.show()
        total_pressure = (sum(Simulation.pressure_hits))/num_frames
        print(total_pressure)
        '''
        boltzmann = plt.figure()
        ax1 = boltzmann.add_subplot(1, 1, 1)
        ax1.plot(bolts)       
        plt.title('Total momentum of the gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()'''
        
I = Simulation()
I.run(100, True)

#%%
for i in range(self._NoB): 
    np.histogram(np.sqrt(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[i].pos())))

