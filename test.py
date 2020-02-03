#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 20:57:14 2020

@author: adam
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



#%%
import numpy as np
import pylab as pl
import random as rand
from classes import Ball

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
                if tb != [] and tb[i-k-1] < t and tb[i-k-1] < tc:
                    t = tb[i-k-1]
                    ball1_id = k
                    ball2_id = i
                    Container_collision = False
                    print(k,i)
                elif tc < t:
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
I.run(3, True)

#%%
for i in range(self._NoB):
    print(Simulation.ball_list[i].pos())