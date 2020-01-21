# -*- coding: utf-8 -*-
"""
Adam Brzosko

Ball class for Y2 Proj
"""
import numpy as np
import pylab as pl

class Ball:
    ballNo = 0
    
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), Patch = True):
        self._mass = m
        self._rad = R
        self._pos = np.array(r)
        self._vel = np.array(v)
        Ball.ballNo +=1
        if Patch == True:
            patch = pl.Circle(r, R, fc='r')
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_patch(patch)
        
    def __repr__(self):
        return 'Ball of mass %s, radius %s, velocity %s, at %s' % (self._mass, self._rad, self._vel, self._pos)
    
    def pos(self):
        return self._pos
    
    def vel(self):
        return self._vel
    
    def move(self,dt):
        self._pos = self._pos + dt*self._vel
        return self._pos
        
    def time_to_collision(self, other):
        r = self._pos - other._pos
        v = self._vel - other._vel 
        R = self._rad + other._rad
        a = np.dot(v,v)
        b = 2*np.dot(r,v)
        c = (np.dot(r,r) - (R**2))
        if other._mass == float('inf'):
            return ((-b)+np.sqrt((b**2)-(4*a*c)))/(2*a)
        elif b**2 - (4*a*c) == 0: 
            return (-b)/(2*a)
        elif b**2 - (4*a*c) > 0 and b < 0 :
            return  ((-b)-np.sqrt((b**2)-(4*a*c)))/(2*a)
        else: raise Exception ("No collision")
    
    def collide(self, other):
        '''change vel components parallel to the axis
        containing centres of the balls in both balls'''
        pos_self_vector = other._pos - self._pos  #vectors along the axis
        pos_other_vector = - pos_self_vector          
      
        paral_vel_self = ((np.dot(self._vel,pos_self_vector))* \
(pos_self_vector))/(np.dot(pos_self_vector,pos_self_vector))   #component splitting self
        perp_vel_self = self._vel - paral_vel_self      
        
        if other._mass == float('inf'):
            paral_vel_other = 0
            perp_vel_other = 0
        else:
            paral_vel_other = ((np.dot(other._vel,pos_other_vector))* \
(pos_other_vector))/(np.dot(pos_other_vector,pos_other_vector))    #component splitting other
        perp_vel_other = other._vel - paral_vel_other
        
        if other._mass == float('inf'):
            new_paral_vel_self = - paral_vel_self
            new_paral_vel_other = 0
        else:
            new_paral_vel_self = - (((self._mass - other._mass)/ \
(self._mass + other._mass))*paral_vel_self) + \
(((2*other._mass)/(self._mass + other._mass))*paral_vel_other)  #calculate new component self
            new_paral_vel_other = - (((2*self._mass)/ \
(self._mass + other._mass))*paral_vel_self) + \
(((other._mass-self._mass)/(self._mass + other._mass))*paral_vel_other)   #calculate new component other
        
        new_vel_self = perp_vel_self + new_paral_vel_self
        new_vel_other = perp_vel_other + new_paral_vel_other
        self._vel = new_vel_self
        other._vel = new_vel_other      #change the componenets
        
        
class Simulation(Ball):
    
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), Patch = True):
        self._mass = m
        self._rad = R
        self._pos = np.array(r)
        self._vel = np.array(v)
        Ball.__init__(self,m,R,r,v)
        if Patch == True:
            patch = pl.Circle([0., 0.], 10, ec='b', fill=False, ls='solid')
            ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
            ax.add_patch(patch)
    
    def next_collision(self, other): 
        t = self.time_to_collision(other)
        self.move(t)
        other.move(t)
        self.collide(other)
