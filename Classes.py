# -*- coding: utf-8 -*-
"""
Adam Brzosko
Ball class for Y2 Proj
"""
import numpy as np
import pylab as pl
import random as rand

class Ball:
    ballNo = 0
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False):
        self._mass = m
        self._rad = R
        self._pos = np.array(r)
        self._vel = np.array(v)
        Ball.ballNo +=1
        col1 = rand.randint(0,255)
        col2 = rand.randint(0,255)
        col3 = rand.randint(0,255)
        col = tuple([col1,col2,col3])
        self._patch = pl.Circle(self._pos, self._rad, color=str(col1))
        if container == True:
            Ball.ballNo -=1
            self._mass = float('inf')
            self._rad = -Simulation.volume
            self._pos = np.array([0,0])
            self._vel = np.array([0,0])
            self._patch = pl.Circle([0., 0.], Simulation.volume, ec='b', fill=False, ls='solid')
        
    def __repr__(self):
        return 'Ball of mass %s, radius %s, velocity %s, at %s' % (self._mass, self._rad, self._vel, self._pos)
    
    def mass(self):
        return self._mass
    
    def pos(self):
        return self._pos
    
    def vel(self):
        return self._vel
    
    def move(self,dt):
        self._pos = self._pos + dt*self._vel# - 1e-14
        self._patch.center = self._pos
        return self._pos
        
    def time_to_collision(self, other):
        r = self._pos - other._pos
        v = self._vel - other._vel 
        R = self._rad + other._rad
        a = np.dot(v,v)
        b = 2*np.dot(r,v)
        c = (np.dot(r,r) - (R**2))
        det = b**2 - (4*a*c)
        if a == 0: 
            return float('inf')
        elif other._mass == float('inf'):
            return ((-b)+np.sqrt(det))/(2*a)   #not comparing to zero cause of floating point accuracy
        elif det == 0 and b < 0: 
            return (-b)/(2*a)
        elif det > 0 and b < 0:
            return  ((-b)-np.sqrt(det))/(2*a)
        else: return float('inf')
    
    def collide(self, other):
        '''change vel components parallel to the axis
        containing centres of the balls in both balls'''
        pos_self_vector = other._pos-self._pos  #vectors along the axis
        pos_other_vector = -pos_self_vector          
      
        paral_vel_self = ((np.dot(self._vel,pos_self_vector))*(pos_self_vector))/(np.dot(pos_self_vector,pos_self_vector))   #component splitting self
        perp_vel_self = self._vel-paral_vel_self      
        
        if other._mass == float('inf'):
            new_paral_vel_other = 0
            perp_vel_other = 0
            new_paral_vel_self = -paral_vel_self
        else:
            paral_vel_other = ((np.dot(other._vel,pos_other_vector))*(pos_other_vector))/(np.dot(pos_other_vector,pos_other_vector))    #component splitting other
            perp_vel_other = other._vel-paral_vel_other
            new_paral_vel_self = (((self._mass-other._mass)/(self._mass+other._mass))*paral_vel_self)+(((2*other._mass)/(self._mass+other._mass))*paral_vel_other)  #calculate new component self
            new_paral_vel_other = (((2*self._mass)/(self._mass+other._mass))*paral_vel_self)+(((other._mass-self._mass)/(self._mass+other._mass))*paral_vel_other)   #calculate new component other
        self._vel = perp_vel_self+new_paral_vel_self
        other._vel = perp_vel_other+new_paral_vel_other    #change the componenets
        
    def get_patch(self):
        return self._patch


class Simulation:
    ball_list = []
    volume = 10
    
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 3
        self._container = Ball(container=True)
        x = []
        for i in range(-10,11,2):    #sqrt of half of the square of radius
            x.append(i)
            y = x.copy()
            vx = x.copy()
            vy = x.copy()
        for i in range(self._NoB):
            pos_x = rand.choice(x)
            pos_y = rand.choice(y)
            Simulation.ball_list.append(Ball(r= np.array([pos_x,pos_y]), v=np.array([rand.choice(vx),rand.choice(vy)])))
            x.remove(pos_x)
            y.remove(pos_y)

    def next_collision(self):
        tb = []
        t = Simulation.ball_list[self._NoB-1].time_to_collision(self._container)
        ball1_id = self._NoB-1    #set last ball as default to collide with a container
        Container_collision = True
        for k in range(self._NoB-1):     #calculate all other balls times to collisions
            tc = Simulation.ball_list[k].time_to_collision(self._container)  #first collision time with container
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
            tb = []
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t)
        if Container_collision == False:
            Simulation.ball_list[ball1_id].collide(Simulation.ball_list[ball2_id])
        else:
           Simulation.ball_list[ball1_id].collide(self._container)
        
    def run(self, num_frames, animate=False):
        if animate:
            ax = pl.axes(xlim=(-Simulation.volume, Simulation.volume), ylim=(-Simulation.volume, Simulation.volume))
            ax.add_artist(self._container.get_patch())
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
                    pl.pause(0.5)

