# -*- coding: utf-8 -*-
"""
Adam Brzosko
Ball class for Y2 Proj
"""
import numpy as np
import pylab as pl
import random as rand
import math

class Ball:
    ballNo = 0
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False):
        self._mass = m
        self._rad = R
        self._pos = np.array(r)
        self._vel = np.array(v)
        Ball.ballNo +=1
        col1 = (rand.randint(0,10))/10
        col2 = (rand.randint(0,10))/10
        col3 = (rand.randint(0,10))/10
        col = tuple([col1,col2,col3])
        self._patch = pl.Circle(self._pos, self._rad, color=col)
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
            t = ((-b)+np.sqrt(det))/(2*a)   #not comparing to zero cause of floating point accuracy
            if t > 1e-14:           #choose one order of magnitude above floating point accuracy
                return t
            else: return float('inf')
        elif det >= 0 and b < 0:
            t =((-b)-np.sqrt(det))/(2*a)
            if t > 1e-14:
                return t
            else: return float('inf')
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
    volume = 20 #this is actually the radius
    pressure_hits = []
    T = 0
    ball1_id = 0  #collision ball ids
    ball2_id = 0
    Container_collision = True
    tc = []  #container collisions list
    tb = []  #ball collisions list
    cb = []  #collision ball ids
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 50  #maximum is (f^2)/(3^2), which is just a smaller square than half of square of radius / (3)^2 
        self._container = Ball(container=True)
        f = math.floor(np.sqrt(0.5*(Simulation.volume**2)))-1
        size = range(-f,f,3*R)  #sqrt of half of the square of (radius minus one) (size of possible grid)
        vx = []
        for i in size:       #velocities are random, but container size dependent
            vx.append(i)
            vy = vx.copy()
        a = []
        b = []
        for i in size:      #generate grid of all possible spawn locations
            a.append(i)
        for j in size:
            c = a.copy()
            c.append(j/(3*R))
            b.append(c)
        for k in range(self._NoB):
            x = rand.choice(b) 
            while len(x) == 1:
                x = rand.choice(b)
            pos_x = b.index(x)
            x.pop()
            pos_y = rand.choice(x)
            x.remove(pos_y)
            x.append(pos_x)
            Simulation.ball_list.append(Ball(r= np.array([(3*R*pos_x)-f,pos_y]), v=np.array([rand.choice(vx),rand.choice(vy)])))
    def next_collision_initial(self):
        for k in range(self._NoB):
            tck = Simulation.ball_list[k].time_to_collision(self._container)
            Simulation.tc.append(tck)
            tbk = []
            for i in range(self._NoB):
                time = Simulation.ball_list[k].time_to_collision(Simulation.ball_list[i])
                tbk.append(time)
            tbk_min = min(tbk) #shortest time to another ball
            Simulation.tb.append(tbk_min)
            Simulation.cb.append(tbk.index(tbk_min))
        t_cont_min = min(Simulation.tc)
        t_ball_min = min(Simulation.tb)
        if t_cont_min < t_ball_min:
            t = t_cont_min
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont_min)
        else:
            t = t_ball_min
            Simulation.Container_collision = False
            Simulation.ball1_id = Simulation.tb.index(t_ball_min)
            Simulation.ball2_id = Simulation.cb[Simulation.ball1_id]
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t) #move all the balls
            Simulation.tb[i] = Simulation.tb[i] - t #increment all the times
            Simulation.tc[i] = Simulation.tc[i] - t
        if Simulation.Container_collision == False:
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
            print('ball', Simulation.ball1_id, 'collides with', Simulation.ball2_id)
            Simulation.pressure_hits.append(0)
        else:
           Simulation.ball_list[Simulation.ball1_id].collide(self._container)
           Simulation.pressure_hits.append((2*Simulation.ball_list[Simulation.ball1_id]._mass*np.sqrt(abs(np.dot(Simulation.ball_list[Simulation.ball1_id]._vel,Simulation.ball_list[Simulation.ball1_id]._vel))))/(2*np.pi*Simulation.volume))
           print('ball', Simulation.ball1_id, 'collides with container')    

    def next_collision_main(self):
        tb1 = []
        tb2 = []
        if Simulation.Container_collision == True:
            for i in range(self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].time_to_collision(Simulation.ball_list[i])) #calculate shortest collision time with ball
            t = min(tb1)
            Simulation.tb[Simulation.ball1_id] = t   #replace shortest collision time with ball
            Simulation.cb[Simulation.ball1_id] = tb1.index(t) #which ball
            Simulation.tc[Simulation.ball1_id] = Simulation.ball_list[Simulation.ball1_id].time_to_collision(self._container) #replace time to collision with container
        else:
            for i in range (self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].time_to_collision(Simulation.ball_list[i])) #calculate shortest collision times with balls
                tb2.append(Simulation.ball_list[Simulation.ball2_id].time_to_collision(Simulation.ball_list[i]))
            t1 = min(tb1)
            t2 = min(tb2)
            Simulation.tb[Simulation.ball1_id] = t1
            Simulation.cb[Simulation.ball1_id] = tb1.index(t1)
            Simulation.tb[Simulation.ball2_id] = t2
            Simulation.cb[Simulation.ball2_id] = tb2.index(t2)
            Simulation.tc[Simulation.ball1_id] = Simulation.ball_list[Simulation.ball1_id].time_to_collision(self._container)
            Simulation.tc[Simulation.ball2_id] = Simulation.ball_list[Simulation.ball2_id].time_to_collision(self._container)
        t_ball = min(Simulation.tb)  #choose shortest times from the general list
        t_cont = min(Simulation.tc)   
        if t_ball < t_cont:
            Simulation.Container_collision = False
            Simulation.ball1_id = Simulation.tb.index(t_ball)   #get coliding ball numbers
            Simulation.ball2_id = Simulation.cb[Simulation.ball1_id]
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_ball) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_ball #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_ball
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
            print('ball', Simulation.ball1_id, 'collides with', Simulation.ball2_id)
        else:
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont)
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_cont) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_cont #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_cont
            Simulation.ball_list[Simulation.ball1_id].collide(self._container)
            print('ball', Simulation.ball1_id, 'collides with container')  
            
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
        self.next_collision_initial()
        pl.pause(0.000001)
        for frame in range(num_frames-1):
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
            self.next_collision_main()

            if animate:
                    pl.pause(0.000001)
