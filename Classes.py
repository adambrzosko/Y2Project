# -*- coding: utf-8 -*-
"""
Adam Brzosko
File containing classes (ball and simulation)
for Y2 Project
"""

import numpy as np
import pylab as pl
import random as rand
import math
import pickle

class Ball:
    '''A class for initialising Ball objects (particles and a container).
    
    Attributes
    ----------
    m: Mass of a the object
    R: Radius of the object
    r: 2D numpy array containing initial position of the object
    v: 2D numpy array containing initial velocity of the object
    container: (True/False) 'True' if container, 'False' if particle
    '''
    ballNo = 0
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]),\
                 container=False):
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
            self._patch = pl.Circle([0., 0.], Simulation.volume, ec='b',\
                                    fill=False, ls='solid')
        
    def __repr__(self):
        return 'Ball of mass %s, radius %s, velocity %s, at %s' %\
                (self._mass, self._rad, self._vel, self._pos)
    
    def mass(self):
        '''Returns the mass of a Ball object.'''
        return self._mass
    
    def pos(self):
        '''Returns the position of a Ball object.'''
        return self._pos
    
    def vel(self):
        '''Returns the velocity of a Ball object.'''
        return self._vel
    
    def move(self,dt):
        '''Moves a Ball object by time dt, according to its velocity.
        
        Arguments
        ---------
        dt: Time through which to move the object.
        '''
        self._pos = self._pos + dt*self._vel
        self._patch.center = self._pos
        return self._pos
        
    def time_to_collision(self, other):
        '''Calculates time to the next collision with a Ball object.
        
        Arguments
        ---------
        other: Object, which a collision time to is calculated.
        '''
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
            t = ((-b)+np.sqrt(det))/(2*a)   
            if t > 1e-14:
                return t
            else: return float('inf')
        elif det >= 0 and b < 0:
            t =((-b)-np.sqrt(det))/(2*a)
            if t > 1e-14:
                return t
            else: return float('inf')
        else: return float('inf')
    
    def collide(self, other):
        '''Collides the object with another Ball object.
        
        Arguments
        ---------
        other: Object which to collide with
        '''
        pos_self_vector = other.pos()-self.pos()  #relative position vectors
        pos_other_vector = -pos_self_vector          
      
        paral_vel_self = ((np.dot(self.vel(),pos_self_vector))*\
                          (pos_self_vector))/\
                          (np.dot(pos_self_vector,pos_self_vector))
        perp_vel_self = self.vel()-paral_vel_self      #splitting components
        
        if other._mass == float('inf'): #collision with a container
            new_paral_vel_other = 0
            perp_vel_other = 0
            new_paral_vel_self = -paral_vel_self
        else:
            paral_vel_other = ((np.dot(other.vel(),pos_other_vector))*\
                               (pos_other_vector))/\
                               (np.dot(pos_other_vector,pos_other_vector))
            perp_vel_other = other._vel-paral_vel_other
            new_paral_vel_self = (((self.mass()-other.mass())/\
                                 (self.mass()+other.mass()))*paral_vel_self)+\
                                 (((2*other.mass())/(self.mass()\
                                 +other.mass()))*paral_vel_other)
            new_paral_vel_other = (((2*self.mass())/(self.mass()+other.mass()))\
                                  *paral_vel_self)+(((other.mass()-self.mass())\
                                                /(self.mass()+other.mass()))\
                                                *paral_vel_other)   
        self._vel = perp_vel_self+new_paral_vel_self #recombination of comps
        other._vel = perp_vel_other+new_paral_vel_other    
        
    def get_patch(self):
        '''Returns a patch attribute (draws a circular patch).'''
        return self._patch


class Simulation:
    '''A class for running a simulation of a number of Ball objects.
    
    Attributes
    ----------
    m: Mass of the object
    R: Radius of the object
    r: 2D numpy array containing initial position of the object
    v: 2D numpy array containing initial velocity of the object
    container: (True/False) True if container, False if particle
    NoB: Number of objects in the simulation
    '''
    ball_list = []
    simulation_time = []
    volume = 20 #this is actually the radius of the container
    pressure_hits = []
    ball1_id = 0  #collision ball ids
    ball2_id = 0
    Container_collision = True
    tc = []  #container collisions list
    tb = []  #ball collisions list
    cb = []  #next collision ball ids list
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]),\
                 container=False, NoB=1):
        self._NoB = NoB  #maximum is ((2f)^2)/(3^2) 
        self._container = Ball(container=True)
        f = math.floor(np.sqrt(0.5*(Simulation.volume**2)))-1
        size = range(-f,f,3*R)  #(size of possible grid)
        vx = []
        for i in size:    #velocities are random, but container size dependent
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
            Simulation.ball_list.append(Ball(r=np.array([(3*R*pos_x)-f,pos_y]),\
                                             v=np.array([rand.choice(vx),\
                                                         rand.choice(vy)])))
    def next_collision_initial(self):
        '''Performs the first collision, including calculation of the shortest
        collision times for each ball.
        '''
        for k in range(self._NoB):
            tck = Simulation.ball_list[k].time_to_collision(self._container)
            Simulation.tc.append(tck)
            tbk = []
            for i in range(self._NoB):
                time = Simulation.ball_list[k].\
                        time_to_collision(Simulation.ball_list[i])
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
        Simulation.simulation_time.append(t)
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t) #move all the balls
            Simulation.tb[i] = Simulation.tb[i] - t #increment all the times
            Simulation.tc[i] = Simulation.tc[i] - t
        if Simulation.Container_collision == False:
            Simulation.ball_list[Simulation.ball1_id].\
                collide(Simulation.ball_list[Simulation.ball2_id])
            print('Particle', Simulation.ball1_id, 'collides with',\
                  Simulation.ball2_id)
            Simulation.pressure_hits.append(0)
        else:
           Simulation.ball_list[Simulation.ball1_id].collide(self._container)
           Simulation.pressure_hits.append((\
            Simulation.ball_list[Simulation.ball1_id].mass()*\
            np.sqrt(abs(np.dot(Simulation.ball_list[Simulation.ball1_id].vel(),\
            Simulation.ball_list[Simulation.ball1_id].vel()))))/\
            (np.pi*Simulation.volume))
           print('Particle', Simulation.ball1_id, 'collides with container')    

    def next_collision_main(self):
        '''Performs all collisions subsequent to the first, along with
        calculation of shortest collision time for particles which just 
        collided.
        '''
        tb1 = []
        tb2 = []
        if Simulation.Container_collision == True:
            for i in range(self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].\
                           time_to_collision(Simulation.ball_list[i])) 
                #calculate shortest collision time with ball
            t = min(tb1)
            Simulation.tb[Simulation.ball1_id] = t   
            #replace shortest collision time with ball
            Simulation.cb[Simulation.ball1_id] = tb1.index(t) 
            #which ball
            Simulation.tc[Simulation.ball1_id] = Simulation.\
              ball_list[Simulation.ball1_id].time_to_collision(self._container) 
              #replace time to collision with container
        else:
            for i in range (self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].\
                           time_to_collision(Simulation.ball_list[i])) 
                #calculate shortest collision times with balls
                tb2.append(Simulation.ball_list[Simulation.ball2_id].\
                           time_to_collision(Simulation.ball_list[i]))
            t1 = min(tb1)
            t2 = min(tb2)
            Simulation.tb[Simulation.ball1_id] = t1
            Simulation.cb[Simulation.ball1_id] = tb1.index(t1)
            Simulation.tb[Simulation.ball2_id] = t2
            Simulation.cb[Simulation.ball2_id] = tb2.index(t2)
            Simulation.tc[Simulation.ball1_id] = Simulation.\
              ball_list[Simulation.ball1_id].time_to_collision(self._container)
            Simulation.tc[Simulation.ball2_id] = Simulation.\
              ball_list[Simulation.ball2_id].time_to_collision(self._container)
        t_ball = min(Simulation.tb)#choose shortest times from the general list
        t_cont = min(Simulation.tc)   
        if t_ball < t_cont:
            Simulation.Container_collision = False
            Simulation.pressure_hits.append(0)
            Simulation.ball1_id = Simulation.tb.index(t_ball)
            #get coliding ball numbers
            Simulation.ball2_id = Simulation.cb[Simulation.ball1_id]
            Simulation.simulation_time.append(t_ball)
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_ball) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_ball 
                #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_ball
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.\
                                ball_list[Simulation.ball2_id])
            print('Particle', Simulation.ball1_id, 'collides with',\
                  Simulation.ball2_id)
        else:
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont)
            Simulation.simulation_time.append(t_cont)
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_cont) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_cont 
                #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_cont
            Simulation.pressure_hits.append((Simulation.\
                ball_list[Simulation.ball1_id].mass()*np.sqrt(abs(np.\
                dot(Simulation.ball_list[Simulation.ball1_id].vel(),Simulation.\
                ball_list[Simulation.ball1_id].vel()))))/\
                (np.pi*Simulation.volume))
            Simulation.ball_list[Simulation.ball1_id].collide(self._container)
            print('Particle', Simulation.ball1_id, 'collides with container')  
            
    def run(self, num_frames, animate=False):
        '''Executes the simulation and saves data to external files. 
        
        Arguments
        ---------
        num_frames: Specifies how many collisions to perform
        animate: (True/False) If 'True' displays the animation
        '''
        if self._NoB > 80:
            raise Exception('Too many particles')
        if animate:
            ax = pl.axes(xlim=(-Simulation.volume, Simulation.volume),\
                         ylim=(-Simulation.volume, Simulation.volume))
            ax.add_artist(self._container.get_patch())
            p = []
            for i in range(self._NoB):
                p.append(ax.add_patch(Simulation.ball_list[i].get_patch()))
        velocities = []
        masses = []
        self.next_collision_initial()
        for frame in range(num_frames-1):
            for i in range(self._NoB):
                masses.append(Simulation.ball_list[i].mass())
                velocities.append(Simulation.ball_list[i].vel())
            self.next_collision_main()
            if animate:
                    pl.pause(0.000001)

        with open('masses.data', 'wb') as filehandle:
            pickle.dump(masses, filehandle)
        with open('velocities.data', 'wb') as filehandle:
            pickle.dump(velocities, filehandle)
        with open('instantenous_pressure.data', 'wb') as filehandle:
            pickle.dump(Simulation.pressure_hits, filehandle)
        with open('time.data', 'wb') as filehandle:
            pickle.dump(Simulation.simulation_time, filehandle)


        