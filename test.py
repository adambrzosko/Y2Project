"""
Created on Tue Jan 28 20:57:14 2020
@author: Adam Brzosko
Y2 project
This is a test file, containing snippets of code that were
used to build andtest the code used in the main part 
"""
#%%
'''This cell runs the simulation thay was last loaded
with a given number of frames. True specifies that
 the animation is shown.'''
from Classes import Simulation as s

I = s()

I.run(100, True)

#%%
'''Iniitial idea for spawning the balls randomly.
Two identical lists are created. The values are 
chosen in pairs from both lists at random and then
removed.'''
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
'''Creating the method for collliding balls.'''
from Classes import Ball
import numpy as np
def collide(self, other):
    '''change vel components parallel to the axis
    containing centres of the balls in both balls'''
    pos_self_vector = other._pos - self._pos  #vectors along the axis
    pos_other_vector = - pos_self_vector       
    
    paral_vel_self = ((np.dot(self._vel,pos_self_vector))*(pos_self_vector))/(np.dot(pos_self_vector,pos_self_vector))   #component splitting self
    perp_vel_self = self._vel - paral_vel_self      
        
    if other._mass == float('inf'):   #container collision case 
        new_paral_vel_other = 0   #only reverse the colliding ball's velocity
        perp_vel_other = 0          
        new_paral_vel_self = - paral_vel_self
    else:
        paral_vel_other = ((np.dot(other._vel,pos_other_vector))*(pos_other_vector))/(np.dot(pos_other_vector,pos_other_vector))    #component splitting other
        perp_vel_other = other._vel - paral_vel_other
        new_paral_vel_self = (((self._mass - other._mass)/(self._mass + other._mass))*paral_vel_self) + (((2*other._mass)/(self._mass + other._mass))*paral_vel_other)  #calculate new component self
        new_paral_vel_other = (((2*self._mass)/(self._mass + other._mass))*paral_vel_self) + (((other._mass-self._mass)/(self._mass + other._mass))*paral_vel_other)   #calculate new component other
    self._vel = perp_vel_self + new_paral_vel_self    #change the componenets
    other._vel = perp_vel_other + new_paral_vel_other    


b1 = Ball(r=np.array([4,4]), v=np.array([-1,-1]))  #test on two balls
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
'''This cell test a simulation of a number of balls colliding,
including the animation.'''
import numpy as np
import pylab as pl
import random as rand
from Classes import Ball

class Simulation:
    ball_list = [] #all the balls apart from container are stored in this list
    volume = 20  #this is the radius of the container
    
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
            if animate:
                for i in range(self._NoB):
                    print(p[i])
                    p[i].center = (Simulation.ball_list[i]._pos)
                    print(p[i])
                    pl.pause(0.001)

I = Simulation()
I.run(20, True)
#%%
'''Improved simulation class. The initilisaion method is more efficient.
The run method now fixed. Added plots for kinetic energy,
momentum, distances from centre and other balls and comparison
to maxwell-boltzmann distribution.'''
import matplotlib.pyplot as plt

class Simulation:
    ball_list = []
    volume = 20
    pressure_hits = []
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 11
        self._container = Ball(container=True)
        x = []
        for i in range(-10,11,2):    
            x.append(i)
        y = x.copy()
        vx = x.copy()
        vy = x.copy()
        for i in range(self._NoB):
            pos_x = rand.choice(x)
            pos_y = rand.choice(y)
            Simulation.ball_list.append(Ball(r= np.array([pos_x,pos_y]), v=np.array([2*rand.choice(vx),2*rand.choice(vy)])))
            x.remove(pos_x)
            y.remove(pos_y)

    def next_collision(self):
        tb = []
        t = Simulation.ball_list[self._NoB-1].time_to_collision(self._container)
        print('t is', t)
        ball1_id = self._NoB-1    
        Container_collision = True
        for k in range(self._NoB-1):     
            tc = Simulation.ball_list[k].time_to_collision(self._container)  
            for i in range(k+1, self._NoB): 
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
        tot_mx_list = []
        tot_my_list = []
        vel_squares = []
        for frame in range(num_frames):
            energy_list = []
            mx = []
            my = []
            for i in range(self._NoB):
                centr_list.append(np.sqrt(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[i].pos())))
                energy_list.append((Simulation.ball_list[i]._mass*np.dot(Simulation.ball_list[i]._vel, Simulation.ball_list[i]._vel))/2)
                mx.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[0]))
                my.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[1]))
                vel_squares.append(abs(np.dot(Simulation.ball_list[i]._vel,Simulation.ball_list[i]._vel)))
                for k in range(self._NoB):
                    if k != i:
                        dist_list.append(np.sqrt(abs(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[k].pos()))))
            tot_energy_list.append(sum(energy_list))
            tot_mx_list.append(sum(mx))
            tot_my_list.append(sum(my))
            tot_m_list = [a-b for a,b in zip(tot_mx_list,tot_mx_list)]
            self.next_collision()
            if animate:
                    pl.pause(0.001)
        
        kB = 1
        T = [(x*2)/(3*kB*self._NoB) for x in tot_energy_list]
        average_pressure = (sum(Simulation.pressure_hits))/num_frames
        print(average_pressure)

        bolts = []
        for i in range(num_frames):
            for j in range(i*self._NoB, (i+1)*self._NoB):
                bolts.append(vel_squares[j]*np.exp(-(0.5*(vel_squares[j]**2))/(kB*T[i])))
       
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
        plt.title('Total kinetic energy of the gas molecules')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Energy')
        plt.show()
        
        temperature = plt.figure()
        ax1 = temperature.add_subplot(1, 1, 1)
        ax1.plot(T)       
        plt.title('Temperature of the gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Temperature')
        plt.show()
        
        momentumx = plt.figure()
        ax1 = momentumx.add_subplot()
        ax1.plot(tot_mx_list)       
        plt.title('Maxwell-Boltzmann distribution for our gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()
        
        momentumy = plt.figure()
        ax1 = momentumy.add_subplot(1, 1, 1)
        ax1.plot(tot_my_list)       
        plt.title('Total momentum of the gas molecules in y direction')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()
        
        momentum = plt.figure()
        ax1 = momentum.add_subplot(1, 1, 1)
        ax1.plot(tot_m_list)       
        plt.title('Total momentum of the gas molecules')
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
        
        boltzmann = plt.figure()
        ax1 = boltzmann.add_subplot(1, 1, 1)
        plt.hist(bolts)       
        plt.title('Total momentum of the gas')
        ax1.set_xlabel('Energy of a molecule')
        ax1.set_ylabel('Number of instances')
        plt.show() 
        
I = Simulation()
I.run(100, True)

#%%
'''An attempt at optimisng the calculation of times to collision.
The main idea is calculating the collision times as previously
for firt frame and then after every collision only the ones for
the balls that just collided and negatively incrementing other
times by the smallest time each time. This significantly reduces
the complexity of calculation.'''
class Simulation:
    ball_list = []
    volume = 20
    pressure_hits = []
    T = 0
    ball1_id = 0  #collision ball ids
    ball2_id = 0
    Container_collision = True
    tb = []  #ball collision times
    tc = []  #container collision times
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 11
        self._container = Ball(container=True)
        x = []
        for i in range(-10,11,2):    
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
            
    def next_collision_initial(self):
        t = Simulation.ball_list[self._NoB-1].time_to_collision(self._container)
        print('t is', t)
        Simulation.tc.append(t) #append this time to the container collision times
        Simulation.ball1_id = self._NoB-1    #set last ball as default to collide with a container
        Simulation.Container_collision = True
        for k in range(self._NoB-1):     #calculate all other balls times to collisions
            tck = Simulation.ball_list[k].time_to_collision(self._container)  #first calculate collision time with container
            Simulation.tc.append(tck) #append it
            print('tck is:', tck)
            for i in range(k+1, self._NoB):  #then with all balls after it on the list
                Simulation.tb.append(Simulation.ball_list[k].time_to_collision(Simulation.ball_list[i])) #append them
                if Simulation.tb != [] and Simulation.tb[i-k-1] < t and Simulation.tb[i-k-1] < tck and Simulation.tb[i-k-1] > 1e-15:
                    t = Simulation.tb[i-k-1]
                    Simulation.ball1_id = k
                    Simulation.ball2_id = i
                    Simulation.Container_collision = False
                    print(k,i)
                elif tck < t and tck > 1e-15:
                    t = tck
                    Simulation.ball1_id = k
                    Simulation.Container_collision = True
            print('tb is:', Simulation.tb)
        for i in range(self._NoB):
            Simulation.ball_list[i].move(t) #move all the balls
            Simulation.tb[i] = Simulation.tb[i] - t #increment all the times
            Simulation.tc[i] = Simulation.tc[i] - t
        Simulation.tc.append(Simulation.tc.pop(0))  #move the first entry to the end to obtain an ordered list
        if Simulation.Container_collision == False:
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
            print('ball', Simulation.ball1_id, 'collides with', Simulation.ball2_id)
            Simulation.pressure_hits.append(0)
        else:
           Simulation.ball_list[Simulation.ball1_id].collide(self._container)
           Simulation.pressure_hits.append((2*Simulation.ball_list[Simulation.ball1_id]._mass*np.sqrt(abs(np.dot(Simulation.ball_list[Simulation.ball1_id]._vel,Simulation.ball_list[Simulation.ball1_id]._vel))))/(2*np.pi*Simulation.volume))
           print('ball', Simulation.ball1_id, 'collides with container')    
           print(Simulation.tc)
           print(Simulation.tb)
    def next_collision_main(self):
        tb1 = []
        tb2 = []
        if Simulation.Container_collision == True:
            for i in range(self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].time_to_collision(Simulation.ball_list[i]))
            t1 = min(tb1)
            t1c = Simulation.ball_list[Simulation.ball1_id].time_to_collision(self._container)
            if t1 > 1e-15:
                Simulation.tb[Simulation.ball1_id] = min(tb1)
            if t1c > 1e-15:
                Simulation.tc[Simulation.ball1_id] = t1c
        else:
            for i in range (self._NoB):
                tb1.append(Simulation.ball_list[Simulation.ball1_id].time_to_collision(Simulation.ball_list[i]))
                tb2.append(Simulation.ball_list[Simulation.ball2_id].time_to_collision(Simulation.ball_list[i]))
            t1 = min(tb1)
            t2 = min(tb2)
            if t1 > 1e-15:
                Simulation.tb[Simulation.ball1_id] = min(tb1)
            if t2 > 1e-15:
                Simulation.tb[Simulation.ball2_id] = min(tb2)
            t1c = Simulation.ball_list[Simulation.ball1_id].time_to_collision(self._container)
            t2c = Simulation.ball_list[Simulation.ball2_id].time_to_collision(self._container)
            if t1c > 1e-15:
                Simulation.tc[Simulation.ball1_id] = t1c
            if t2c > 1e-15:
                Simulation.tc[Simulation.ball2_id] = t2c
        t_ball = min(Simulation.tb)
        t_cont = min(Simulation.tc)
        k = -1
        n = 0
        if t_ball < t_cont:
            Simulation.Container_collision = False
            while n < Simulation.tb.index(t_ball):
                n = n + (self._NoB-(k+1))
                k += 1
            Simulation.ball1_id = k
            Simulation.ball2_id = self._NoB-(n-Simulation.tb.index(t_ball))
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_ball) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_ball #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_ball
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
        else:
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont)
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_cont) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_cont #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_cont
            Simulation.ball_list[Simulation.ball1_id].collide(self._container)
            
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
                    pl.pause(0.001)
        
I = Simulation()
I.run(100, True)
#%%
'''Now properly implemented the time calculation optimisation.'''
import math
class Simulation:
    ball_list = []
    volume = 20
    pressure_hits = []
    T = 0
    ball1_id = 0  #collision ball ids
    ball2_id = 0
    Container_collision = True
    Container_collision = True
    tc = []  #container collisions list
    tb = []  #ball collisions list
    cb = []  #collision ball ids
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 40  #maximum is (f^2)/(3^2), which is just a smaller square than half of square of radius / (3)^2 
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
        print('cb is', Simulation.cb)
        print('tc is', Simulation.tc)
        print('tb is', Simulation.tb)
        t_cont_min = min(Simulation.tc)
        t_ball_min = min(Simulation.tb)
        print('t_cont_min is', t_cont_min)
        print('t_ball_min is', t_ball_min)
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
        else:
           Simulation.ball_list[Simulation.ball1_id].collide(self._container)
           print('ball', Simulation.ball1_id, 'collides with container')    
           print(Simulation.tc)
           print(Simulation.tb)
           
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
            Simulation.pressure_hits.append(0)
            Simulation.ball1_id = Simulation.tb.index(t_ball)   #get coliding ball numbers
            Simulation.ball2_id = Simulation.cb[Simulation.ball1_id]
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_ball) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_ball #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_ball
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
        else:
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont)
            Simulation.pressure_hits.append((2*Simulation.ball_list[Simulation.ball1_id]._mass*np.sqrt(abs(np.dot(Simulation.ball_list[Simulation.ball1_id]._vel,Simulation.ball_list[Simulation.ball1_id]._vel))))/(2*np.pi*Simulation.volume))
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_cont) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_cont #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_cont
            Simulation.ball_list[Simulation.ball1_id].collide(self._container)
            
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
        tot_mx_list = []
        tot_my_list = []
        vel_squares = []
        self.next_collision_initial()
        pl.pause(0.001)
        for frame in range(num_frames-1):
            energy_list = []
            mx = []
            my = []
            for i in range(self._NoB):
                centr_list.append(np.sqrt(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[i].pos())))
                energy_list.append((Simulation.ball_list[i]._mass*np.dot(Simulation.ball_list[i]._vel, Simulation.ball_list[i]._vel))/2)
                mx.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[0]))
                my.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[1]))
                vel_squares.append(abs(np.dot(Simulation.ball_list[i]._vel,Simulation.ball_list[i]._vel)))
                for k in range(self._NoB):
                    if k != i:
                        dist_list.append(np.sqrt(abs(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[k].pos()))))
            tot_energy_list.append(sum(energy_list))
            tot_mx_list.append(sum(mx))
            tot_my_list.append(sum(my))
            tot_m_list = [a-b for a,b in zip(tot_mx_list,tot_mx_list)]
            self.next_collision_main()

            if animate:
                    pl.pause(0.001)
            
        kB = 1
        T = [(x*2)/(3*kB*self._NoB) for x in tot_energy_list]
        average_pressure = (sum(Simulation.pressure_hits))/num_frames
        print(average_pressure)

        bolts = []
        for i in range(num_frames-1):
            for j in range(i*self._NoB, (i+1)*self._NoB):
                bolts.append(vel_squares[j]*np.exp(-(0.5*(vel_squares[j]**2))/(kB*T[i])))
         
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
        plt.title('Total kinetic energy of the gas molecules')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Energy')
        plt.show()
        
        temperature = plt.figure()
        ax1 = temperature.add_subplot(1, 1, 1)
        ax1.plot(T)       
        plt.title('Temperature of the gas')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Temperature')
        plt.show()
        
        momentumx = plt.figure()
        ax1 = momentumx.add_subplot()    
        ax1.plot(tot_mx_list)       
        plt.title('Total momentum of the gas molecules in x direction')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()
        
        momentumy = plt.figure()
        ax1 = momentumy.add_subplot(1, 1, 1)
        ax1.plot(tot_my_list)       
        plt.title('Total momentum of the gas molecules in y direction')
        ax1.set_xlabel('Number of collisions from the beginning')
        ax1.set_ylabel('Momentum')
        plt.show()
        
        momentum = plt.figure()
        ax1 = momentum.add_subplot(1, 1, 1)
        ax1.plot(tot_m_list)       
        plt.title('Total momentum of the gas molecules')
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
        
        boltzmann = plt.figure()
        ax1 = boltzmann.add_subplot(1, 1, 1)
        plt.hist(bolts, bins = 50)       
        plt.title('Maxwell-Boltzmann distribution for our gas')
        ax1.set_xlabel('Energy of a molecule')
        ax1.set_ylabel('Number of instances')
        plt.show()
        
I = Simulation()
I.run(100, True)
#%%
'''Looking for a better ball initialisation method
with meshgrids.'''
x = np.arange(-10,10,3)

y = np.meshgrid(x,x)
print(y[0])
for i in range(5):
    a = (rand.choice(y[0]))
    a = a.tolist()
    b = rand.choice(a)
    a.pop(a.index(b))
    print(a, b)
print(y[0])
'''meshgrid hard to implement, better try lists'''
#%%
'''2D list to assign random coordinates.
Create a list of lists of a given range. 
each one of them also has an index, which is the last number.
Choose the list at random (pos_x), check if it has elements,
remove the index, choose an element at random (pos_y), remove that element
from the list, re-attach the index.
Repeat until all balls assigned coordinates.'''
size = range(-14,15,3)
a = []
b = []
for i in size:
    a.append(i)
for j in size:
    c = a.copy()
    c.append(j/3)
    b.append(c)
print(b)
for k in range(10):
    x = rand.choice(b) 
    print(x)
    while len(x) == 1:
        x = rand.choice(b)
    pos_x = b.index(x)
    x.pop()
    pos_y = rand.choice(x)
    x.remove(pos_y)
    x.append(pos_x)
    print(3*pos_x, pos_y)
    print(b)
'''A factor of 3*radius of a ball has to be introduced to
avoid overlap, since the original lists in b are indexed 
ordinarily (step 1). The index of each c list has to be divided
by 3*radius as well to bring them down to ordinary indexing.'''

#%%
'''Investigating storing large amount of data in a txt file.'''
import pickle

class Simulation:
    ball_list = []
    volume = 20
    pressure_hits = []
    T = 0
    ball1_id = 0  #collision ball ids
    ball2_id = 0
    Container_collision = True
    Container_collision = True
    tc = []  #container collisions list
    tb = []  #ball collisions list
    cb = []  #collision ball ids
    def __init__(self,m=1,R=1,r=np.array([0,0]),v=np.array([0,0]), container=False, NoB=1):
        self._NoB = 40  #maximum is (f^2)/(3^2), which is just a smaller square than half of square of radius / (3)^2 
        self._container = Ball(container=True)
        f = math.floor(np.sqrt(0.5*(Simulation.volume**2)))-1
        size = range(-f,f,3*R)  #sqrt of half of the square of (radius minus one) (size of possible grid)
        vx = []
        for i in range(-2*f,2*f):       #velocities are random, but container size dependent
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
            Simulation.ball_list.append(Ball(r= np.array([(3*R*pos_x)-f,pos_y]), v=np.array([2*rand.choice(vx),2*rand.choice(vy)])))

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
        print('cb is', Simulation.cb)
        print('tc is', Simulation.tc)
        print('tb is', Simulation.tb)
        t_cont_min = min(Simulation.tc)
        t_ball_min = min(Simulation.tb)
        print('t_cont_min is', t_cont_min)
        print('t_ball_min is', t_ball_min)
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
        else:
           Simulation.ball_list[Simulation.ball1_id].collide(self._container)
           print('ball', Simulation.ball1_id, 'collides with container')    
           print(Simulation.tc)
           print(Simulation.tb)
           
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
            Simulation.pressure_hits.append(0)
            Simulation.ball1_id = Simulation.tb.index(t_ball)   #get coliding ball numbers
            Simulation.ball2_id = Simulation.cb[Simulation.ball1_id]
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_ball) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_ball #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_ball
            Simulation.ball_list[Simulation.ball1_id].collide(Simulation.ball_list[Simulation.ball2_id])
        else:
            Simulation.Container_collision = True
            Simulation.ball1_id = Simulation.tc.index(t_cont)
            Simulation.pressure_hits.append((2*Simulation.ball_list[Simulation.ball1_id]._mass*np.sqrt(abs(np.dot(Simulation.ball_list[Simulation.ball1_id]._vel,Simulation.ball_list[Simulation.ball1_id]._vel))))/(2*np.pi*Simulation.volume))
            for i in range(self._NoB):
                Simulation.ball_list[i].move(t_cont) #move all the balls
                Simulation.tb[i] = Simulation.tb[i] - t_cont #increment all the times
                Simulation.tc[i] = Simulation.tc[i] - t_cont
            Simulation.ball_list[Simulation.ball1_id].collide(self._container)
            
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
        tot_mx_list = []
        tot_my_list = []
        vel_squares = []
        self.next_collision_initial()
        pl.pause(0.001)
        for frame in range(num_frames-1):
            energy_list = []
            mx = []
            my = []
            for i in range(self._NoB):
                centr_list.append(np.sqrt(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[i].pos())))
                energy_list.append((Simulation.ball_list[i]._mass*np.dot(Simulation.ball_list[i]._vel, Simulation.ball_list[i]._vel))/2)
                mx.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[0]))
                my.append(abs(Simulation.ball_list[i]._mass*Simulation.ball_list[i]._vel[1]))
                vel_squares.append(abs(np.dot(Simulation.ball_list[i]._vel,Simulation.ball_list[i]._vel)))
                for k in range(self._NoB):
                    if k != i:
                        dist_list.append(np.sqrt(abs(np.dot(Simulation.ball_list[i].pos(),Simulation.ball_list[k].pos()))))
            tot_energy_list.append(sum(energy_list))
            tot_mx_list.append(sum(mx))
            tot_my_list.append(sum(my))
            tot_m_list = [a-b for a,b in zip(tot_mx_list,tot_mx_list)]
            self.next_collision_main()

            if animate:
                    pl.pause(0.001)
            
        kB = 1
        T = [(x*2)/(3*kB*self._NoB) for x in tot_energy_list]
        average_pressure = (sum(Simulation.pressure_hits))/num_frames
        print(average_pressure)

        bolts = []
        for i in range(num_frames-1):
            for j in range(i*self._NoB, (i+1)*self._NoB):
                bolts.append((np.sqrt(vel_squares[j]))*np.exp(-(0.5*(vel_squares[j]))/(kB*T[i])))
        print(len(bolts))
        with open('listfile.data', 'wb') as filehandle:
            # store the data as binary data stream
            pickle.dump(bolts, filehandle)
I = Simulation()
I.run(1000, True)
#%%
'''Now loading the data and plotting a graph.'''
with open('listfile.data', 'rb') as filehandle:
    # read the data as binary data stream
    bolts = pickle.load(filehandle)

del bolts[:900*40]  #ignore first 200 frames for all 40 balls

boltzmann = plt.figure()
ax1 = boltzmann.add_subplot(1, 1, 1)
plt.hist(bolts, bins = 20)       
plt.title('Maxwell-Boltzmann distribution for our gas')
ax1.set_xlabel('Energy of a molecule')
ax1.set_ylabel('Number of instances')
plt.show()

