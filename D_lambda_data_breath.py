#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 20:51:08 2023

@author: Alexander
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import scipy.linalg # library needed for matrices
from SQR import SQR



# perform initial set up for number of points in time and space
L=1.0; T=3.0
Nx_spaces = 50; Nt_gaps = 150; 
Nx_points = Nx_spaces +  1 ; Nt_points = Nt_gaps + 1
x = np.linspace(0, L, Nx_points)   # mesh points in space
dx = x[1] - x[0]

t = np.linspace(0, T, Nt_points) # mesh points in time
dt = t[1] - t[0]
#diffusivity of disease



def D_func(x,x_screen_1,x_screen_2,D1,D2):
    if x < x_screen_1:
        return D1
    elif x >= x_screen_1:
        if x < x_screen_2:
            return D2
        else:
            return D1

def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

def loop_U(x_screen_1,x_screen_2,D1,D2,lambd):
    
    #index of position of person breathing in room
    n = x.size
    x_person =  (n/4)*dx
    
    
    
    
    D = np.zeros(Nx_points)
    i = 0
    for x_i in x:
        D[i] = D_func(x_i,x_screen_1,x_screen_2,D1,D2)
        i += 1



    # set up structures to hold U and U_ex2 and interim arrays
    u   = np.zeros(Nx_points)
    u_old = np.zeros(Nx_points)
    U = np.zeros((Nx_points,Nt_points))
    
    
    # Data structures for the linear system
    A = np.zeros((Nx_points, Nx_points))
    b = np.zeros(Nx_points)
    
    # set up the matrix A
    for i in range(1, Nx_points-1): # rows from 1 to Nx-2
        A[i,i-1] = -D[i]*dt/(dx**2) + ((D[i+1]-D[i-1])*dt)/(4*dx**2)
        A[i,i+1] = -D[i]*dt/(dx**2) - ((D[i+1]-D[i-1])*dt)/(4*dx**2)
        A[i,i] = 1 + 2*D[i]*dt/(dx**2) + lambd*dt
        
    A[0,0] = 1+2*D[0]*dt/(dx**2)+lambd*dt
    A[0,1]= -2*D[0]*dt/(dx**2)
    A[Nx_points-1,Nx_points-1] = 1+2*D[Nx_points-1]*dt/(dx**2)+lambd*dt
    A[Nx_points-1,Nx_points-2] = -2*D[Nx_points-1]*dt/(dx**2)
    


    
    
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = 0.05*gaussian(x,x_person,sig) # no dirichlet boundary conditions in this example
    
    # initialise matrices U and U_ex2 for first time step
    U[:,0] = u_old[:]
    '''
    n_breath = [4,8,12]
    n = x.size
    pos_breath = round(n/4)
    n_inhale = [6,10,14]
    '''
    #perform time=stepping
    for n in range(1, Nt_points): # timestep for 1 to t = T-1 so last step finishes on t=T
        # Compute b and solve linear system
        b[:] = u_old[:]
        u[:] = SQR(A,b)
        '''
        if n in n_breath:
            u[pos_breath] = u[pos_breath] + 0.5
        elif n in n_inhale:
            u[pos_breath] = u[pos_breath] - 0.1
            
        '''
        # Update u_1 before next step
        u_old = u
        U[:,n] = u
        
    return U

def loop_U_breath(x_screen_1,x_screen_2,D1,D2,lambd):
    
    #index of position of person breathing in room
    n = x.size
    x_person =  (n/4)*dx
    
    index_person = round((n/4))

    
    
    
    D = np.zeros(Nx_points)
    i = 0
    for x_i in x:
        D[i] = D_func(x_i,x_screen_1,x_screen_2,D1,D2)
        i += 1



    # set up structures to hold U and U_ex2 and interim arrays
    u   = np.zeros(Nx_points)
    u_old = np.zeros(Nx_points)
    U = np.zeros((Nx_points,Nt_points))
    
    
    # Data structures for the linear system
    A = np.zeros((Nx_points, Nx_points))
    b = np.zeros(Nx_points)
    
    # set up the matrix A
    for i in range(1, Nx_points-1): # rows from 1 to Nx-2
        A[i,i-1] = -D[i]*dt/(dx**2) + ((D[i+1]-D[i-1])*dt)/(4*dx**2)
        A[i,i+1] = -D[i]*dt/(dx**2) - ((D[i+1]-D[i-1])*dt)/(4*dx**2)
        A[i,i] = 1 + 2*D[i]*dt/(dx**2) + lambd*dt
        
    A[0,0] = 1+2*D[0]*dt/(dx**2)+lambd*dt
    A[0,1]= -2*D[0]*dt/(dx**2)
    A[Nx_points-1,Nx_points-1] = 1+2*D[Nx_points-1]*dt/(dx**2)+lambd*dt
    A[Nx_points-1,Nx_points-2] = -2*D[Nx_points-1]*dt/(dx**2)
    


    
    
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = 0.05*gaussian(x,x_person,sig) # no dirichlet boundary conditions in this example
    
    # initialise matrices U and U_ex2 for first time step
    U[:,0] = u_old[:]
    
    n_breath = [25,50,75,100,125]
    n_inhale = [13,37,62,87,112]
    
    #perform time=stepping
    for n in range(1, Nt_points): # timestep for 1 to t = T-1 so last step finishes on t=T
        # Compute b and solve linear system
        b[:] = u_old[:]
        u[:] = SQR(A,b)
        
        if n in n_breath:
            u[:] = u[:] + 0.01*gaussian(x,x_person,sig)
        elif n in n_inhale:
            u[:] = u[:] - 0.01*gaussian(x,x_person,sig)
            
        # Update u_1 before next step
        u_old = u
        U[:,n] = u
        
    return U

def loop_pos(x_screen_1,x_screen_2,D1,D2,lambd,x_distance,t_arr):
    
    U = loop_U(x_screen_1,x_screen_2,D1,D2,lambd)
    #x position of the screen's centre
    x_screen_centre = (x_screen_2 + x_screen_1)/2
    #x position distance x_distance from x_screen_centre
    x_target = x_screen_centre + x_distance
    
    
    
    # find index of x pos which is closest to x_target
    x_abs = abs(x - x_target)
    x_index = np.where(x_abs == x_abs.min())[0]


    
    
    
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    U_t = np.zeros(n_t)
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        U_t[i] = U[x_index,t_index]
        
        i += 1
        
    return U_t
    
def loop_area(x_screen_1,x_screen_2,D1,D2,lambd,t_arr):

    U = loop_U(x_screen_1,x_screen_2,D1,D2,lambd)
    
    total_virus = np.sum(U[:,0])
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    area_arr = np.zeros(n_t)
    
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        area_arr[i] = area_arr[i] + np.sum(U[x>= x_screen_2,t_index])/total_virus
        i += 1
    return area_arr


def loop_area_breath(x_screen_1,x_screen_2,D1,D2,lambd,t_arr):

    U = loop_U_breath(x_screen_1,x_screen_2,D1,D2,lambd)
    
    total_virus = np.sum(U[:,0])
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    area_arr = np.zeros(n_t)
    
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        area_arr[i] = area_arr[i] + np.sum(U[x>= x_screen_2,t_index])/total_virus
        
        i += 1
    return area_arr

def U_loop_times(x_screen_1,x_screen_2,D1,D2,lambd,t_arr):
    
    U = loop_U(x_screen_1,x_screen_2,D1,D2,lambd)
    n_t = t_arr.shape[0]
    
    U_t = np.zeros([Nx_points,n_t])
    #find indices of times that are closest to those in t_arr, find U and append to U_t
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr

        t_abs = abs(t - t_i)

        t_index = np.where(t_abs == t_abs.min())[0]

        U_t[:,i] = U[:,t_index][:,0]
        i += 1
    return U_t

#1st test_data
'''
L=1.0; T=3.0
Nx_spaces = 50; Nt_gaps = 150; 
Nx_points = Nx_spaces +  1 ; Nt_points = Nt_gaps + 1
x = np.linspace(0, L, Nx_points)   # mesh points in space
dx = x[1] - x[0]

t = np.linspace(0, T, Nt_points) # mesh points in time
dt = t[1] - t[0]

sig = 0.04

D1 = 0.5
D2_arr = np.array([0.1,0.4])
x_screen_1 = 0.7
x_screen_2 = 0.8
lambd = 1.0
x_distance = 0.1
t_arr = np.array([0.02,0.6,1.0,2.0,2.9])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,L,0,L])

colours = ['r','b'] # make comparison easy


for i in range(len(colours)):
    print(D2_arr[i])
    U = U_loop_times(x_screen_1,x_screen_2,D1,D2_arr[i],lambd,t_arr)

    for t_i in range(len(t_arr)):
        label = "D2 = "+str(D2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.plot(x,U[:,t_i],linestyle = ':',color = colours[i], label=label)

plt.show()

'''
'''
#curves at different times for two diff values of D2 (graph 1)
sig = 0.04
   
D1 = 0.5 
D2_arr = np.array([0.1,0.4])
x_screen_1 = 0.7
x_screen_2 = 0.8
lambd = 1.0
x_distance = 0.1
t_arr = np.array([0.02,0.6,1.0,2.0,2.9])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,L,0,L])

colours = ['r','b'] # make comparison easy


for i in range(len(colours)):
    print(D2_arr[i])
    U = U_loop_times(x_screen_1,x_screen_2,D1,D2_arr[i],lambd,t_arr)

    for t_i in range(len(t_arr)):
        label = "D2 = "+str(D2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.plot(x,U[:,t_i],linestyle = ':',color = colours[i], label=label)

plt.show()
'''
'''
#without breathing data for percentage over time varying D2

sig = 0.04

sig = 0.04
   
D1 = 0.3
D2_arr = np.array([0.08,0.1,0.12,0.15,0.2])
x_screen_1 = 0.6
x_screen_2 = 0.8
lambd = 0.1

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    U_t = loop_area(x_screen_1,x_screen_2,D1,D2_arr[i],lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "v2 = "+str(D2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()
'''
'''
#with breathing data for percentage over time varying D2


sig = 0.04
   
D1 = 0.3
D2_arr = np.array([0.08,0.1,0.12,0.15,0.2])
x_screen_1 = 0.6
x_screen_2 = 0.8
lambd = 0.1

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    U_t = loop_area_breath(x_screen_1,x_screen_2,D1,D2_arr[i],lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "v2 = "+str(D2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()
'''

'''
#without breathing data for percentage over time varying wall size

sig = 0.04

sig = 0.04
   
D1 = 0.3
D2 = 0.2
x_screen_1 = 0.6

x_dist_arr = np.array([0.05,0.12,0.15,0.2,0.25])
lambd = 0.1

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    x_screen_2 = x_screen_1 + x_dist_arr[i]
    U_t = loop_area(x_screen_1,x_screen_2,D1,D2,lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "x dist = "+str(x_dist_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()
'''

#with breathing data for percentage over time varying wall size

sig = 0.04

sig = 0.04
   
D1 = 0.3
D2 = 0.2
x_screen_1 = 0.6

x_dist_arr = np.array([0.05,0.12,0.15,0.2,0.25])
lambd = 0.1

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    x_screen_2 = x_screen_1 + x_dist_arr[i]
    U_t = loop_area_breath(x_screen_1,x_screen_2,D1,D2,lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "x dist = "+str(x_dist_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()

'''
U = loop_U(x_screen_1,x_screen_2,D1,D2,lambd)

def animate(i):
    l.set_data(x, U[:,i])

    
ax.axis([0,L,0,L])
l, = ax.plot([],[],':r')

ani2 = matplotlib.animation.FuncAnimation(fig, animate, frames=Nt_points)

from IPython.display import HTML
HTML(ani2.to_jshtml())
'''
