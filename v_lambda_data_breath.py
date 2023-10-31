#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 20:16:55 2023

@author: Alexander
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import scipy.linalg # library needed for matrices



# perform initial set up for number of points in time and space
L=1.0; T=3.0
Nx_spaces = 100; Nt_gaps = 150;
Nx_points = Nx_spaces +  1 ; Nt_points = Nt_gaps + 1
x = np.linspace(0, L, Nx_points)   # mesh points in space
dx = x[1] - x[0]

t = np.linspace(0, T, Nt_points) # mesh points in time
dt = t[1] - t[0]
#diffusivity of disease
D = 0.1


def v_func(x,x_screen,v1,v2):
    if x < x_screen:
        return v1
    elif x >= x_screen:
        
        return v2


def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)



def loop_U(x_screen,v1,v2,lambd):

    
    # position of person breathing in room
    n = x.size
    x_person =  (n/4)*dx
    

    v = np.zeros(Nx_points)
    i = 0
    for x_i in x:
        v[i] = v_func(x_i,x_screen,v1,v2)
        i += 1

    
    
    # set up structures to hold U and U_ex2 and interim arrays
    u   = np.zeros(Nx_points)
    u_old = np.zeros(Nx_points)
    U = np.zeros((Nx_points,Nt_points))
    '''
    # to hold U_exact
    U_ex2 = np.zeros((Nx_points,Nt_points))
    '''
    
    # Data structures for the linear system
    A = np.zeros((Nx_points, Nx_points))
    b = np.zeros(Nx_points)
    
    # set up the matrix A
    #v has +1 as it has Nx_points+1 points, to incorporate v[-1] case
    for i in range(1, Nx_points-1): # rows from 1 to Nx-2
        A[i,i-1] = -D*dt/(dx**2) - v[i]*dt/(2*dx)
        A[i,i+1] = -D*dt/(dx**2) + v[i]*dt/(2*dx)
        A[i,i] = 1 + 2*D*dt/(dx**2) + (dt/(2*dx))*(v[i]-v[i-1]) + lambd*dt

    A[0,0] =  (v[0]**2)*dt/D + 1 + 2*D*dt/(dx**2) + (dt/dx)*(v[1]+v[0])+ lambd*dt
    A[0,1]= -2*D*dt/(dx**2)
    A[Nx_points-1,Nx_points-1] = 1 + 2*D*dt/(dx**2) - (dt/dx)*(v[Nx_points-1] + v[Nx_points-2]) + (v[Nx_points-1]**2)*dt/D + lambd*dt
    A[Nx_points-1,Nx_points-2] = -2*D*dt/(dx**2)
    
    
    
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = 0.05*gaussian(x,x_person,sig)  # no dirichlet boundary conditions in this example
    
    # initialise matrices U and U_ex2 for first time step
    U[:,0] = u_old[:]
    '''
    U_ex2[:,0]=U_exact2(x,0)
    '''
    
    #perform time=stepping
    for n in range(1, Nt_points): # timestep for 1 to t = T-1 so last step finishes on t=T
        # Compute b and solve linear system
        b[:] = u_old[:]
        u[:] = np.linalg.solve(A,b)
        # Update u_1 before next step
        u_old = u
        U[:,n] = u
        
        
    return U

def loop_U_breath(x_screen,v1,v2,lambd):

    
    # position of person breathing in room
    n = x.size
    x_person =  (n/4)*dx
    index_person = round((n/4))

    v = np.zeros(Nx_points)
    i = 0
    for x_i in x:
        v[i] = v_func(x_i,x_screen,v1,v2)
        i += 1

    
    
    # set up structures to hold U and U_ex2 and interim arrays
    u   = np.zeros(Nx_points)
    u_old = np.zeros(Nx_points)
    U = np.zeros((Nx_points,Nt_points))
    '''
    # to hold U_exact
    U_ex2 = np.zeros((Nx_points,Nt_points))
    '''
    
    # Data structures for the linear system
    A = np.zeros((Nx_points, Nx_points))
    b = np.zeros(Nx_points)
    
    # set up the matrix A
    #v has +1 as it has Nx_points+1 points, to incorporate v[-1] case
    for i in range(1, Nx_points-1): # rows from 1 to Nx-2
        A[i,i-1] = -D*dt/(dx**2) - v[i]*dt/(2*dx)
        A[i,i+1] = -D*dt/(dx**2) + v[i]*dt/(2*dx)
        A[i,i] = 1 + 2*D*dt/(dx**2) + (dt/(2*dx))*(v[i]-v[i-1]) + lambd*dt

    A[0,0] =  (v[0]**2)*dt/D + 1 + 2*D*dt/(dx**2) + (dt/dx)*(v[1]+v[0])+ lambd*dt
    A[0,1]= -2*D*dt/(dx**2)
    A[Nx_points-1,Nx_points-1] = 1 + 2*D*dt/(dx**2) - (dt/dx)*(v[Nx_points-1] + v[Nx_points-2]) + (v[Nx_points-1]**2)*dt/D + lambd*dt
    A[Nx_points-1,Nx_points-2] = -2*D*dt/(dx**2)
    
    
    
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = 0.05*gaussian(x,x_person,sig)  # no dirichlet boundary conditions in this example
    
    # initialise matrices U and U_ex2 for first time step
    U[:,0] = u_old[:]
    
    n_breath = [25,50,75,100,125]
    n_inhale = [13,37,62,87,112]
    
    #perform time=stepping
    for n in range(1, Nt_points): # timestep for 1 to t = T-1 so last step finishes on t=T
        # Compute b and solve linear system
        b[:] = u_old[:]
        u[:] = np.linalg.solve(A,b)
        
        if n in n_breath:
            u[:] = u[:] + 0.01*gaussian(x,x_person,sig)
        elif n in n_inhale:
            u[:] = u[:] - 0.01*gaussian(x,x_person,sig)
            
        
        # Update u_1 before next step
        u_old = u
        U[:,n] = u
        
        
    return U


def loop_pos(x_screen,v1,v2,lambd,x_distance,t_arr):
    
    U = loop_U(x_screen,v1,v2,lambd)
    #x position of the screen's centre
    #x position distance x_distance from x_screen_centre
    x_target = x_screen + x_distance
    
    
    
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

def loop_area(x_screen,v1,v2,lambd,t_arr):

    U = loop_U(x_screen,v1,v2,lambd)
    
    total_virus = np.sum(U[:,0])
    
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    area_arr = np.zeros(n_t)
    
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        area_arr[i] = area_arr[i] + np.sum(U[x>= x_screen,t_index])/total_virus
        
        i += 1
    return area_arr

def loop_area_breath(x_screen,v1,v2,lambd,t_arr):

    U = loop_U_breath(x_screen,v1,v2,lambd)
    
    total_virus = np.sum(U[:,0])
    
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    area_arr = np.zeros(n_t)
    
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        area_arr[i] = area_arr[i] + np.sum(U[x>= x_screen,t_index])/total_virus
        
        i += 1
    return area_arr

def U_loop_times(x_screen,v1,v2,lambd,t_arr):
    
    U = loop_U(x_screen,v1,v2,lambd)
    
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

'''
photo 2 data


# perform initial set up for number of points in time and space
L=1.0; T=3.0
Nx_spaces = 100; Nt_gaps = 150;
Nx_points = Nx_spaces +  1 ; Nt_points = Nt_gaps + 1
x = np.linspace(0, L, Nx_points)   # mesh points in space
dx = x[1] - x[0]

t = np.linspace(0, T, Nt_points) # mesh points in time
dt = t[1] - t[0]
#diffusivity of disease
D = 0.1

sig = 0.04

x_screen = 0.7

v1 = 0.4
v2_arr = np.array([0.3,0.1])
lambd = 0.1
x_distance = 0.2

t_arr = np.array([0.02,0.6,1.0,2.0,2.9])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,L,0,L])

colours = ['r','b'] # make comparison easy


for i in range(len(colours)):
    U_t = U_loop_times(x_screen,v1,v2_arr[i],lambd,t_arr)

    for t_i in range(len(t_arr)):
        label = "v2 = "+str(v2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.plot(x,U_t[:,t_i],linestyle = ':',color = colours[i], label=label)

plt.show()



'''

#without breathing data for percentage over time

sig = 0.04

x_screen = 0.7

v1 = 0.4
v2_arr = np.array([0.02,0.1,0.2,0.3,0.4])
lambd = 0.1
x_distance = 0.2

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    U_t = loop_area(x_screen,v1,v2_arr[i],lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "v2 = "+str(v2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()
'''
#breath data area percentage
sig = 0.04

x_screen = 0.7

v1 = 0.4
v2_arr = np.array([0.02,0.1,0.2,0.3,0.4])
lambd = 0.1
x_distance = 0.2

t_arr = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])


# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,t_arr[-1],0,L])

colours = ['r','g','b','purple','yellow'] # make comparison easy

for i in range(len(colours)):
    U_t = loop_area_breath(x_screen,v1,v2_arr[i],lambd,t_arr)
    for t_i in range(len(t_arr)):
        print(t_arr)
        label = "v2 = "+str(v2_arr[i])+" t=" + "%0.3f" % (t_arr[t_i])
        ax.scatter(t_arr,U_t,linestyle = ':',color = colours[i], label=label)

plt.show()
'''
'''
sig = 0.04

x_screen = 0.7
lambd = 0.1
v1 = 0.4
v2 = 0.03
U = loop_U_breath(x_screen,v1,v2,lambd)
fig, ax = plt.subplots()
ax.axis([0,L,0,L])
def animate(i):
    l.set_data(x, U[:,i])

    
ax.axis([0,L,0,L])
l, = ax.plot([],[],':r')

ani2 = matplotlib.animation.FuncAnimation(fig, animate, frames=Nt_points)

from IPython.display import HTML
HTML(ani2.to_jshtml())
'''
