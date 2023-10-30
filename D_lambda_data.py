#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 16:03:51 2023

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

D1 = 0.4



def D_func(x,x_screen_1,x_screen_2,D1,D2):
    if x < x_screen_1:
        return D1
    elif x >= x_screen_1:
        if x < x_screen_2:
            return D2
        else:
            return D1

def loop_U(x_screen_1,x_screen_2,D1,D2,lambd):
    

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
    


    # function for setting the initial condition in space  I(x)
    def I2(x):
        n = x.size
        I2_arr = np.zeros(n)
        I2_arr[round(n/4)] = 1
        return I2_arr
    
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = I2(x) # no dirichlet boundary conditions in this example
    
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
    
    
    # number of times in t_arr
    n_t = t_arr.shape[0]
    
    area_arr = np.zeros(n_t)
    
    
    i = 0
    for t_i in t_arr:
        
        #find index of discretised time value closest to t_i in t_arr
        
        t_abs = abs(t - t_i)
        t_index = np.where(t_abs == t_abs.min())[0]
        
        area_arr[i] = area_arr[i] + np.sum(U[x>= x_screen_2,t_index])
        
        i += 1
    return area_arr
    
    
    

D2 = 0.3
x_screen_1 = 0.7
x_screen_2 = 0.75
lambd = 0.1
x_distance = 0.1
t_arr = np.array([0.2,0.5,0.7,0.9,1.2,1.9,2.9])

print(loop_area(x_screen_1,x_screen_2,D1,D2,lambd,t_arr))

# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,L,0,L])
'''
def animate(i):
    l.set_data(x, U[:,i])
    
    
ax.axis([0,L,0,L])
l, = ax.plot([],[],':r')

ani2 = matplotlib.animation.FuncAnimation(fig, animate, frames=Nt_points)

from IPython.display import HTML
HTML(ani2.to_jshtml())
'''