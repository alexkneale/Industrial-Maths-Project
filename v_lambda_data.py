#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 17:17:18 2023

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


def v_func(x,x_screen_1,x_screen_2,v1,v2):
    if x < x_screen_1:
        return v1
    elif x >= x_screen_1:
        if x < x_screen_2:
            return v2
        else:
            return v1





def loop_U(x_screen_1,x_screen_2,v1,v2,lambd):

    v = np.zeros(Nx_points)
    i = 0
    for x_i in x:
        v[i] = v_func(x_i,x_screen_1,x_screen_2,v1,v2)
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
    
    
    # function for setting the initial condition in space  I(x)
    def I2(x):
        n = x.size
        I2_arr = np.zeros(n)
        
        I2_arr[round(n/2)] = 1
        return I2_arr
    
    
    # Set initial condition u(x,0) = I(x)
    u_old = I2(x) # no dirichlet boundary conditions in this example
    
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
        print(np.sum(u_old))
        
        
    return U

def loop_pos(x_screen_1,x_screen_2,v1,v2,lambd,x_distance,t_arr):
    
    U = loop_U(x_screen_1,x_screen_2,v1,v2,lambd)
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

def loop_area(x_screen_1,x_screen_2,v1,v2,lambd,t_arr):

    U = loop_U(x_screen_1,x_screen_2,v1,v2,lambd)
    
    
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
    



x_screen_1 = 0.6
x_screen_2 = 0.8
v1 = 0.3
v2 = 0.1
lambd = 0.05
x_distance = 0.2
t_arr = np.array([0.2,0.5,2.0])

U_t = loop_pos(x_screen_1,x_screen_2,v1,v2,lambd,x_distance,t_arr)
print(U_t)

A_t = loop_area(x_screen_1,x_screen_2,v1,v2,lambd,t_arr)
print(A_t)



