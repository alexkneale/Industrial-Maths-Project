#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 11:55:07 2023

@author: Alexander
"""

'''
This script finds a numerical approximation for the solution to
 the following 1D heat equation with Neumann conditions:
   u_t = u_xx  for  x \in (0,L),  t \in (0,T),
   u(x,t=0) = f(x),
   u'(x=0,t) = 0, u'(x=L,t) = 0,
 using a backward Euler scheme.
''' 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
import scipy.linalg # library needed for matrices
import scipy.integrate as integrate





# perform initial set up for number of points in time and space
L=1.0; T=4.0
Nx_spaces = 70; Nt_gaps = 200; 
Nx_points = Nx_spaces +  1 ; Nt_points = Nt_gaps + 1
x = np.linspace(0, L, Nx_points)   # mesh points in space
dx = x[1] - x[0]

t = np.linspace(0, T, Nt_points) # mesh points in time
dt = t[1] - t[0]
#diffusivity of disease
D = 0.1
C = D*dt/dx**2

def U_exact(x,t,x_centre):
    M = np.size(x)
    u_ex = np.zeros(M)  
    #calculated manually
    
    
    limit1 = x_centre - dx/2
    limit2 = x_centre + dx/2
    c_0 = limit2-limit1
    u_ex = u_ex + c_0

    for n in range(1,10000):
        

        c_n = (2/(n*np.pi))*(np.sin(n*np.pi*limit2)-np.sin(n*np.pi*limit1))

        u_ex = u_ex + c_n*np.cos(n*np.pi*x/L)*np.exp(-1*(n*np.pi/L)**2*D*t)
    
    return u_ex


# set up structures to hold U and U_ex2 and interim arrays
u   = np.zeros(Nx_points)
u_old = np.zeros(Nx_points)
U = np.zeros((Nx_points,Nt_points))
U_ex = np.zeros((Nx_points,Nt_points))

'''
# to hold U_exact
U_ex2 = np.zeros((Nx_points,Nt_points))
'''

# Data structures for the linear system
A = np.zeros((Nx_points, Nx_points))
b = np.zeros(Nx_points)

# set up the matrix A
for i in range(1, Nx_points-1): # rows from 1 to Nx-2
    A[i,i-1] = -C
    A[i,i+1] = -C
    A[i,i] = 1 + 2*C
    
A[0,0] = 1+2*C  ; A[0,1]= -2*C 
A[Nx_points-1,Nx_points-1] = 1+2*C 
A[Nx_points-1,Nx_points-2] = -2*C


# function for setting the initial condition in space  I(x)
def I2(x):
    n = x.size
    I2_arr = np.zeros(n)
    I2_arr[round(n/2)] = 1
    return I2_arr



# Set initial condition u(x,0) = I(x)
u_old = I2(x) # no dirichlet boundary conditions in this example
n = x.size
x_centre = x[round(n/2)]

# initialise matrices U and U_ex2 for first time step
U[:,0] = u_old[:]


U_ex[:,0]=U_exact(x,0,x_centre)

#perform time=stepping
for n in range(1, Nt_points): # timestep for 1 to t = T-1 so last step finishes on t=T
    # Compute b and solve linear system
    b[:] = u_old[:]
    u[:] = np.linalg.solve(A,b)
    # Update u_1 before next step
    u_old = u
    U[:,n] = u
    
    U_ex[:,n]=U_exact(x,t[n],x_centre)
    

# set up animation plots
fig, ax = plt.subplots()
ax.axis([0,L,0,L])
def animate(i):
    l.set_data(x, U[:,i])

    m.set_data(x,U_ex[:,i])
    
    
ax.axis([0,L,0,L])
l, = ax.plot([],[],':r')

m, = ax.plot([],[],'-.b')

ani2 = matplotlib.animation.FuncAnimation(fig, animate, frames=Nt_points)

from IPython.display import HTML
HTML(ani2.to_jshtml())