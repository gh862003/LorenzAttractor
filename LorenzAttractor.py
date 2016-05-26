# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:53:37 2016

@author: christiana
"""
import numpy as np

def lorenz_derivatives(x,y,z,lorenzParam):
    '''
    takes the x y and z values and params and solves each
    of the lorenz derivatives for use in the Runge Kutta update scheme
    '''
    sigma = lorenzParam['sigma']
    rho = lorenzParam['rho']
    beta = lorenzParam['beta']
    
    dx = -sigma*(x-y)
    dy = rho*x -y -x*z
    dz = x*y - beta*z
    
    
    return dx,dy,dz
    
def Euler_step(x,y,z,lorenzParam):
    dt = lorenzParam['dt']/4.
    dx,dy,dz = lorenz_derivatives(x,y,z,lorenzParam)
    x_new = x + dt*dx
    y_new = y +dt*dy
    z_new = z +dt*dz
    return x_new,y_new,z_new
    
def CN(lorenzParam):
    '''
    runs with a qurter of the timestep of rungeKutta
    '''
        
    #Access needed params
    nt = 10*lorenzParam['nt']
    dt = lorenzParam['dt']/10.
    x_o = lorenzParam['x_o']
    y_o = lorenzParam['y_o']
    z_o = lorenzParam['z_o']
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    x[0]=x_o
    y[0]=y_o
    z[0]=z_o
    
    for n in range(0,nt):
        
        x_eul,y_eul,z_eul= Euler_step(x[n],y[n],z[n],lorenzParam)
        
        dx_n,dy_n,dz_n = lorenz_derivatives(x[n],y[n],z[n],lorenzParam)
        
        dx_n1,dy_n1,dz_n1 = lorenz_derivatives(x_eul,y_eul,z_eul,lorenzParam)
        
        x[n+1] = x[n] +(dx_n +dx_n1)*dt/2.
        y[n+1] = y[n] +(dy_n +dy_n1)*dt/2.
        z[n+1] = z[n] +(dz_n +dz_n1)*dt/2.
    
    return x,y,z
    
    

    
def RungeKutta(lorenzParam,factor,x_o,y_o,z_o,dt):                 
    '''
    Takes the input paramaters and preforms a 4 stage runge kutta method
    for the lorenz eqns returns the values of x,y and z at all timesteps
    including intiial conditions
    '''
    
    #Access needed params
    nt = int(lorenzParam['nt']*factor)

    dt = lorenzParam['dt']/factor

    
    #nt+1 as we want space to store information from nt timesteps PLUS the initial conditions
    x = np.zeros(nt+1)
    y = np.zeros(nt+1)
    z = np.zeros(nt+1)
    x[0]=x_o
    y[0]=y_o
    z[0]=z_o
    r = np.zeros(nt+1)
    v = np.zeros(nt+1)
    w = np.zeros(nt+1)
    r[0] =0
    v[0] =0

    
    for n in range(0,nt):

    
        #At each timestep we must calculate k(1-4) for x(k),y(l) and z(m) 
        k1,l1,m1 = lorenz_derivatives(x[n],y[n],z[n],lorenzParam)
        
        k2,l2,m2 = lorenz_derivatives(x[n]+(dt/2.)*k1,y[n]+(dt/2.)*l1,z[n]+(dt/2.)*m1,\
                  lorenzParam)
                
        k3,l3,m3 = lorenz_derivatives(x[n]+(dt/2.)*k2, y[n]+(dt/2.)*l2,z[n]+(dt/2.)*m2,\
                  lorenzParam)
                
        k4,l4,m4 = lorenz_derivatives(x[n]+dt*k3,y[n]+dt*l3,z[n]+dt*m3,\
                 lorenzParam)
        
        #Update Eqns
        x[n+1] = x[n] + (dt/6.)*(k1+2*k2+2*k3+k4)
        y[n+1] = y[n] + (dt/6.)*(l1+2*l2+2*l3+l4)
        z[n+1] = z[n] + (dt/6.)*(m1+2*m3+2*m3+m4)
        #Calculate the speed to find the distance 
        v= (k1+2*k2+2*k3+k4)/6.
        u= (l1+2*l2+2*l3+l4)/6.
        w= (m1+2*m3+2*m3+m4)/6.
        V = np.sqrt(v**2 +u**2 +w**2)
        r[n+1] = r[n] + V*dt
        
    return x,y,z,r