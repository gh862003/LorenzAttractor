# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:53:01 
@author: christiana
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np




def errorPlot(diff, lorenzParam,name,dttype,nt):
    dt = lorenzParam['dt']
  
    duration = dt*nt
    
    t = np.linspace(0,duration,len(diff))
    plt.plot(t,diff, 'k')
    plt.xlabel('duration, t')
    plt.ylabel('distance travelled error:'+str(name)+' and RK4 High Res')
    plt.savefig('LongRunPlot_' + str(name) + '.pdf')
    plt.show()
    

def errorCompar(d_const,d_var,dt,nt,name,t_var):
    duration = dt*nt
    
    t = np.linspace(0,duration,len(d_const))
    plt.plot(t,d_const, 'k')
    plt.plot(t_var,d_var,'firebrick')
    plt.xlabel('duration, t')
    plt.ylabel('distance error wtr RK: '+str(name)+' and '+str(name)+ 'variable')
    plt.show()

    
    

def lorentzPlotting(x,y,z,lorenzParam,nt):
    sigma = lorenzParam['sigma']
    rho = lorenzParam['rho']
    beta = lorenzParam['beta']
    dt = lorenzParam['dt']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x,y,z)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Lorenz attractor. Duration:%.2f, timestep:%.4f,rho:%d,beta:%.2f,sigma :%d ' %(nt*dt,dt,rho,beta,sigma))
    fig.show()
    
    
    

