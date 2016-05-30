# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:53:01 
@author: christiana
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np




def errorPlot(diff, lorenzParam,name,dttype):
    dt = lorenzParam['dt']
    nt = lorenzParam['nt']
    duration = dt*nt
    t = np.linspace(0,duration,nt)
    plt.plot(t,diff, 'k')
    plt.xlabel('duration, t')
    plt.ylabel('distance travelled error:RK_'+str(name)+' and RK'+str(dttype))
    plt.show()
    

def lorentzPlotting(x,y,z,lorenzParam):
    sigma = lorenzParam['sigma']
    rho = lorenzParam['rho']
    beta = lorenzParam['beta']
    nt = lorenzParam['nt']
    dt = lorenzParam['dt']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x,y,z)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Lorenz attractor. Duration:%.2f, timestep:%.4f,rho:%d,beta:%.2f,sigma :%d ' %(nt*dt,dt,rho,beta,sigma))
    fig.show()
    
    
    

