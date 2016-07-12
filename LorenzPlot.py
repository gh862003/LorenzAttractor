# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:53:01 
Contains various plotting functions for solutions and errors of lorenz attractor
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np




def errorPlot(diff, lorenzParam,name,dttype,nt):
    '''
    Plots the error for 1 (constant time-step) case
    '''
    dt = lorenzParam['dt']
  
    duration = dt*nt
    diff = [np.log(x) for x in diff]
    t = np.linspace(0,duration,len(diff))

    plt.plot(t,diff, 'k')
    plt.xlabel('duration, t')
    plt.ylabel('distance travelled error:'+str(name)+' and RK4 High Res')
    plt.savefig('LongRunPlot_' + str(name) + '.pdf')
    plt.show()
    

def errorCompar(d_const,d_var,dt,nt,name,dt_var,var_type):
    '''
    plots two errors one constant one variable
    d_const is the errors of the constant case, d_var - erors of variable case
    dt_var is the variable case time steps
    '''
    
    duration = dt*nt
    #const case time axis
    t = np.linspace(0,duration,len(d_const))
    t_var = np.cumsum(dt_var)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.plot(t,d_const, 'k')
    ax1.plot(t_var,d_var,'firebrick')
    ax1.set_xlabel('duration, t')
    ax1.set_ylabel('distance error wtr RK: '+str(name)+' and '+str(name)+ 'variable')
    
    
    #second x axis that shows the size of the timestep
    second_x = np.array(dt_var)
    x_tick_labels = ["%.3f" % z for z in dt_var]
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(t_var)
    ax2.set_xticklabels(x_tick_labels)
    plt.savefig('ComparisionPlot_'+str(name)+'variable_'+str(var_type)+'.pdf')
    plt.show()

    
    

def lorentzPlotting(x,y,z,lorenzParam,nt):
    '''
    Plots the solution of the lorenz equations
    '''
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
    fig.savefig('Attractor.pdf')
    fig.show()
    
def ICplot(x,y,z,lorenzParam,xs,ys,zs):
    '''
    plots the chosen 100 initial contions x,y,z out of the larger section xs,ys,zs
    xyz are points on a xs ys zs solution
    '''
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xs,ys,zs,c = 'royalblue',alpha = 0.6)
    ax.scatter(x,y,z,'o', alpha =1.0)
    
    plt.gca().patch.set_facecolor('white')
    
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')  
    fig.savefig('IC_distribution_plot.pdf')
    fig.show()

