# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 13:44:19 2016
Here, the order of acc test
@author: christiana
"""

import matplotlib.pyplot as plt
import numpy as np
import lorenzParams
import LorenzAttractor


def orderOfAccuracy():
    
   #Set up dictionary of input parameters from lorenzParams file
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 100.
   
   d_err_ICs =[]
   d_err_dts =[]
   
   d_RKA_av =[]
   d_RKAs =[]
   l2s =[]
   l2Err =[]
   

   

   dt_range = [0.00075,0.0015,0.003]
   d_err_dts = [[],[],[]]

   
   #loop over the number of initial conditions to use

       
   for i in range(0,10):
   
       #create random intial conditions between -15 and 15
       x_o = -15 + 30 * np.random.random()
       y_o = -15 + 30 * np.random.random()
       z_o = -15 + 30 * np.random.random()
       '''
       x_o =ICs[i]
       y_o =ICs[(i+3)%len(ICs)]
       z_o =ICs[(2*i)%len(ICs)]
       '''
       
       
       for index,dt in enumerate(dt_range):
           
           lorenzParam['nt_RK'] = int(round(7./dt))
           lorenzParam['nt_Eu'] = int(round(7./dt))
           nt = lorenzParam['nt_RK']
           
           print dt
           #Get values
           x_Eu,y_Eu,z_Eu,d_Eu= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
           x_r,y_r,z_r,d_r = LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
     
           x_hr,y_hr,z_hr,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
           
           #Cut the higher res run so it is comparable to the lower
           RK_mark = int(RK_factor)
           d_RKA =  d_RKh[1::RK_mark]
           d_RKA = np.array([0.0] + np.ndarray.tolist(d_RKA))
    
           # calculate the error a
           dError = abs(d_RKA[-1] -d_Eu[-1])
           d_RKA = np.array(d_RKA)
           d_err_dts[index].append(dError)

   

   d_err_dts = [np.average(d_err_dts[i]) for i in range(len(d_err_dts))]

   plt.loglog(dt_range,d_err_dts, linestyle = '--', marker ='x', color = 'firebrick')

   xlog = np.log(dt_range)
   #y1log=np.log(l2_dts)
   y2log = np.log(d_err_dts)
   #n1, intercept = np.polyfit(xlog,y1log,1)
   n2, intercept = np.polyfit(xlog,y2log,1)
   print n2
   
   
orderOfAccuracy()   