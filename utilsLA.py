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
   

   

   dt_range = [0.0001,0.0003,0.0005,0.001,0.003,0.005]
   d_err_dts = []*len(dt_range)
   l2_dts = []*len(dt_range)
   d_av_err = [0]*len(dt_range)
   
   #loop over the number of initial conditions to use

       
  # for i in range(0,1):
   
   #create random intial conditions between -15 and 15
   x_o = -15 + 30 * np.random.random()
   y_o = -15 + 30 * np.random.random()
   z_o = -15 + 30 * np.random.random()
   '''
   x_o =ICs[i]
   y_o =ICs[(i+3)%len(ICs)]
   z_o =ICs[(2*i)%len(ICs)]
   '''
   
   
   for dt in dt_range:
       lorenzParam['nt'] = int(0.5/dt)
       nt = lorenzParam['nt']
       print dt
       #Get values
       d= LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)
       d_r = LorenzAttractor.RungeKutta(lorenzParam,1,x_o,y_o,z_o,dt)
 
       d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)

       #Cut the higher res run so it is comparable to the lower
       RK_mark = int(RK_factor)
       d_RKA =  d_RKh[1::RK_mark]
       d_RKA = np.array([0.0] + np.ndarray.tolist(d_RKA))

       # calculate the error a
       dError = abs(d_RKA[-1] -d[-1])
       d_RKA = np.array(d_RKA)
       l2 = np.sqrt((dError**2)/(d_RKA[-1]**2))
       l2_dts.append(l2)
       d_err_dts.append(dError)
       t = [dt]*nt
       time = np.cumsum(t)
       
       plt.plot(time,abs(d_RKA[1:] -d[1:]),'r')
       plt.plot(time,(d_RKA[1:]-d_r[1:]),'k')

           
   '''       
   print d_err_dts[0][0]
   for x in range(len(dt_range)): 
       d_av_err[x] = [np.average(d_err_dts[x][a]) for a in range(len(d_err_dts[x]))]
   l2_av_err = [np.average(l2_dts[a]) for a in range(len(l2_dts))]
   '''
   #plt.loglog(dt_range,d_err_dts, linestyle = '--', marker ='x', color = 'g')
   #plt.loglog(dt_range,l2_dts, linestyle = '--', marker ='x', color = 'g')
   xlog = np.log(dt_range)
   #y1log=np.log(l2_dts)
   y2log = np.log(d_err_dts)
   #n1, intercept = np.polyfit(xlog,y1log,1)
   n2, intercept = np.polyfit(xlog,y2log,1)
   print n2
   
   
orderOfAccuracy()   