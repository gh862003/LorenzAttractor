# -*- coding: utf-8 -*-
"""
Order of accuracy
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
   

   #use non-cumulative error measure

   dt_range = [0.00075,0.0015,0.003]
   #dt_range = [0.00625,0.0125,0.025]
   d_err_dts = [[],[],[]]
   str_coords = [line.rstrip('\n') for line in open('IC_selected.txt')]
   coords = []

   for str_coord in str_coords:
       float_coords = [ float(x) for x in str_coord.split(',')]
       coords.append(float_coords)
   
   #loop over the number of initial conditions to use

       
   for i in range(0,2):
   
       #create random intial conditions between -15 and 15
       x_o = coords[i][0]
       y_o = coords[i][1]
       z_o = coords[i][2]
       
       for index,dt in enumerate(dt_range):
           #make sure they run for the same duration
           lorenzParam['nt_RK'] = int(round(1.0/dt))
           lorenzParam['nt_FB'] = int(round(1.0/dt))
           lorenzParam['nt_Eu'] = int(round(1.0/dt))
           
           print dt
           #Get values
           x_Eu,y_Eu,z_Eu,d_Eu= LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)
           x_r,y_r,z_r,d_r = LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
           x_FB,y_FB,z_FB,d_FB = LorenzAttractor.ForwardBack(x_o,y_o,z_o, lorenzParam,1.,dt)
           x_hr,y_hr,z_hr,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
           
           #Cut the higher res run so it is comparable to the lower
           RK_mark = int(RK_factor)
           d_RKA =  d_RKh[1::RK_mark]
           d_RKA = np.array([d_RKh[0]] + np.ndarray.tolist(d_RKA))
      
           
           # calculate the error a
           d_diff_RK = abs(d_RKA - d_Eu)
           d_diff_RK = np.cumsum(np.array(d_diff_RK)**2)
           
           d_diff_RK = [np.sqrt(r/ind) for ind,r in enumerate(d_diff_RK)]
           d_diff_RK[0]=0 #Ics are the same
           dError = d_diff_RK[-1]
           print dError
           d_RKA = np.array(d_RKA)

           d_err_dts[index].append(dError)
           '''
           #try different error measure
           error = abs(d_r[-1] - d_RKA[-1])
           d_err_dts[index].append(error)
           '''
   
   print d_err_dts
   d_err_dts  = np.mean(np.array(d_err_dts),axis = 1)
   print len(d_err_dts)
   print d_err_dts

   plt.loglog(dt_range,d_err_dts, linestyle = '--', marker ='x', color = 'firebrick')

   xlog = np.log(dt_range)
   y2log = np.log(d_err_dts)
   n2, intercept = np.polyfit(xlog,y2log,1)
   print n2
   
   
orderOfAccuracy()   