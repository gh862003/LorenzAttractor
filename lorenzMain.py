# -*- coding: utf-8 -*-
"""
Main to exeute the solving of the lorenz equations using a 4th order runge kutta scheme 
and CN
"""

import LorenzAttractor
import LorenzPlot
import numpy as np

def main():
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 10.
   
   #Get values
   x,y,z = LorenzAttractor.RungeKutta(lorenzParam,1.)
   x_CN,y_CN,z_CN = LorenzAttractor.CN(lorenzParam)
   x_RKh,y_RKh,z_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor)
   

   #take evry tenth after the first
   x_A = x_CN[1::10]
   y_A = y_CN[1::10]
   z_A = z_CN[1::10]
   r_A = np.sqrt(x_A**2 +y_A**2 +z_A**2)
   
   #Have also done a run of RK at a factor RK_factor of the timestep
   RK_mark = int(RK_factor)
   x_RKA = x_RKh[1::RK_mark]
   y_RKA = y_RKh[1::RK_mark]
   z_RKA = z_RKh[1::RK_mark]
   r_RK = np.sqrt(x_RKA**2 +y_RKA**2 +z_RKA**2)
   
   #remove IC for easier comparison for now
   x = x[1:]
   y = y[1:]
   z = z[1:]
   r = np.sqrt(x**2+y**2+z**2)
   
   diff_CN = abs(r_A -r)
   diff_RK = abs(r_RK - r)

   LorenzPlot.errorPlot(diff_CN,lorenzParam,'CN')
   LorenzPlot.errorPlot(diff_RK,lorenzParam,'RK'+str(RK_factor))
 


   
   LorenzPlot.lorentzPlotting(x,y,z,lorenzParam)
   
main()