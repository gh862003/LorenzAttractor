# -*- coding: utf-8 -*-
"""
Main to exeute the solving of the lorenz equations using a 4th order runge kutta scheme 


"""

import LorenzAttractor
import LorenzPlot
import numpy as np

def main():
   #Set up dictionary of input parameters from lorenzParams file
   lorenzParam = {}
   execfile("lorenzParams.py", lorenzParam)
   #access time-step size
   dt = float(lorenzParam['dt'])
   #access list of variable step sizes
   dt_variable = lorenzParam['variable_dt']
   
   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 10.
   
   RK_dt= dt/RK_factor
   #st_varible = 
   
   #Initialise list to store the total-distance-travelled errors
   #for each set of initial conditions and one for the averages
   distance_errors =[]
   var_dist_errors=[]
   av_distance_errors = [] 
   var_av_dist_erros =[]
   
   
   #loop over the number of initial conditions to use
   for i in range(0,30):
       
       #create random intial conditions between -15 and 15
       x_o =-15 + 30 * np.random.random()
       y_o =-15 + 30 * np.random.random()
       z_o =-15 + 30 * np.random.random()
   
   
   
   
       #Get values
       x,y,z,d= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
       #x_CN,y_CN, z_CN = LorenzAttractor.CN(lorenzParam)
       x_RKh,y_RKh,z_RKh,d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
       x_v,y_v,z_v,d_v = LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt_variable)
   
       #take evry tenth after the first
       #x_A = x_CN[1::10]
       #y_A = y_CN[1::10]
       #z_A = z_CN[1::10]

   
       #Cut the higher res run so it is comparable to the lower
       RK_mark = int(RK_factor)
       x_RKA = x_RKh[1::RK_mark]
       y_RKA = y_RKh[1::RK_mark]
       z_RKA = z_RKh[1::RK_mark]
       d_RKA = d_RKh[1::RK_mark]

      
       #remove IC for easier comparison for now
       x = x[1:]
       y = y[1:]
       z = z[1:]
       d = d[1:]
       d_v = d_v[1:]
   

   
       #pos_diff_CN = np.sqrt((x_A-x)**2 +(y_A-y)**2 +(z_A-z)**2)
       #pos_diff_RK = np.sqrt((x_RKA-x)**2 +(y_RKA-y)**2 +(z_RKA-z)**2)
       dist_diff_RK = abs(d_RKA - d)
       dist_diff_Variable = abs(d_RKA-d_v)
       
       distance_errors.append(dist_diff_RK)
       var_dist_errors.append(dist_diff_Variable)

   
   for a in range(0,len(distance_errors[0])):
       
       #want to average the errors at each step ie the first element of each sublist etc
       #get ath element of each sublist
       errors_at_aRK = [err[a] for err in distance_errors]
       errors_at_aVar =[err[a] for err in var_dist_errors]
       av_distance_errors.append(np.average(errors_at_aRK))
       var_av_dist_erros.append(np.average(errors_at_aVar))
       

   

   #LorenzPlot.errorPlot(pos_diff_CN,lorenzParam,'CN')
   #LorenzPlot.errorPlot(pos_diff_RK,lorenzParam,'RK'+str(RK_factor))
   LorenzPlot.errorPlot(av_distance_errors,lorenzParam,str(RK_factor),'const dt')
   LorenzPlot.errorPlot(var_av_dist_erros,lorenzParam,1,'varible dt')
 


   
   LorenzPlot.lorentzPlotting(x,y,z,lorenzParam)

   
main()