# -*- coding: utf-8 -*-
"""
Main to exeute the solving of the lorenz equations using a 4th order runge kutta scheme 
Adaptve time-stepping so far has been in constant multiplication

-- need to implement the constant addition time step growth 

"""
import matplotlib.pyplot as plt
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
   #time axes for variable timestep
   var_time = lorenzParam['t']
   time = lorenzParam['tconst']
   nt - lorenzParam['nt']
   


   
   #Will be the factor to divide the timestep by for higher res
   RK_factor = 100.
   highres_dt = [dt/RK_factor]*int(nt*RK_factor)
   highres_t = np.cumsum(highres_dt)
   
   #RK_dt= dt/RK_factor
   
   
   #Initialise list to store the total-distance-travelled errors
   #for each set of initial conditions and one for the averages
   d_errorRK =[]
   d_errorRKv=[]
   av_d_errorsRK = [] 
   av_d_errorsRKv =[]
   d_errorEu =[]
   d_errorEuv =[]
   av_d_errorsEu =[]
   av_d_errorsEuv =[]


   
   
   #loop over the number of initial conditions to use
   for i in range(5):
       
       #create random intial conditions between -15 and 15
       x_o = -15 + 30 * np.random.random()
       y_o = -15 + 30 * np.random.random()
       z_o = -15 + 30 * np.random.random()
   
   
       #Get values
       d= LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt)
       print 'fixed run'
       d_RKh = LorenzAttractor.RungeKutta(lorenzParam,RK_factor,x_o,y_o,z_o,dt)
       print 'high res run'
       d_v = LorenzAttractor.RungeKutta(lorenzParam,1.,x_o,y_o,z_o,dt_variable)
       print 'variable run'
       d_eu = LorenzAttractor.Euler(lorenzParam,1.,x_o,y_o,z_o,dt)


   
       #Cut the higher res run so it is comparable to the lower
       RK_mark = int(RK_factor)
       d_RKA =  d_RKh[1::RK_mark]
       d_RKA = np.array([0.0] + np.ndarray.tolist(d_RKA))
   

      
       #find the nearest indicies to compare high res and variable time-stepping 
       if i ==0:
           print 'getting indicies'
           indicies = LorenzAttractor.find_nearest_index(highres_t, var_time)
           print 'got indicies'
       
       
       
       #Get the difference between the high res and the variable solutions
       d_diff_RKv =[abs(d_RKh[indicies[g]] - d_v[g]) for g in range(len(indicies))]
       d_diff_Euv = [abs(d_RKh[indicies[g]] - d_v[g]) for g in range(len(indicies))]
           

      
   

   
       #pos_diff_CN = np.sqrt((x_A-x)**2 +(y_A-y)**2 +(z_A-z)**2)
       #pos_diff_RK = np.sqrt((x_RKA-x)**2 +(y_RKA-y)**2 +(z_RKA-z)**2)
       d_diff_RK = abs(d_RKA - d)
       d_diff_eu = abs(d_RKA - d_eu)

       
       d_errorRK.append(d_diff_RK)
       d_errorEu.append(d_diff_eu)
       d_errorRKv.append(d_diff_RKv)
       d_errorEuv.append(d_diff_Euv)
       print i
   print 'doing calcs' 
   '''
   
   av_d_errorsRK  = np.mean(np.array(d_errorRK),axis=0)
   av_d_errorsEu  = np.mean(np.array(d_errorEu),axis =0)
   av_d_errorsRKv = np.mean(np.array(d_errorRKv),axis =0)
   av_d_errorsEuv = np.mean(np.array(d_errorEuv),axis =0)
   
   
   '''
   for a in range(0,len(d_errorRK[0])):
       
       
       #want to average the errors at each step ie the first element of each sublist etc
       #get ath element of each sublist
       errors_at_aRK = [err[a] for err in d_errorRK]
       errors_at_aEu = [er[a] for er in d_errorEu]
       av_d_errorsRK.append(np.average(errors_at_aRK))
       av_d_errorsEu.append(np.average(errors_at_aEu))
   
     
   for a in range(0,len(d_errorRKv[0])): 
       errors_at_aVar =[err[a] for err in d_errorRKv]
       av_d_errorsRKv.append(np.average(errors_at_aVar))
    
     
       
   
   


   LorenzPlot.errorCompar(av_d_errorsRK,av_d_errorsRKv,\
                           lorenzParam,'RK4',time)
   #LorenzPlot.errorCompar(av_d_errorsEu,av_d_errorsEuv,\
                           #lorenzParam,'EulerF',var_time)

 


   
   #LorenzPlot.lorentzPlotting(x[1:],y[1:],z[1:],lorenzParam)
   #LorenzPlot.lorentzPlotting(x_eu[1:],y_eu[1:],z_eu[1:],lorenzParam)
   '''
   plt.clf()
   plt.plot(range(len(d_v)),d_v,'k')
   plt.plot(range(len(d_RKA)),d_RKA,'g')
   plt.plot(range(len(d)),d,'r')
   '''
   
main()

